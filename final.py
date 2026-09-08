import re
import time
import requests
import streamlit as st

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

CLASS_COLORS = {
    "pathogenic": "#b23b3b",
    "likely pathogenic": "#c9645a",
    "pathogenic/likely pathogenic": "#b23b3b",
    "uncertain significance": "#c98a1a",
    "benign": "#2a7f4f",
    "likely benign": "#4f9f6f",
    "benign/likely benign": "#2a7f4f",
    "conflicting classifications of pathogenicity": "#8a5fb2",
}


def extract_c_notation(title: str) -> str:
    """Pull the c. HGVS notation out of a ClinVar title like
    'NM_001048174.2(MUTYH):c.1142G>A (p.Trp381Ter)' -> 'c.1142G>A'."""
    match = re.search(r":(c\.[^\s(]+)", title)
    return match.group(1).strip() if match else ""


def query_clinvar(gene: str, variant: str, retries: int = 3):
    """Search ClinVar for a gene + HGVS variant, return classification records
    filtered down to an EXACT match on the c. notation the user typed
    (ClinVar's own search is fuzzy and returns other variants in the same gene)."""
    term = f'{gene}[gene] AND "{variant}"[Variant name]'
    params = {"db": "clinvar", "term": term, "retmode": "json", "retmax": 20}

    ids = []
    last_err = None
    for attempt in range(retries):
        try:
            r = requests.get(ESEARCH_URL, params=params, timeout=15)
            r.raise_for_status()
            ids = r.json().get("esearchresult", {}).get("idlist", [])
            last_err = None
            break
        except Exception as e:
            last_err = e
            time.sleep(1)
    if last_err is not None:
        raise last_err

    if not ids:
        params["term"] = f"{gene}[gene] AND {variant}"
        r = requests.get(ESEARCH_URL, params=params, timeout=15)
        r.raise_for_status()
        ids = r.json().get("esearchresult", {}).get("idlist", [])

    if not ids:
        return []

    sum_params = {"db": "clinvar", "id": ",".join(ids), "retmode": "json"}
    r = requests.get(ESUMMARY_URL, params=sum_params, timeout=15)
    r.raise_for_status()
    result = r.json().get("result", {})

    variant_clean = variant.strip().lower()
    records = []
    for uid in ids:
        rec = result.get(uid)
        if not rec:
            continue
        title = rec.get("title", "")

        # Exact-match filter: only keep this record if its own c. notation
        # matches what the user typed, character-for-character.
        if extract_c_notation(title).lower() != variant_clean:
            continue

        germline = rec.get("germline_classification", {}) or {}
        last_evaluated_raw = germline.get("last_evaluated", "")
        # Keep only the date part, drop the time (e.g. "2026/01/26 00:00" -> "2026/01/26")
        last_evaluated_date = last_evaluated_raw.split(" ")[0] if last_evaluated_raw else ""
        records.append(
            {
                "accession": rec.get("accession", uid),
                "title": title,
                "classification": germline.get("description", "Not provided"),
                "last_evaluated": last_evaluated_date,
                "link": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{uid}/",
            }
        )
    return records


def badge(classification: str) -> str:
    color = CLASS_COLORS.get(classification.strip().lower(), "#5a6675")
    return (
        f'<span style="background:{color};color:white;padding:2px 10px;'
        f'border-radius:12px;font-size:0.85em;font-weight:600">{classification}</span>'
    )


def needs_reclassification(classification: str) -> bool:
    """Only Pathogenic/Likely Pathogenic and Benign/Likely Benign are settled.
    Everything else (VUS, conflicting, not provided, not found) still needs work."""
    c = classification.strip().lower()
    is_settled = ("pathogenic" in c or "benign" in c) and "conflicting" not in c and "uncertain" not in c
    return not is_settled


st.set_page_config(page_title="VUS Reclassification", page_icon="🧬", layout="centered")

# Which "slide" is currently showing: "clinvar" (the check) or "reclass" (the intake form)
if "view" not in st.session_state:
    st.session_state["view"] = "clinvar"


# =====================================================================
# SLIDE 1 — ClinVar check
# =====================================================================
def render_clinvar_slide():
    st.title("ClinVar Verification")
    st.caption("Check the current ClinVar classification")

    with st.form("intake_form"):
        c1, c2 = st.columns(2)
        with c1:
            gene = st.text_input("Gene", placeholder="e.g. EYS")
        with c2:
            variant = st.text_input("HGVS variant (c. notation)", placeholder="e.g. c.9270_9271insA")
        submitted = st.form_submit_button("Check ClinVar")

    if submitted:
        if not gene.strip() or not variant.strip():
            st.warning("Enter both gene symbol and variant.")
        else:
            with st.spinner(f"Querying ClinVar for {gene} {variant}..."):
                try:
                    results = query_clinvar(gene.strip(), variant.strip())
                    error = None
                except Exception as e:
                    results, error = None, e

            if error is not None:
                st.error(f"ClinVar lookup failed: {error}")
            elif not results:
                st.info(
                    f"No ClinVar record found for {gene} {variant}. "
                    "It may be novel or not yet submitted, proceed to Reclassification"
                )
                st.session_state["start"] = {"gene": gene, "variant": variant, "status": "not_found"}
                st.session_state.pop("reclass_variants", None)
            else:
                st.success(f"Found {len(results)} ClinVar record(s) for {gene} {variant}")
                for rec in results:
                    st.markdown(
                        f"**{rec['title']}**  \n"
                        f"{badge(rec['classification'])} &nbsp;·&nbsp; "
                        f"Last evaluated: {rec['last_evaluated'] or '—'}  \n"
                        f"[View in ClinVar]({rec['link']}) &nbsp;·&nbsp; Accession: `{rec['accession']}`",
                        unsafe_allow_html=True,
                    )
                    st.divider()

                primary = results[0]["classification"].strip().lower()

                if needs_reclassification(primary):
                    if "conflicting" in primary:
                        st.warning("Conflicting classifications in ClinVar. Proceed to Reclassification.")
                    else:
                        st.warning("Uncertain significance / not provided. Proceed to Reclassification.")
                else:
                    label = "Pathogenic" if "pathogenic" in primary else "Benign"
                    st.success(
                        f"Already classified as {label} based on current ClinVar evidence — "
                        f"reclassification is not needed."
                    )

                st.session_state["start"] = {
                    "gene": gene,
                    "variant": variant,
                    "status": primary,
                    "records": results,
                }
                st.session_state.pop("reclass_variants", None)

    if "start" in st.session_state:
        st.divider()
        status = st.session_state["start"]["status"]
        reclass_enabled = needs_reclassification(status)

        if st.button("Reclassification →", disabled=not reclass_enabled, use_container_width=True):
            st.session_state["view"] = "reclass"
            st.rerun()


# =====================================================================
# SLIDE 2 — Variant intake (reclassification)
# =====================================================================
def render_reclass_slide():
    if st.button("← Back to ClinVar check"):
        st.session_state["view"] = "clinvar"
        st.rerun()

    st.title("Variant Intake")
    st.caption("Enter the variant details required for database/tool lookup")

    # Pre-fill gene + HGVS coding from what was already typed in the ClinVar check,
    # so the user doesn't have to retype what's already known.
    prior = st.session_state.get("start", {})

    with st.form("variant_form"):

        col1, col2 = st.columns(2)
        with col1:
            gene_symbol = st.text_input("Gene symbol", value=prior.get("gene", ""), placeholder="ARL6")
        with col2:
            transcript_id = st.text_input("Transcript ID", placeholder="NM_177976.3")

        col3, col4 = st.columns(2)
        with col3:
            chromosome = st.text_input("Chromosome", placeholder="3")
        with col4:
            position = st.text_input("Position (GRCh38)", placeholder="97784972")

        col5, col6 = st.columns(2)
        with col5:
            ref_allele = st.text_input("Reference allele", placeholder="T")
        with col6:
            alt_allele = st.text_input("Alternate allele", placeholder="C")

        col7, col8 = st.columns(2)
        with col7:
            hgvs_c = st.text_input("HGVS coding (c.)", value=prior.get("variant", ""), placeholder="c.272T>C")
        with col8:
            hgvs_p = st.text_input("HGVS protein (p.)", placeholder="p.Ile91Thr")

        variant_submitted = st.form_submit_button("Submit variant")

    if variant_submitted:
        required = {
            "Gene symbol": gene_symbol,
            "Chromosome": chromosome,
            "Position": position,
            "Reference allele": ref_allele,
            "Alternate allele": alt_allele,
        }
        missing = [name for name, val in required.items() if not val.strip()]

        if missing:
            st.error(f"Missing required field(s): {', '.join(missing)}")
        else:
            variant_record = {
                "gene_symbol": gene_symbol.strip(),
                "transcript_id": transcript_id.strip(),
                "chromosome": chromosome.strip(),
                "position": position.strip(),
                "ref_allele": ref_allele.strip().upper(),
                "alt_allele": alt_allele.strip().upper(),
                "hgvs_c": hgvs_c.strip(),
                "hgvs_p": hgvs_p.strip(),
                "hgvs_genomic": f"chr{chromosome.strip()}:g.{position.strip()}{ref_allele.strip().upper()}>{alt_allele.strip().upper()}",
            }

            st.success("Variant recorded.")
            st.json(variant_record)
            if "reclass_variants" not in st.session_state:
                st.session_state.reclass_variants = []
            st.session_state.reclass_variants.append(variant_record)

    if st.session_state.get("reclass_variants"):
        st.divider()
        st.subheader(f"Variants entered this session ({len(st.session_state.reclass_variants)})")
        st.table(st.session_state.reclass_variants)


# --- Render whichever slide is active. Only one is ever on screen at a time. ---
if st.session_state["view"] == "clinvar":
    render_clinvar_slide()
else:
    render_reclass_slide()