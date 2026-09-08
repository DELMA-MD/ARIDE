import os
import re
import time
import requests
import streamlit as st

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
MYVARIANT_URL = "https://myvariant.info/v1/variant/{hgvs}"

# IndiGen has no public API - point this at a local, bgzip+tabix-indexed copy
# of the IndiGen VCF (download from https://clingen.igib.res.in/indigen/).
# Leave as-is if you don't have it yet; the lookup will just report "unavailable".
INDIGEN_VCF_PATH = os.environ.get("INDIGEN_VCF_PATH", "data/indigen.vcf.gz")

# ACMG population-frequency thresholds (PM2/BS1/BA1)
BA1_FAF = 0.05   # popmax AF above this -> Benign Standalone (too common for a rare Mendelian disease)
BS1_FAF = 0.01   # popmax AF above this -> Benign Strong

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


def fetch_gnomad(hgvs_genomic: str) -> dict:
    """Query myvariant.info for current gnomAD exome + genome popmax / South Asian AF."""
    params = {
        "assembly": "hg38",
        "fields": "gnomad_exome.af.af_popmax,gnomad_exome.af.af_sas,"
                  "gnomad_genome.af.af_popmax,gnomad_genome.af.af_sas,"
                  "exac.af,dbnsfp.1000gp3.af,dbsnp.rsid",
    }
    r = requests.get(MYVARIANT_URL.format(hgvs=hgvs_genomic), params=params, timeout=20)
    r.raise_for_status()
    return r.json()


def popmax_af(data: dict) -> float:
    """Higher of gnomAD exome/genome popmax AF. 0.0 means not observed (novel)."""
    afs = []
    for key in ("gnomad_exome", "gnomad_genome"):
        block = data.get(key, {}).get("af", {})
        val = block.get("af_popmax")
        if isinstance(val, (int, float)):
            afs.append(val)
    return max(afs) if afs else 0.0


def sas_af(data: dict) -> float:
    """Higher of gnomAD exome/genome South Asian AF. 0.0 means not observed."""
    afs = []
    for key in ("gnomad_exome", "gnomad_genome"):
        block = data.get(key, {}).get("af", {})
        val = block.get("af_sas")
        if isinstance(val, (int, float)):
            afs.append(val)
    return max(afs) if afs else 0.0


def exac_af(data: dict):
    """ExAC overall AF, if present (ExAC is deprecated/merged into gnomAD v2 - kept for
    completeness since the original clinical report cited it separately)."""
    val = data.get("exac", {}).get("af")
    return float(val) if isinstance(val, (int, float)) else None


def thousand_genomes_af(data: dict):
    """1000 Genomes Phase 3 overall AF, via dbNSFP."""
    val = data.get("dbnsfp", {}).get("1000gp3", {}).get("af")
    if isinstance(val, list):
        val = max(v for v in val if isinstance(v, (int, float))) if val else None
    return float(val) if isinstance(val, (int, float)) else None


def dbsnp_rsid(data: dict):
    """dbSNP rsID, if the variant has one on record (dbSNP is primarily an ID registry,
    not a curated frequency source, so we report presence/rsID rather than an AF)."""
    return data.get("dbsnp", {}).get("rsid")


def fetch_indigen(chrom: str, pos: str, ref: str, alt: str):
    """Look up allele frequency in a local, tabix-indexed IndiGen VCF.
    Returns None if the file isn't available (IndiGen has no public API)."""
    if not os.path.exists(INDIGEN_VCF_PATH):
        return None, "IndiGen VCF not found locally - download + tabix-index it to enable this check."
    try:
        import pysam  # only needed if IndiGen lookups are actually used
    except ImportError:
        return None, "pysam not installed (pip install pysam) - required to query the local IndiGen VCF."

    try:
        vcf = pysam.VariantFile(INDIGEN_VCF_PATH)
        chrom_query = chrom if chrom.startswith("chr") else f"chr{chrom}"
        for rec in vcf.fetch(chrom_query, int(pos) - 1, int(pos)):
            if rec.ref == ref and alt in rec.alts:
                af = rec.info.get("AF")
                af = af[0] if isinstance(af, (list, tuple)) else af
                return float(af) if af is not None else 0.0, None
        return 0.0, None  # position covered, variant not seen -> effectively absent/novel in IndiGen
    except Exception as e:
        return None, f"IndiGen lookup failed: {e}"


def classify_by_population(faf: float) -> dict:
    """Apply the ACMG/AMP population-frequency rules (BA1 / BS1 / PM2) to a popmax AF.
    Returns the ACMG code triggered and the resulting classification call."""
    if faf > BA1_FAF:
        return {
            "code": "BA1",
            "classification": "Benign",
        }
    if faf > BS1_FAF:
        return {
            "code": "BS1",
            "classification": "Likely Benign",
        }
    if faf > 0:
        return {
            "code": "PM2 (weak)",
            "classification": "Uncertain Significance",
        }
    return {
        "code": "PM2",
        "classification": "Uncertain Significance",
    }


def classify_across_databases(variant_record: dict) -> dict:
    """Query gnomAD, 1000 Genomes, ExAC, dbSNP (all via myvariant.info) plus IndiGen
    (local file), classify each individually against ACMG population-frequency rules,
    and combine into one overall call.

    ACMG combining logic: BA1 from ANY single database is stand-alone sufficient for
    Benign. Otherwise BS1 from any database -> Likely Benign. Otherwise, agreement
    across multiple databases that the variant is rare/absent still only supports
    PM2 (Moderate) - ACMG does not let repeated "it's rare" observations stack into
    something stronger than PM2 alone.
    """
    hgvs_genomic = variant_record["hgvs_genomic"]
    chrom = variant_record["chromosome"]
    pos = variant_record["position"]
    ref = variant_record["ref_allele"]
    alt = variant_record["alt_allele"]

    per_db = {}

    # --- gnomAD, 1000G, ExAC, dbSNP: one myvariant.info call covers all four ---
    try:
        data = fetch_gnomad(hgvs_genomic)
        gnomad_popmax = popmax_af(data)
        gnomad_sas = sas_af(data)
        per_db["gnomAD"] = {"af": gnomad_popmax, **classify_by_population(gnomad_popmax)}
        per_db["gnomAD (South Asian)"] = {"af": gnomad_sas, **classify_by_population(gnomad_sas)}

        exac = exac_af(data)
        per_db["ExAC"] = (
            {"af": exac, **classify_by_population(exac)} if exac is not None
            else {"af": None, "code": "n/a", "classification": "Not in ExAC (superseded by gnomAD - kept for reference only)"}
        )

        kg = thousand_genomes_af(data)
        per_db["1000 Genomes"] = (
            {"af": kg, **classify_by_population(kg)} if kg is not None
            else {"af": None, "code": "n/a", "classification": "Not observed in 1000 Genomes"}
        )

        rsid = dbsnp_rsid(data)
        per_db["dbSNP"] = {"af": None, "code": "n/a", "classification": f"rsID: {rsid}" if rsid else "No dbSNP entry"}

    except Exception as e:
        per_db["myvariant.info (gnomAD/1000G/ExAC/dbSNP)"] = {"af": None, "code": "error", "classification": str(e)}

    # --- IndiGen: separate local lookup ---
    indigen_af, indigen_msg = fetch_indigen(chrom, pos, ref, alt)
    if indigen_af is not None:
        per_db["IndiGen"] = {"af": indigen_af, **classify_by_population(indigen_af)}
    else:
        per_db["IndiGen"] = {"af": None, "code": "unavailable", "classification": indigen_msg}

    # --- Combine into one overall call ---
    frequency_based = [v for v in per_db.values() if isinstance(v.get("af"), (int, float))]
    if any(v["code"] == "BA1" for v in frequency_based):
        overall = "Benign (BA1 triggered by at least one population database)"
    elif any(v["code"] == "BS1" for v in frequency_based):
        overall = "Likely Benign (BS1 triggered by at least one population database)"
    else:
        overall = "Uncertain Significance (PM2 only - rare/absent across all queried population databases; insufficient alone for Likely Pathogenic)"

    return {"per_database": per_db, "overall": overall}


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

            # --- Query all population databases and apply ACMG classification per database ---
            with st.spinner("Checking gnomAD, 1000 Genomes, ExAC, dbSNP, IndiGen..."):
                try:
                    result = classify_across_databases(variant_record)
                    db_error = None
                except Exception as e:
                    result, db_error = None, e

            if db_error is not None:
                st.error(f"Database lookup failed: {db_error}")
            else:
                st.divider()
                st.subheader("Population database check")

                table_rows = []
                for db_name, info in result["per_database"].items():
                    af_display = f"{info['af']:.4%}" if isinstance(info.get("af"), (int, float)) else "—"
                    table_rows.append({
                        "Database": db_name,
                        "Allele frequency": af_display,
                        "ACMG code": info.get("code", "n/a"),
                        "Note / classification": info.get("classification", ""),
                    })
                st.table(table_rows)

                overall = result["overall"]
                if "Benign" in overall and "Uncertain" not in overall:
                    st.success(f"**Overall classification: {overall}**")
                else:
                    st.warning(f"**Overall classification: {overall}**")

                variant_record["population_check"] = result

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
