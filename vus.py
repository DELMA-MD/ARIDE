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
    match = re.search(r":(c\.[^\s(]+)", title)
    return match.group(1).strip() if match else ""


def query_clinvar(gene: str, variant: str, retries: int = 3):
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

        if extract_c_notation(title).lower() != variant_clean:
            continue

        germline = rec.get("germline_classification", {}) or {}
        last_evaluated_raw = germline.get("last_evaluated", "")
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
    c = classification.strip().lower()
    is_settled = ("pathogenic" in c or "benign" in c) and "conflicting" not in c and "uncertain" not in c
    return not is_settled


st.set_page_config(page_title="VUS Reclassification", page_icon="🧬", layout="centered")
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
                    st.warning("Uncertain significance in ClinVar. Proceed to Reclassification.")
            else:
                label = "Pathogenic" if "pathogenic" in primary else "Benign"
                st.success(
                    f"Already classified as {label} based on current ClinVar. "
                    f"Reclassification is not needed."
                )

            st.session_state["start"] = {
                "gene": gene,
                "variant": variant,
                "status": primary,
                "records": results,
            }

if "start" in st.session_state:
    st.divider()
    status = st.session_state["start"]["status"]
    reclass_enabled = needs_reclassification(status)

    if st.button("Reclassification", disabled=not reclass_enabled):
        st.info("Under processing...")