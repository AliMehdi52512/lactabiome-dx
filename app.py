import streamlit as st
import pandas as pd
import numpy as np

# Page configuration
st.set_page_config(
    page_title="LactaBiome-Dx | let-7g-5p qPCR Decision Support",
    page_icon="🥛",
    layout="wide"
)

# Header
st.title("🥛 LactaBiome-Dx: RT-qPCR Decision Support Platform")
st.markdown("""
**Diagnostic stratification of low milk supply risk via milk-derived microRNA quantification.**  
*Translational reference model based on maternal lactation biomarkers (Hicks, Kelleher et al., 2023).*
""")
st.divider()

# Sidebar - Configuration and reference
st.sidebar.header("🔬 Assay Setup & Controls")
ref_gene = st.sidebar.selectbox(
    "Endogenous Normalizer Control",
    ["U6 snRNA", "miR-16-5p", "miR-30b-5p"],
    index=0
)
sample_source = st.sidebar.selectbox(
    "Milk Fraction Analyzed",
    ["Skim Milk / Milk Exosomes (EVs)", "Whole Milk Colostrum", "Milk Fat Globule Membrane (MFGM)"]
)

st.sidebar.markdown("---")
st.sidebar.subheader("Baseline Reference Values")
calibrator_dct = st.sidebar.number_input(
    "Calibrator Baseline ΔCt (Healthy Control)",
    min_value=-5.0,
    max_value=10.0,
    value=3.20,
    step=0.10,
    help="Default 3.20 reflects baseline ΔCt (let-7g-5p - U6) established in normolactating cohorts."
)

st.sidebar.markdown("---")
st.sidebar.info("Developed for computational and translational screening workflows.")

# Main layout
tab1, tab2 = st.tabs(["📊 Individual Sample Analysis", "📑 Batch Cohort Evaluation"])

# TAB 1: Individual Sample Analysis
with tab1:
    col1, col2 = st.columns([1, 1], gap="large")

    with col1:
        st.subheader("1. Enter PCR Cycle Thresholds (Ct)")
        sample_id = st.text_input("Patient / Specimen ID", value="PATIENT-MLK-001")
        
        target_ct = st.number_input(
            "let-7g-5p Target Ct",
            min_value=10.0,
            max_value=40.0,
            value=23.4,
            step=0.1,
            help="Cycle threshold for let-7g-5p (lower Ct means higher abundance)."
        )
        
        ref_ct = st.number_input(
            f"{ref_gene} Reference Ct",
            min_value=10.0,
            max_value=40.0,
            value=26.8,
            step=0.1,
            help="Cycle threshold for the selected reference control."
        )

        # Mathematical calculations (Livak & Schmittgen 2^-ddCt method)
        delta_ct = target_ct - ref_ct
        ddct = delta_ct - calibrator_dct
        fold_change = 2.0 ** (-ddct)

    with col2:
        st.subheader("2. Biomarker & Molecular Risk Assessment")

        # Metric displays
        m1, m2 = st.columns(2)
        m1.metric("Sample ΔCt", f"{delta_ct:.2f}")
        m2.metric("Relative Fold-Change (2^-ΔΔCt)", f"{fold_change:.2f}x")

        # Clinical stratification logic
        if fold_change >= 2.50:
            st.error("🚨 HIGH RISK: Lactation Failure / Secretory Insufficiency")
            st.markdown(f"""
            - **Status:** Marked elevation of `let-7g-5p` detected ({fold_change:.2f}-fold relative to healthy baseline).
            - **Molecular Mechanism:** High `let-7g-5p` strongly represses **PRLR** (Prolactin Receptor) translation in mammary epithelial cells.
            - **Translational Action:** Candidate for targeted galactagogue therapy (e.g., botanical polyphenol screening) or hormone pathway support.
            """)
        elif 1.40 <= fold_change < 2.50:
            st.warning("⚠️️ MODERATE RISK: Delayed Secretory Activation (Lactogenesis II)")
            st.markdown(f"""
            - **Status:** Borderline elevation of `let-7g-5p` ({fold_change:.2f}-fold).
            - **Molecular Mechanism:** Partial suppression of prolactin receptor downstream cascades (Jak2/Stat5).
            - **Translational Action:** Re-sample milk at 48–72 hours postpartum and observe feeding efficiency.
            """)
        else:
            st.success("✅ LOW RISK: Normal Lactogenesis Profile")
            st.markdown(f"""
            - **Status:** Baseline physiological expression of `let-7g-5p` ({fold_change:.2f}-fold).
            - **Molecular Mechanism:** Normal prolactin receptor signaling capacity in alveolar epithelial lactocytes.
            - **Translational Action:** Standard postpartum lactation support.
            """)

# TAB 2: Batch Cohort Evaluation
with tab2:
    st.subheader("Simulate Cohort Testing")
    st.write("Demonstration data evaluating a cohort of 5 post-delivery specimens:")

    cohort_data = pd.DataFrame({
        "Sample_ID": ["MLK-001", "MLK-002", "MLK-003", "MLK-004", "MLK-005"],
        "Target_Ct (let-7g-5p)": [22.8, 26.5, 21.0, 25.1, 28.0],
        "Reference_Ct (U6)": [26.0, 26.2, 25.8, 25.4, 26.0],
    })
    
    # Batch calculate
    cohort_data["ΔCt"] = cohort_data["Target_Ct (let-7g-5p)"] - cohort_data["Reference_Ct (U6)"]
    cohort_data["ΔΔCt"] = cohort_data["ΔCt"] - calibrator_dct
    cohort_data["Fold_Change"] = 2.0 ** (-cohort_data["ΔΔCt"])
    cohort_data["Predicted_Status"] = np.where(
        cohort_data["Fold_Change"] >= 2.50, "High Risk",
        np.where(cohort_data["Fold_Change"] >= 1.40, "Moderate Risk", "Normal")
    )

    st.dataframe(cohort_data, use_container_width=True)

st.markdown("---")
st.caption("LactaBiome-Dx | Open-source computational diagnostics tool.")
