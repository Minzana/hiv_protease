import streamlit as st
import pandas as pd
import sys
import os
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

# Add parent directory so we can import candidates
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from candidates import CANDIDATES

# 1. Security Check
if 'authenticated' not in st.session_state or not st.session_state['authenticated']:
    st.warning("Please login on the main page first.")
    st.stop()

st.set_page_config(page_title="ADMET Analysis", page_icon="💊", layout="wide")

st.title("💊 ADMET Analysis")
st.markdown("### Pharmacokinetic & Toxicity Profiling")
st.markdown("""
**ADMET** stands for **Absorption, Distribution, Metabolism, Excretion, and Toxicity** — the key 
pharmacokinetic properties that determine whether a drug candidate will succeed or fail in clinical trials.

> ~60% of drug candidates fail due to poor ADMET profiles. This analysis uses **RDKit molecular 
> descriptors** combined with established pharmaceutical rules to predict each property.
""")
st.markdown("---")


# 2. ADMET Computation
@st.cache_data
def compute_admet():
    results = []
    for c in CANDIDATES:
        mol = Chem.MolFromSmiles(c['smiles'])
        if mol is None:
            continue

        mw = round(Descriptors.MolWt(mol), 2)
        logp = round(Descriptors.MolLogP(mol), 2)
        tpsa = round(Descriptors.TPSA(mol), 2)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        aromatic_rings = Descriptors.NumAromaticRings(mol)
        heavy_atoms = Descriptors.HeavyAtomCount(mol)

        # -- Rule-based ADMET predictions --

        # Lipinski Rule of Five (drug-likeness / oral absorption)
        lipinski_violations = sum([
            mw > 500,
            logp > 5,
            hbd > 5,
            hba > 10,
        ])
        lipinski_pass = lipinski_violations <= 1

        # Veber's Rules (oral bioavailability)
        veber_pass = tpsa <= 140 and rotatable_bonds <= 10

        # GI Absorption (based on Lipinski + Veber)
        gi_absorption = "High" if (lipinski_pass and veber_pass) else "Low"

        # BBB Penetration (Egan egg model: LogP vs TPSA)
        # BBB permeable: TPSA < 90 and LogP < 3
        bbb_penetration = "Yes" if (tpsa < 90 and logp < 3) else "No"

        # CYP2D6 Inhibition risk (rough heuristic: basic nitrogen + aromatic rings)
        # This is a simplified rule — not a formal QSAR model
        has_basic_nitrogen = 'N' in c['smiles'] and aromatic_rings >= 2
        cyp2d6_risk = "Likely" if (has_basic_nitrogen and logp > 2) else "Unlikely"

        # Drug-likeness score (QED - Quantitative Estimate of Drug-likeness)
        from rdkit.Chem import QED
        qed_score = round(QED.qed(mol), 3)

        # Bioavailability Score (based on Lipinski violations)
        if lipinski_violations == 0:
            bioavailability = "0.85 (Good)"
        elif lipinski_violations == 1:
            bioavailability = "0.56 (Moderate)"
        else:
            bioavailability = "0.17 (Poor)"

        results.append({
            "Candidate": f"#{c['id']}",
            "SMILES": c['smiles'],
            "MW": mw,
            "LogP": logp,
            "TPSA": tpsa,
            "HBD": hbd,
            "HBA": hba,
            "Rot. Bonds": rotatable_bonds,
            "Lipinski": "✅ Pass" if lipinski_pass else "❌ Fail",
            "GI Absorb.": gi_absorption,
            "BBB": bbb_penetration,
            "CYP2D6": cyp2d6_risk,
            "QED": qed_score,
            "Bioavail.": bioavailability,
        })

    return results


if st.button("🧬 Run ADMET Analysis"):
    st.session_state['admet_computed'] = True

if st.session_state.get('admet_computed', False):
    with st.spinner("Computing ADMET properties..."):
        results = compute_admet()

    df = pd.DataFrame(results)

    st.success(f"ADMET Profiling Complete: {len(results)} candidates analyzed.")

    # 3. ADMET Overview Table
    st.subheader("📋 ADMET Profile Overview")
    display_cols = ["Candidate", "MW", "LogP", "TPSA", "HBD", "HBA", "Rot. Bonds",
                    "Lipinski", "GI Absorb.", "BBB", "CYP2D6", "QED", "Bioavail."]
    st.dataframe(
        df[display_cols],
        use_container_width=True,
        hide_index=True,
    )

    # 4. Property Explanations
    st.markdown("---")
    st.subheader("📖 Property Reference Guide")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("#### 🔵 Absorption Properties")
        st.markdown("""
        | Property | Threshold | Meaning |
        |----------|-----------|---------|
        | **MW** | < 500 Da | Smaller molecules absorb better |
        | **LogP** | < 5 | Lipophilicity — how well it dissolves in fats |
        | **TPSA** | < 140 Å² | Polar surface area — affects membrane crossing |
        | **HBD** | < 5 | H-bond donors — affects solubility |
        | **HBA** | < 10 | H-bond acceptors — affects permeability |
        | **GI Absorb.** | High | Gastrointestinal absorption prediction |
        """)

    with col2:
        st.markdown("#### 🟣 Distribution & Safety")
        st.markdown("""
        | Property | Ideal | Meaning |
        |----------|-------|---------|
        | **BBB** | Depends | Blood-Brain Barrier penetration |
        | **CYP2D6** | Unlikely | Risk of liver enzyme inhibition |
        | **QED** | > 0.5 | Quantitative drug-likeness (0-1) |
        | **Bioavail.** | > 0.5 | Fraction absorbed into bloodstream |
        | **Rot. Bonds** | < 10 | Molecular flexibility |
        | **Lipinski** | Pass | Overall drug-likeness rule |
        """)

    # 5. Best/Worst candidates
    st.markdown("---")
    st.subheader("🔑 Key Findings")

    # Best QED
    best_qed = df.loc[df['QED'].idxmax()]
    worst_qed = df.loc[df['QED'].idxmin()]

    col1, col2, col3 = st.columns(3)
    with col1:
        st.success(f"**Best Drug-Likeness (QED):**\nCandidate {best_qed['Candidate']}\nQED = {best_qed['QED']}")
    with col2:
        st.warning(f"**Lowest Drug-Likeness (QED):**\nCandidate {worst_qed['Candidate']}\nQED = {worst_qed['QED']}")
    with col3:
        lipinski_pass_count = len(df[df['Lipinski'] == "✅ Pass"])
        st.info(f"**Lipinski Compliance:**\n{lipinski_pass_count}/{len(df)} candidates pass\nthe Rule of Five")

    # 6. Candidate Cards
    st.markdown("---")
    st.subheader("🧪 Individual Candidate Profiles")
    for i in range(0, len(results), 2):
        cols = st.columns(2)
        for j in range(2):
            idx = i + j
            if idx < len(results):
                r = results[idx]
                mol = Chem.MolFromSmiles(r['SMILES'])
                with cols[j]:
                    with st.container(border=True):
                        c1, c2 = st.columns([1, 2])
                        with c1:
                            img = Draw.MolToImage(mol, size=(250, 250))
                            st.image(img, use_container_width=True)
                        with c2:
                            st.markdown(f"### Candidate {r['Candidate']}")
                            st.write(f"**QED:** {r['QED']} | **Lipinski:** {r['Lipinski']}")
                            st.write(f"**GI Absorption:** {r['GI Absorb.']} | **BBB:** {r['BBB']}")
                            st.write(f"**CYP2D6 Risk:** {r['CYP2D6']}")
                            st.write(f"**Bioavailability:** {r['Bioavail.']}")