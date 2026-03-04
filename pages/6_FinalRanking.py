import streamlit as st
import pandas as pd
import sys
import os
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, QED, RDConfig

# SA Score
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

# Add parent directory so we can import candidates
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from candidates import CANDIDATES

# 1. Security Check
if 'authenticated' not in st.session_state or not st.session_state['authenticated']:
    st.warning("Please login on the main page first.")
    st.stop()

st.set_page_config(page_title="Final Ranking", page_icon="🏆", layout="wide")

st.title("🏆 Final Drug Candidate Ranking")
st.markdown("### Comprehensive Multi-Criteria Evaluation")
st.markdown("""
This page combines **all evaluation metrics** from the previous analyses into a single 
**weighted composite score** to determine the most promising drug candidate for HIV-1 Protease.
""")
st.markdown("---")

# 2. Pre-computed docking results (same as 5_Docking.py)
DOCKING_RESULTS = {
    1: -8.4, 2: -7.1, 3: -9.2, 4: -6.5,
    5: -9.8, 6: -8.7, 7: -7.8, 8: -5.9,
}

# 3. Scoring weights
WEIGHTS = {
    "AI Score": 0.20,
    "Stability (Gap)": 0.10,
    "Docking": 0.30,
    "ADMET (QED)": 0.25,
    "Synthesizability": 0.15,
}


def normalize(value, min_val, max_val, invert=False):
    """Normalize a value to 0-1 range. If invert=True, lower original values get higher scores."""
    if max_val == min_val:
        return 0.5
    norm = (value - min_val) / (max_val - min_val)
    if invert:
        norm = 1.0 - norm
    return max(0.0, min(1.0, norm))


@st.cache_data
def compute_final_ranking():
    raw_data = []
    for c in CANDIDATES:
        mol = Chem.MolFromSmiles(c['smiles'])
        if mol is None:
            continue

        sa_score = sascorer.calculateScore(mol)
        qed_score = QED.qed(mol)
        docking = DOCKING_RESULTS[c['id']]

        raw_data.append({
            "id": c['id'],
            "smiles": c['smiles'],
            "ai_score": c['score'],
            "egap": c['egap'],
            "docking": docking,
            "qed": qed_score,
            "sa_score": sa_score,
        })

    # Compute normalization ranges
    ai_scores = [d['ai_score'] for d in raw_data]
    egaps = [d['egap'] for d in raw_data]
    dockings = [d['docking'] for d in raw_data]
    qeds = [d['qed'] for d in raw_data]
    sa_scores = [d['sa_score'] for d in raw_data]

    results = []
    for d in raw_data:
        # Normalize each metric to 0-1 (higher = better)
        norm_ai = normalize(d['ai_score'], min(ai_scores), max(ai_scores))
        norm_gap = normalize(d['egap'], min(egaps), max(egaps), invert=True)  # smaller gap = more reactive = better
        norm_dock = normalize(d['docking'], min(dockings), max(dockings))  # more negative = better, so min is best
        norm_qed = normalize(d['qed'], min(qeds), max(qeds))
        norm_sa = normalize(d['sa_score'], min(sa_scores), max(sa_scores), invert=True)  # lower SA = better

        # Weighted composite
        composite = (
            WEIGHTS["AI Score"] * norm_ai +
            WEIGHTS["Stability (Gap)"] * norm_gap +
            WEIGHTS["Docking"] * norm_dock +
            WEIGHTS["ADMET (QED)"] * norm_qed +
            WEIGHTS["Synthesizability"] * norm_sa
        )

        results.append({
            "Candidate": f"#{d['id']}",
            "SMILES": d['smiles'],
            "AI Score": f"{int(d['ai_score'] * 100)}%",
            "HOMO-LUMO Gap": f"{d['egap']} eV",
            "Docking": f"{d['docking']} kcal/mol",
            "QED": round(d['qed'], 3),
            "SA Score": round(d['sa_score'], 2),
            "Composite": round(composite * 100, 1),
        })

    df = pd.DataFrame(results)
    df = df.sort_values("Composite", ascending=False).reset_index(drop=True)
    df.index = df.index + 1  # Rank starts from 1
    df.index.name = "Rank"
    return df


# 4. Weight configuration
st.subheader("⚖️ Scoring Weights")
st.markdown("The final composite score is calculated using these weights:")

weight_cols = st.columns(5)
for i, (metric, weight) in enumerate(WEIGHTS.items()):
    with weight_cols[i]:
        with st.container(border=True):
            st.metric(metric, f"{int(weight * 100)}%")


if st.button("🏆 Compute Final Rankings"):
    st.session_state['ranking_computed'] = True

if st.session_state.get('ranking_computed', False):
    with st.spinner("Computing composite scores..."):
        df = compute_final_ranking()

    st.success("Final Rankings Computed!")

    # 5. Composite Score Chart
    st.markdown("---")
    st.subheader("📊 Composite Score Ranking")
    chart_data = df[["Candidate", "Composite"]].set_index("Candidate")
    st.bar_chart(chart_data, color="#FF9800", horizontal=True)

    # 6. Full Rankings Table
    st.subheader("📋 Comprehensive Rankings")
    display_cols = ["Candidate", "AI Score", "HOMO-LUMO Gap", "Docking", "QED", "SA Score", "Composite"]
    st.dataframe(
        df[display_cols],
        use_container_width=True,
        column_config={
            "Composite": st.column_config.ProgressColumn(
                "Composite Score",
                min_value=0,
                max_value=100,
                format="%.1f",
            ),
        }
    )

    # 7. Winner Announcement
    st.markdown("---")
    winner = df.iloc[0]
    runner_up = df.iloc[1]

    st.subheader("🥇 Recommended Lead Compound")
    col1, col2 = st.columns([1, 2])
    with col1:
        mol = Chem.MolFromSmiles(winner['SMILES'])
        if mol:
            img = Draw.MolToImage(mol, size=(400, 400))
            st.image(img, use_container_width=True)
    with col2:
        st.markdown(f"## Candidate {winner['Candidate']}")
        st.markdown(f"**Composite Score: {winner['Composite']}%**")
        st.markdown(f"""
        | Metric | Value |
        |--------|-------|
        | AI Score | {winner['AI Score']} |
        | HOMO-LUMO Gap | {winner['HOMO-LUMO Gap']} |
        | Docking Affinity | {winner['Docking']} |
        | Drug-Likeness (QED) | {winner['QED']} |
        | SA Score | {winner['SA Score']} |
        """)
        st.success("✅ This candidate shows the best overall balance of binding affinity, drug-likeness, and synthesizability.")

    # 8. Runner up
    st.subheader("🥈 Runner-Up")
    col1, col2 = st.columns([1, 2])
    with col1:
        mol2 = Chem.MolFromSmiles(runner_up['SMILES'])
        if mol2:
            img2 = Draw.MolToImage(mol2, size=(400, 400))
            st.image(img2, use_container_width=True)
    with col2:
        st.markdown(f"## Candidate {runner_up['Candidate']}")
        st.markdown(f"**Composite Score: {runner_up['Composite']}%**")
        st.markdown(f"""
        | Metric | Value |
        |--------|-------|
        | AI Score | {runner_up['AI Score']} |
        | HOMO-LUMO Gap | {runner_up['HOMO-LUMO Gap']} |
        | Docking Affinity | {runner_up['Docking']} |
        | Drug-Likeness (QED) | {runner_up['QED']} |
        | SA Score | {runner_up['SA Score']} |
        """)

    # 9. Conclusion
    st.markdown("---")
    st.subheader("📝 Conclusion")
    st.markdown(f"""
    After comprehensive evaluation across **five key metrics** — AI prediction confidence, 
    electronic stability (HOMO-LUMO gap), molecular docking affinity against HIV-1 Protease, 
    ADMET drug-likeness (QED), and synthetic accessibility — **Candidate {winner['Candidate']}** 
    emerges as the most promising lead compound for further development.
    
    **Recommended next steps:**
    1. Experimental validation through *in vitro* protease inhibition assays
    2. Lead optimization to improve any weak properties
    3. Scale-up synthesis feasibility study
    4. Extended toxicity profiling
    """)