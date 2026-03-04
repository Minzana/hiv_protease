import streamlit as st
import pandas as pd
import sys
import os
from rdkit import Chem
from rdkit.Chem import Draw, RDConfig

# Add the SA_Score contrib module to path
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

# Add parent directory so we can import candidates
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from candidates import CANDIDATES

# 1. Security Check
if 'authenticated' not in st.session_state or not st.session_state['authenticated']:
    st.warning("Please login on the main page first.")
    st.stop()

st.set_page_config(page_title="Synthesizability Analysis", page_icon="🏭", layout="wide")

st.title("🏭 Synthesizability Analysis")
st.markdown("### Synthetic Accessibility (SA) Score")
st.markdown("""
The **SA Score** predicts how easy or difficult it would be to actually **synthesize (manufacture)** 
each drug candidate in a chemistry laboratory. It is calculated using RDKit's fragment-based approach.

- **Score Range**: 1 (very easy) → 10 (very difficult)
- A good drug candidate should ideally have an SA Score **below 5**.
""")
st.markdown("---")

# 2. Compute SA Scores
@st.cache_data
def compute_sa_scores():
    results = []
    mols = []
    for c in CANDIDATES:
        mol = Chem.MolFromSmiles(c['smiles'])
        if mol:
            sa_score = sascorer.calculateScore(mol)
            results.append({
                "Candidate": f"#{c['id']}",
                "SMILES": c['smiles'],
                "SA Score": round(sa_score, 2),
                "Difficulty": classify_difficulty(sa_score),
                "AI Score": f"{int(c['score'] * 100)}%",
            })
            mols.append(mol)
    return results, mols

def classify_difficulty(score):
    if score <= 3:
        return "🟢 Easy"
    elif score <= 5:
        return "🟡 Moderate"
    elif score <= 7:
        return "🟠 Difficult"
    else:
        return "🔴 Very Difficult"

if st.button("🔬 Analyze Synthesizability"):
    st.session_state['sa_computed'] = True

if st.session_state.get('sa_computed', False):
    with st.spinner("Computing Synthetic Accessibility Scores..."):
        results, mols = compute_sa_scores()

    df = pd.DataFrame(results)
    df_sorted = df.sort_values("SA Score").reset_index(drop=True)

    st.success(f"Analysis Complete: {len(results)} candidates evaluated.")

    # 3. Bar Chart
    st.subheader("📊 SA Score Ranking (Lower = Better)")
    chart_data = df_sorted[["Candidate", "SA Score"]].set_index("Candidate")
    st.bar_chart(chart_data, color="#4CAF50", horizontal=True)

    # 4. Detailed Results Table
    st.subheader("📋 Detailed Results")
    st.dataframe(
        df_sorted,
        use_container_width=True,
        hide_index=True,
        column_config={
            "SA Score": st.column_config.ProgressColumn(
                "SA Score",
                min_value=0,
                max_value=10,
                format="%.2f",
            ),
        }
    )

    # 5. Show molecules with their scores
    st.subheader("🧪 Candidate Structures")
    for i in range(0, len(results), 4):
        cols = st.columns(4)
        for j in range(4):
            idx = i + j
            if idx < len(results):
                r = df_sorted.iloc[idx]
                mol = Chem.MolFromSmiles(r['SMILES'])
                with cols[j]:
                    with st.container(border=True):
                        st.markdown(f"**Candidate {r['Candidate']}**")
                        img = Draw.MolToImage(mol, size=(300, 300))
                        st.image(img, use_container_width=True)
                        st.metric("SA Score", f"{r['SA Score']}")
                        st.write(f"**{r['Difficulty']}**")

    # 6. Key Findings
    st.markdown("---")
    st.subheader("🔑 Key Findings")
    best = df_sorted.iloc[0]
    worst = df_sorted.iloc[-1]
    col1, col2 = st.columns(2)
    with col1:
        st.success(f"**Easiest to Synthesize:** Candidate {best['Candidate']} (SA Score: {best['SA Score']})")
    with col2:
        st.error(f"**Hardest to Synthesize:** Candidate {worst['Candidate']} (SA Score: {worst['SA Score']})")

    avg_score = df_sorted['SA Score'].mean()
    st.info(f"**Average SA Score:** {avg_score:.2f} — {'Most candidates are reasonably synthesizable.' if avg_score < 5 else 'Some candidates may be challenging to synthesize.'}")