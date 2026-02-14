import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, AllChem
from stmol import showmol
import py3Dmol

# 1. Security Check
if 'authenticated' not in st.session_state or not st.session_state['authenticated']:
    st.warning("Please login on the main page first.")
    st.stop()

st.set_page_config(page_title="Discovery Engine", page_icon="⚡", layout="wide")

# Initialize session states
if 'show_results' not in st.session_state:
    st.session_state['show_results'] = False
if 'view_3d_smiles' not in st.session_state:
    st.session_state['view_3d_smiles'] = None

# 2. Candidate Data
candidates = [
    {"id": 1, "smiles": "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5", "score": 0.82, "egap": 8.12},
    {"id": 2, "smiles": "CC(C)(C)C1=CC=C(C=C1)C(=O)NC2=CC=CC=C2C(=O)O", "score": 0.75, "egap": 7.45},
    {"id": 3, "smiles": "CN1CCN(CC1)CC2=CC=C(C=C2)C(=O)NC3=CC=C(C=C3)C4=CSC=N4", "score": 0.88, "egap": 9.20},
    {"id": 4, "smiles": "CC1=CC=C(C=C1)S(=O)(=O)NC2=CC=CC=C2", "score": 0.69, "egap": 6.80},
    {"id": 5, "smiles": "C1=CC=C(C=C1)C2=C(C3=CC=CC=C3N2)C(=O)NC4=CC=NC=C4", "score": 0.96, "egap": 10.49},
    {"id": 6, "smiles": "CC1=NC2=C(N1)C=C(C=C2)C(=O)NC3=CC=CC=C3F", "score": 0.84, "egap": 8.90},
    {"id": 7, "smiles": "CNC(=O)C1=CC=CC=C1SC2=CC=C(C=C2)C(=O)C3=CC=CC=C3", "score": 0.79, "egap": 7.15},
    {"id": 8, "smiles": "C1=CC=C(C=C1)NC(=O)N2CCCCC2", "score": 0.72, "egap": 6.50},
]

def calculate_lipinski(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return {
            "MW": round(Descriptors.MolWt(mol), 2),
            "LogP": round(Descriptors.MolLogP(mol), 2),
            "HBD": Descriptors.NumHDonors(mol),
            "HBA": Descriptors.NumHAcceptors(mol)
        }
    return None

st.title("⚡ AI Discovery Engine: Lead Candidates")

# Trigger button sets a state variable
if st.button("Generate New Drug Candidates"):
    st.session_state['show_results'] = True

if st.session_state['show_results']:
    st.success("Analysis Complete: 8 Candidates Identified.")
    
    # 3D VIEWER AREA (Top of results if active)
    if st.session_state['view_3d_smiles']:
        with st.container(border=True):
            col_a, col_b = st.columns([3, 1])
            with col_a:
                st.subheader("Interactive 3D Molecular Analysis")
                xyzview = py3Dmol.view(width=700, height=400)
                mol = Chem.MolFromSmiles(st.session_state['view_3d_smiles'])
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol)
                mblock = Chem.MolToMolBlock(mol)
                xyzview.addModel(mblock, 'mol')
                xyzview.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
                xyzview.addSurface(py3Dmol.VDW, {'opacity': 0.5, 'color': 'white'})
                xyzview.zoomTo()
                showmol(xyzview, height=400, width=700)
            with col_b:
                st.write("**3D Controls**")
                st.info("Rotate: Left Click\nZoom: Scroll\nPan: Right Click")
                if st.button("❌ Close 3D View"):
                    st.session_state['view_3d_smiles'] = None
                    st.rerun()

    # THE GRID
    for i in range(0, len(candidates), 4):
        cols = st.columns(4)
        for j in range(4):
            idx = i + j
            if idx < len(candidates):
                c = candidates[idx]
                lp = calculate_lipinski(c['smiles'])
                with cols[j]:
                    with st.container(border=True):
                        st.markdown(f"### Candidate #{c['id']}")
                        if c['id'] == 5: st.success("⭐ TOP PERFORMER")
                        
                        # 2D Image
                        mol_2d = Chem.MolFromSmiles(c['smiles'])
                        img = Draw.MolToImage(mol_2d, size=(300, 300))
                        st.image(img, use_container_width=True)
                        
                        st.metric("AI Score", f"{int(c['score']*100)}%")
                        st.metric("Stability", f"{c['egap']} eV")
                        
                        # Button to update session state
                        if st.button(f"Analyze 3D #{c['id']}", key=f"btn_{c['id']}"):
                            st.session_state['view_3d_smiles'] = c['smiles']
                            st.rerun()
                        
                        with st.expander("Lipinski Rules"):
                            st.write(f"MW: {lp['MW']}")
                            st.write(f"LogP: {lp['LogP']}")