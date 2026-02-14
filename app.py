import streamlit as st

# 1. Page Configuration
st.set_page_config(
    page_title="HIV-1 Research Login",
    page_icon="🧪",
    layout="centered",
    initial_sidebar_state="collapsed" # Hide sidebar until logged in
)

# 2. Custom CSS for the "Lab" aesthetic
st.markdown("""
    <style>
    .main {
        background-color: #f0f2f6;
    }
    .stButton>button {
        width: 100%;
        border-radius: 5px;
        height: 3em;
        background-color: #004b91;
        color: white;
    }
    </style>
    """, unsafe_allow_html=True)

# 3. Initialize Authentication State
if 'authenticated' not in st.session_state:
    st.session_state['authenticated'] = False

def login_page():
    st.title("🧪 HIV-1 Discovery Portal")
    st.markdown("### Secure Researcher Access")
    
    with st.container():
        # Create a clean login box
        with st.form("login_form"):
            res_id = st.text_input("Researcher ID", placeholder="Enter ID (e.g., ADMIN)")
            password = st.text_input("Security Key", type="password", placeholder="••••••••")
            
            submit = st.form_submit_button("Authorize Access")
            
            if submit:
                # We are hardcoding credentials for now since we have no DB yet
                if res_id == "ADMIN" and password == "PROTEASE2026":
                    st.session_state['authenticated'] = True
                    st.success("Access Granted! Proceed to 'Home' in the sidebar.")
                    st.balloons()
                else:
                    st.error("Access Denied: Invalid Researcher ID or Key.")

# 4. Logic to show content or login
if not st.session_state['authenticated']:
    login_page()
    # This keeps the other pages hidden/locked in a real deployment 
    # if you add logic to check auth at the top of every file.
else:
    st.success("Logged in successfully.")
    st.info("👈 Please select '1_Home' from the sidebar to begin the briefing.")