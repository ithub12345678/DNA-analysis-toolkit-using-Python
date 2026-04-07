import streamlit as st
import matplotlib.pyplot as plt
import random
from Analysis import *

st.set_page_config(page_title="GC Analysis", layout="centered")

st.title("👩🏻‍🔬GC Content and ORF finding tool")

# User guide to run the app 
st.markdown("""### How to use this app: for GC Content Analysis:
1. Paste your DNA sequence in the text area below to Analyze GC Content and see the results by selecting in dropdown GC content analysis".
2. The app will display the GC content percentage and a bar graph of the results.
### Note:
- The input sequence should only contain the characters A, C, G, and T (case insensitive).
- The app will automatically clean the input sequence by removing any whitespace or newline characters and converting it to uppercase.
- If the sequence is too short (less than 5 nucleotides), a warning will be displayed, but the analysis can still be performed.""")

st.markdown("""### How to use this app: for ORF Finder Tool:
1. Paste your DNA sequence in the text area below to find ORFs and see the results by selecting in dropdown "ORF Finder Tool".
2. The app will display three options for ORF sequences found forward reading frames or reverse reading frames or both values simultaneously, with a tabular data
            and a bar graph of the results for deeper understanding.
""")

# -----------------------------
# 1. Input Section
# -----------------------------

if "seq" not in st.session_state:
    st.session_state.seq = ""



st.subheader("Input DNA Sequence",text_alignment="center")
example_seq = Test_seq = "".join(random.choices(["A", "C", "G", "T"], k=50))
if st.button("Example Sequence"):
    st.session_state.seq = example_seq

seq = st.text_area(
    "Paste your DNA sequence below:",
    value=st.session_state.seq,
    placeholder="Example: ATGCGCGTAACG...", 
    height=100
)

st.session_state.seq = seq

try :
    sliding_window_size = st.slider(
        "Sliding Window Size for GC Content Analysis",
        min_value=3,
        max_value=max(5,len(st.session_state.seq)),
        step=1,
        key="window_size"
        )
    st.subheader("Sequence Analysis Options",text_alignment="center")
    id = st.selectbox("Select a sequence to view its GC content:", list(["GC Content Analysis","ORF Finder Tool"]), key="ID_select",)
    if id == "GC Content Analysis":
            seq_calculations(seq,sliding_window_size)

    if id == "ORF Finder Tool":
            ORF_finder(seq)


except Exception as e:
    st.error(f"An error occurred: {e}")
    st.stop()


