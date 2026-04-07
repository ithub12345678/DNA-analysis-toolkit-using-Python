import streamlit as st
from Analysis import *

st.set_page_config(page_title="FASTA File GC Content Analyzer", layout="centered")
st.title("📁 FASTA File GC Content Analyzer")

st.write("This tool allows you to upload a FASTA file and analyze its GC content.")
st.write("Let's have a quick look how to use this tool:")

st.write("Before that lets have a quick qna session to understand what a FASTA file is and how it is structured.")
st.markdown("#### Q1. So, what Exactly is a FASTA file?")
st.write("ANS. A FASTA file is a widely used, plain text bioinformatics format designed to store nucleotide or protein sequences. It begins with a header line starting with a greater-than (" \
           "symbol, followed by lines of raw sequence data. It supports multiple sequences per file, often saved as .fasta, .fa, .fna, or .faa..")

st.markdown("#### Q2. What is the structure of a FASTA file?")
st.image("assets/6.png", width=400)
st.write("ANS. A FASTA file consists of one or more entries, each containing a header line and a sequence. The header line starts with a '>' character, followed by an identifier and optional description. " \
"The sequence lines contain the actual nucleotide or amino acid sequences, which can be split across multiple lines for readability.")

st.markdown("#### Q3. How can I analyze the GC content of a FASTA file using this tool?")
st.write("ANS. To analyze the GC content of a FASTA file using this tool, simply upload your FASTA file using the file uploader. The tool will read the file, extract the sequences, and calculate the GC content for each sequence.")
st.write("The results will be displayed in a clear format, showing the GC content percentage for each sequence, along with any relevant statistics or visualizations to help you understand the composition of your sequences.")


st.subheader("Let's Analyze Your FASTA File",text_alignment="center")
st.write("--"*20)
file = st.file_uploader("Upload FastA file: ", type=["fasta", "fa","fas","fsa","ffn","txt","csv"])
if file:
    st.write("File uploaded successfully! Processing the file...")
    # You can add code here to read and process the FASTA file using your existing functions from gc_analysis.py and parser.py
    read_file = file.read().decode("utf-8")
    if st.button("View File Content"):
        st.text_area("File Content:", value = read_file, height=150, disabled=True)
    
    if st.button("Analyze GC Content from FASTA File"):
        analyze_fasta(read_file)
else:   
    st.info("Please upload a FASTA file to analyze its GC content.")
