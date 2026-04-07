import streamlit as st

# Page config
st.set_page_config(page_title="GC Content Analyzer", layout="centered")

import streamlit as st
st.write(st.experimental_get_pages())

# Title
st.title("🧬 GC-Content Analyzer")

# Intro text
st.markdown("### Welcome to GC-Content Analyzer")

# Q&A Section
st.markdown("#### Q1. So, what does this app do?")
st.write("ANS. The GC Content Analyzer is a dedicated tool designed to " \
"process user-inputted nucleotide sequences, counting each base " \
"to determine the proportion of Guanine (G) and Cytosine (C) — " \
"together referred to as GC content — and expressing it as a " \
"percentage. Alongside this, the app calculates the counts of " \
"Adenine and Thymine, presenting their combined AT percentage " \
"as well. These metrics give users a clear picture of the " \
"purine-to-pyrimidine ratio within their sequence. " \
"Additionally, the tool offers a graphical representation " \
"of the results, making the data easier to interpret and visually engage with.")

st.markdown("#### Q2. What are nucleotides?")
col1, col2 = st.columns(2)
with col1:
    st.image("Assets/1.png", width=400)
with col2:
    st.image("Assets/2.png", width=400)
st.write("ANS. A nucleotide is the basic building block of nucleic acids "
"(RNA and DNA). A nucleotide consists of a sugar molecule (either ribose in RNA or deoxyribose in DNA) " \
"attached to a phosphate group and a nitrogen-containing base. The bases used in DNA are adenine (A), cytosine (C), " \
"guanine (G) and thymine (T). In RNA, the base uracil (U) takes the place of thymine. DNA and RNA molecules are polymers made up of long chains of nucleotides.")

st.markdown("#### Q3. What is GC content and what is its significance?")
col7, col8 = st.columns(2)
with col7:
    st.image("Assets/6.jpg", width=400)
with col8:
    st.image("Assets/5.png", width=400)
st.write("ANS. The GC-content of a strand of nucleic acid is the percentage of nucleotides in the strand that possess either cytosine or " \
"guanine bases. For example, the GC-content of the RNA string [GAUCG] is 60%. " \
"In double-stranded DNA, every guanine base is complementary to the cytosine base in the opposite strand, " \
"so the GC-content of the two strands will be the same. count the number of G-C pairs in double helix. Because the GC-content " \
"throughout the genome differs between species, GC-content can offer a rough preliminary test of the identity of unknown DNA.\n "\
"\n Significance of GC-content: The GC-content of most species does tend to hover near 50%. However, coding regions of the genome have a tendency to contain a higher " \
"percentage of guanine and cytosine; these areas are called GC-rich, in contrast to areas of GC-content below 50%, which are called GC-poor. " \
"Thus, just as GC-content between species offers a rough test of species identity, testing the GC-content of a snippet of DNA " \
"from a known species can offer insight into whether that DNA may belong to a gene.")

# Transition text
st.markdown("""
Now that we understand the basics of DNA composition and GC content,  
let's move on to analyzing your sequence.""")

st.markdown("""### Choose Your Analysis Method \n
Based on the query type user can choose to analyze GC content from either a raw DNA sequence or from a FASTA file.
""")

Ready = st.toggle("Ready to analyze your sequence?", value=False)
col5, col6 = st.columns(2)
# Buttons
with col5:
    if st.button("☘︎ Analyze GC Content From Fasta file Sequence"):
        st.switch_page("pages/3_FastA_analysis.py")

with col6:
    if st.button("☘︎ Analyze GC Content From Raw Sequence"):
        st.switch_page("pages/2_Raw_analysis.py")
