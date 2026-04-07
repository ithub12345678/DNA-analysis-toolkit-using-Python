Nucleotides = ['A', 'C', 'G', 'T']
import re
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt


codon_dict = {
    "F": ["UUU", "UUC"],
    "L": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
    "I": ["AUU", "AUC", "AUA"],
    "M": ["AUG"],  # Start codon
    "V": ["GUU", "GUC", "GUA", "GUG"],
    
    "S": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
    "P": ["CCU", "CCC", "CCA", "CCG"],
    "T": ["ACU", "ACC", "ACA", "ACG"],
    "A": ["GCU", "GCC", "GCA", "GCG"],
    
    "Y": ["UAU", "UAC"],
    "H": ["CAU", "CAC"],
    "Q": ["CAA", "CAG"],
    "N": ["AAU", "AAC"],
    "K": ["AAA", "AAG"],
    
    "D": ["GAU", "GAC"],
    "E": ["GAA", "GAG"],
    "C": ["UGU", "UGC"],
    "W": ["UGG"],
    
    "R": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "G": ["GGU", "GGC", "GGA", "GGG"],
    
    # Stop codons
    "_": ["UAA", "UAG", "UGA"] 
}

def clean_raw_seq(seq):
    """
    Cleans the raw DNA sequence by:
    - Removing invalid characters
    - Converting to uppercase
    """
    
    seq = seq.upper()  # Convert to uppercase
    seq = "".join(re.findall(r"[^\s\n]+", seq))  # Keep only letters

    for nucleotide in seq:
        if nucleotide not in Nucleotides:
            return False
        
    return seq

"""-------------------------------------------------------------"""

def clean_fasta_seq(file):
    """
    Parses FASTA and returns {id: cleaned_sequence}.
    Removes all whitespace and validates nucleotides.
    """
    try:
        matches = re.findall(r">(.+?)\n([^>]*)", file)

        clean_dict = {}

        for ids, seq_block in matches:
            # Remove ALL whitespace (newline, spaces, tabs)
            seq = re.sub(r"\s+", "", seq_block).upper()

            invalid = set(seq) - set(Nucleotides)
            if invalid:
                raise ValueError(f"Invalid characters in sequence {ids}: {invalid}")

            clean_dict[ids] = seq

        return clean_dict
    except Exception as e:
        st.error(f"An error occurred during FASTA file processing: {e} type error")
        st.stop()

"""-------------------------------------------------------------"""

def reverse_complement(seq):
    """Returns the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement[base] for base in reversed(seq))

"""-------------------------------------------------------------"""

def transcribe_dna_to_rna(seq):
    """Transcribes a DNA sequence to RNA by replacing 'T' with 'U'."""
    return seq.replace('T', 'U')

"""-------------------------------------------------------------"""

def translation(seqs):
        """Translates a list of RNA sequences into their corresponding protein sequences using the codon dictionary."""
        protein_seqs = []
        pro = []
        for orf in seqs:
            for i in range(0, len(orf) , 3):
                codon = orf[i:i+3]
                for key, value in codon_dict.items():
                   for j in value:
                        if j == codon:  
                            pro.append(key) 
                        else: pass
        protein_seqs.append("".join(pro))

        return protein_seqs

"""-------------------------------------------------------------"""


def orf_seq(seq):
    """Finds ORF sequences in a given protein sequence (after translation) by identifying sequences that start with 'M' and end with a stop codon ('_'). Returns a dictionary with the index of the start codon and the corresponding ORF sequence."""
    orfs_dict = {}
    if "M" in seq:
        for match in re.finditer("M",seq):
            matches = seq[match.start():]
            # print(match.group())
            for stop_codon in re.finditer("_",matches):
                stop_seq = matches[:stop_codon.end()]
                # print(stop_codon.end())
                orfs_dict[match.start(),match.start() + stop_codon.end()] = stop_seq.replace("_","")
                break
    return orfs_dict

"""--------------------------------------------------------------"""

def visualize_orfs(orf_seq_info,sign = "+" ,strand = "Forward"):
    """
    Processes the list of 3 dictionaries: 
    [{index: orf, ...}, {index: orf, ...}, {index: orf, ...}]
    """
    stats_data = []
    
    # Iterate through the 3 reading frames (+1, +2, +3)
    for i, frame_dict in enumerate(orf_seq_info):
        frame_label = f"{sign}{i+1}"
        
        # Get all ORF sequences (values) from this frame's dictionary
        orf_sequences = list(frame_dict.values())
        lengths = [len(s) for s in orf_sequences]
        
        count = len(orf_sequences)
        longest = max(lengths) if lengths else 0
        avg = round(sum(lengths) / count, 2) if count > 0 else 0
        
        stats_data.append({
            "Reading Frame": frame_label,
            "Total ORFs": count,
            "Longest (aa)": longest,
            "Avg Length": avg
        })

    # 1. Display the Table
    df = pd.DataFrame(stats_data)
    st.subheader(f"{strand} Strand ORF Statistics")
    st.table(df)

    # 2. Create the Bar Graph
    st.subheader(f"ORF Count Comparison: Forward Frames ({sign}1, {sign}2, {sign}3)")
    fig, ax = plt.subplots()
    
    # Plotting the count for each frame
    bars = ax.bar(df["Reading Frame"], df["Total ORFs"], color=['#3498db', '#9b59b6', '#2ecc71'])
    
    # Adding labels on top of bars
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, yval + 0.1, yval, ha='center', va='bottom')

    ax.set_ylabel("Number of ORFs Found")
    ax.set_xlabel("Reading Frame")
    ax.set_title(f"{strand} Strand ORF Distribution")
    
    st.pyplot(fig)

# Integration: Call this function after 'for' loop finishes
# visualize_orfs(orf_seq_info_forw, sign="+", strand="Forward")

"""--------------------------------------------------------------"""

def visualize_orfs_combined(combined_info):
    """
    Processes both Forward and Reverse ORF lists to create a 
    unified table and a 6-frame bar graph.
    """
    stats_data = []
    
    # Define labels for the 6 frames
    frame_labels = ["+1", "+2", "+3", "-1", "-2", "-3"]
    # Combine the lists into one (total 6 dictionaries)

    for i, frame_dict in enumerate(combined_info):
        label = frame_labels[i]
        
        # Extract sequences and calculate lengths
        orf_sequences = list(frame_dict.values())
        lengths = [len(s) for s in orf_sequences]
        
        count = len(orf_sequences)
        longest = max(lengths) if lengths else 0
        avg = round(sum(lengths) / count, 2) if count > 0 else 0
        
        stats_data.append({
            "Reading Frame": label,
            "Total ORFs": count,
            "Longest (aa)": longest,
            "Avg Length": avg
        })

    # 1. Create and display the table
    df = pd.DataFrame(stats_data)
    st.subheader("All-Strand ORF Statistics (+1 to -3)")
    st.table(df)

    # 2. Create the Bar Graph for all 6 frames
    st.subheader("ORF Count Distribution Across All Frames")
    fig, ax = plt.subplots(figsize=(10, 5))
    
    # Use distinct colors for Forward (+) and Reverse (-)
    colors = ['#3498db']*3 + ['#e74c3c']*3
    bars = ax.bar(df["Reading Frame"], df["Total ORFs"], color=colors)
    
    # Add data labels on top of bars
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, yval + 0.1, int(yval), ha='center', va='bottom')

    ax.set_ylabel("Number of ORFs")
    ax.set_xlabel("Reading Frame")
    st.pyplot(fig)

"""--------------------------------------------------------------"""

if __name__ == "__main__":
    pass
