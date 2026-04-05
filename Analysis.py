from parser import *
from collections import Counter
import matplotlib.pyplot as plt
import streamlit as st
import parser



def seqs_calculation(file):
    """Performs GC content analysis."""
    try :
        seq_dict = parser.clean_fasta_seq(file)
        gc_percent = {}
        at_percent = {}
        nucleotide_counts = {}
        for key,seq in seq_dict.items():
             
            if len(seq) == 0:
                gc_percent[key] = 0
                at_percent[key] = 0
                nucleotide_counts[key] = {}
                continue

            count = Counter(seq)
            nucleotide_counts[key] = count

            gc = ((count['G'] + count['C']) / len(seq)) * 100
            at = ((count['A'] + count['T']) / len(seq)) * 100

            gc_percent[key] = round(gc, 2)
            at_percent[key] = round(at, 2)

    except Exception as e:
        # st.error(f"An error occurred during FASTA file processing: {e}")
        st.stop()

    return nucleotide_counts, gc_percent, at_percent, seq_dict

"""-------------------------------------------------------------"""

def slider(info,id,min_val,max_val,default_window,step_size):
    """Provides a slider for selecting the size of the sliding window for GC content analysis and visualizes the GC content across the sequence splits."""
    N = len(info[3][id])
    sliding_window_size = st.slider("Sliding Window Size for GC Content Analysis", min_value= min_val, max_value = max(5,max_val), step=step_size, value = default_window)

    if sliding_window_size > N:
                    st.error(f"GC Content can't be analyzed for sliding windows of size {sliding_window_size}.")
    else:
            splits = []
            splits_gc_percent = []
            for i in range(0,N,sliding_window_size):
                slide = info[3][id][i:i+sliding_window_size]
                if len(slide) >= sliding_window_size: 
                    splits.append(slide)
                    splits_gc_percent.append(gc_at_Percentage_raw(slide)[1])
                
            fig, ax = plt.subplots()
            plt.plot(range(len(splits)), splits_gc_percent, color='orange')
            ax.set_title(f'GC Content for Sliding Windows of Size {sliding_window_size}')
            ax.set_xlabel('Split Position')
            ax.set_ylabel('GC Content (%)')
            plt.tick_params(axis='x', labelbottom=False)
            st.pyplot(fig)
    
    return sliding_window_size,splits, splits_gc_percent

"""-------------------------------------------------------------"""

@st.fragment
def analyze_fasta(file):
            """Analyzes the GC content of sequences in a FASTA file and visualizes the results."""
            info = seqs_calculation(file)
            st.markdown("### GC content analysis completed! Here are the results:")
            id = st.selectbox("Select a sequence to view its GC content:", list(info[3].keys()), key="ID_select")
            N = len(info[3][id])
            st.write(f"Sequence Length: {N} nucleotides")
            st.write(f"Nucleotide Counts:{', '.join(f'{k}: {v}' for k, v in info[0][id].items())}")
            st.write(f"GC Content: {info[1][id]}%")
            st.write(f"AT Content: {info[2][id]}%")

            # -----------------------------
            # 3. Visualization Section
            # -----------------------------
            st.subheader("Nucleotide Count Visualization")
            fig , ax = plt.subplots()
            for i in range(len(info[0][id].keys())):
                plt.text(i, list(info[0][id].values())[i], list(info[0][id].values())[i], ha='center')
            ax.bar(info[0][id].keys(), info[0][id].values(),color=['blue', 'yellow', 'green', 'red'])
            ax.set_title('Nucleotide Counts')
            ax.set_xlabel('Nucleotide')
            ax.set_ylabel('Count')
            st.pyplot(fig)

            st.subheader("Select the Plot type for GC% visulazation across sequence splits")
            mode = st.selectbox(
                        "Select GC% Visualization Mode",
                        ["Detailed", "Smoothed"]
                    )
            
            if mode == "Detailed":
                if N < 100:
                    slider(info,id,3, N, max(3, N//20), 1)

                elif 100 <= N <= 1000:
                    slider(info,id,3, N//10, N//25, 1)

                else:  # N > 1000
                    slider(info,id,3, N//20, N//50, 1)
            
            elif mode == "Smoothed":
                if N < 100:
                    slider(info,id,5, N, max(5, N//10), 5)

                elif 100 <= N <= 1000:
                    slider(info,id,5, N//10, N//20, 5)

                else:  # N > 1000
                    slider(info,id,5, N//20, N//50, 5)

"""-------------------------------------------------------------"""

def gc_at_Percentage_raw(seq):
    """Percentage of nucleotides in a raw DNA sequence."""
    
    if len(seq) == 0:
        return {}, 0, 0

    count = Counter(seq)

    gc = ((count['G'] + count['C']) / len(seq)) * 100
    at = ((count['A'] + count['T']) / len(seq)) * 100

    return count, round(gc, 2), round(at, 2)

"""-------------------------------------------------------------"""

def reading_frames(n,rna_seq):
    frames = ""
    frames_rna = []
    for i in range(n, len(rna_seq) , 3):
        if len(rna_seq[i:i+3]) >= 3 : frames+=(rna_seq[i:i+3])
            
    frames_rna.append("".join(frames))
    return frames_rna

"""-------------------------------------------------------------"""

@st.fragment
def seq_calculations(seq,window_size):
        """Performs GC content analysis on a raw DNA sequence and visualizes the results."""
        if seq:
            # st.write(f"Your input sequence is: {seq}")
            seq = parser.clean_raw_seq(seq)
            if not seq:
                st.error("Invalid sequence! Please enter a valid DNA sequence containing only A, C, G, and T.")
            else:
                if len(seq) <5 : st.warning(f"Warning: Sequence length is {len(seq)}")
                st.success("Sequence cleaned!")
                st.text_area("Cleaned Sequence:", value=seq, height=50,disabled=True)
                # -----------------------------
                # 2. Analysis Section
                # -----------------------------
                st.subheader("GC Content Analysis Results")
                count, gc_percent, at_percent = gc_at_Percentage_raw(seq)
                st.write(f"Nucleotide Counts:  {', '.join(f'{k}: {v}' for k, v in count.items())}")
                st.write(f"GC Content: {gc_percent}%")
                st.write(f"AT Content: {at_percent}%")

            
                if window_size > len(seq):
                    st.error(f"GC Content can't be analyzed for sliding windows of size {window_size}.")
                else:
                    splits = []
                    splits_gc_percent = []
                    for i in range(0,len(seq),window_size):
                        slide = seq[i:i+window_size]
                        if len(slide) >= window_size: 
                            splits.append(slide)
                            splits_gc_percent.append(gc_at_Percentage_raw(slide)[1])
                    # -----------------------------
                    # 3. Visualization Section
                    # -----------------------------
                    st.subheader("Nucleotide Count Visualization")
                    fig , ax = plt.subplots()
                    for i in range(len(count.keys())):
                        plt.text(i, list(count.values())[i], list(count.values())[i], ha='center')
                    ax.bar(count.keys(), count.values(),color=['blue', 'yellow', 'green', 'red'])
                    ax.set_title('Nucleotide Counts')
                    ax.set_xlabel('Nucleotide')
                    ax.set_ylabel('Count')
                    st.pyplot(fig)

                    st.subheader("GC Content Across Sliding Windows")
                    ax.clear()
                    fig, ax = plt.subplots()
                    for i in range(len(splits)):
                        plt.text(i, splits_gc_percent[i], splits_gc_percent[i], ha='center')
                    ax.bar(splits, splits_gc_percent, color='orange')
                    ax.set_title(f'GC Content in Sliding Windows of Size {window_size}')
                    ax.set_xlabel('Sequence Window')
                    ax.set_ylabel('GC Content (%)')
                    plt.tick_params(axis='x', labelbottom=False)
                    # plt.xticks(rotation=45, ha='right')
                    st.pyplot(fig)

"""-------------------------------------------------------------"""
@st.fragment
def ORF_finder(seq):
    if seq:
            # st.write(f"Your input sequence is: {seq}")
            seq = parser.clean_raw_seq(seq)
            if not seq:
                st.error("Invalid sequence! Please enter a valid DNA sequence containing only A, C, G, and T.")
            else:
                if len(seq) <5 : st.warning(f"Warning: Sequence length is {len(seq)}")
                st.success("Sequence cleaned!")
                st.text_area("Cleaned Sequence:", value=seq, height=50,disabled=True)
                reverse_complement_seq = reverse_complement(seq)
                st.write(f"The reverse complement of the sequence is: {reverse_complement_seq}")
                for_transcibed_seq = parser.transcribe_dna_to_rna(seq)
                rev_transcibed_seq = parser.transcribe_dna_to_rna(reverse_complement_seq)
                st.write(f"The RNA transcript of the forward sequence is: {for_transcibed_seq}")
                st.write(f"The RNA transcript of the reverse sequence is: {rev_transcibed_seq}")

                st.write("--"*20)

                # -----------------------------
                # 2. Analysis Section
                # -----------------------------

                strand = st.radio(
                "Select the strand for ORF analysis:",
                ["***:rainbow[FORWARD STRAND]***", "***:rainbow[REVERSE STRAND]***", "***:rainbow[ BOTH STRANDS ]***"],
                captions=[
                    "Plus strand.",
                    "Minus strand.",
                    "Both strands."])

                if strand == "***:rainbow[FORWARD STRAND]***":

                    st.subheader("The Forward Open Reading Frames (ORFs) with protein sequences are:")
                    orf_seq_info_forw = []
                    for i in range(3):
                        st.write(f"The Forward Reading Frames are: +{i+1} {"".join(reading_frames(i,for_transcibed_seq))}")
                        st.write(f"The Protein sequence is: {"".join(parser.translation(reading_frames(i,for_transcibed_seq)))}")
                        protein_seq = "".join(parser.translation(reading_frames(i,for_transcibed_seq)))
                        st.write(f"The ORF sequence result is specifying the index and the orf: {orf_seq(protein_seq)}")
                        orf_seq_info_forw.append(orf_seq(protein_seq))
                        st.write("--"*20)

                    # -----------------------------
                    # 3. Visualization Section
                    # -----------------------------

                    visualize_orfs(orf_seq_info_forw, sign = "+", strand = "Forward")

                elif strand == "***:rainbow[REVERSE STRAND]***":

                    st.subheader("The Reverse Open Reading Frames (ORFs) with protein sequences are:")
                    orf_seq_info_rev = []
                    for i in range(3):
                        st.write(f"The Reverse Reading Frames are: -{i+1} {"".join(reading_frames(i,rev_transcibed_seq))}")
                        st.write(f"The Protein sequence is: {"".join(parser.translation(reading_frames(i,rev_transcibed_seq)))}")
                        protein_seq = "".join(parser.translation(reading_frames(i,rev_transcibed_seq)))
                        st.write(f"The ORF sequence result is specifying the index and the orf: {orf_seq(protein_seq)}")
                        orf_seq_info_rev.append(orf_seq(protein_seq))
                        st.write("--"*20)

                    # -----------------------------
                    # 3. Visualization Section
                    # -----------------------------

                    visualize_orfs(orf_seq_info_rev, sign = "-", strand = "Reverse")

                else:
                    st.subheader("The Forward and Reverse Open Reading Frames (ORFs) with protein sequences are:")
                    orf_seq_info_forw = []
                    orf_seq_info_rev = []
                    
                    for i in range(3):
                        st.write(f"The Forward Reading Frame is: +{i+1} {''.join(reading_frames(i,for_transcibed_seq))}")
                        st.write(f"The Protein sequence is: {''.join(parser.translation(reading_frames(i,for_transcibed_seq)))}")
                        protein_seq = ''.join(parser.translation(reading_frames(i,for_transcibed_seq)))
                        st.write(f"The ORF sequence result is specifying the index and the orf: {orf_seq(protein_seq)}")
                        orf_seq_info_forw.append(orf_seq(protein_seq))
                        st.write("--"*20)
                    
                    for i in range(3):
                        st.write(f"The Reverse Reading Frames are: -{i+1} {''.join(reading_frames(i,rev_transcibed_seq))}")
                        st.write(f"The Protein sequence is: {''.join(parser.translation(reading_frames(i,rev_transcibed_seq)))}")
                        protein_seq = ''.join(parser.translation(reading_frames(i,rev_transcibed_seq)))
                        st.write(f"The ORF sequence result is specifying the index and the orf: {orf_seq(protein_seq)}")
                        orf_seq_info_rev.append(orf_seq(protein_seq))
                        st.write("--"*20)

                    
                    # -----------------------------
                    # 3. Visualization Section
                    Combined_orf_info = orf_seq_info_forw + orf_seq_info_rev
                    
                    visualize_orfs_combined(Combined_orf_info)
                    All_values = []
                    for i in range(len(Combined_orf_info)):
                        for j in Combined_orf_info[i].values():
                            All_values.append("".join(j))
                    
                    st.write(f"All the unique ORF sequences found in both strands are: {", ".join(set(All_values))}")

"""-------------------------------------------------------------"""

if __name__ == "__main__":
    pass
