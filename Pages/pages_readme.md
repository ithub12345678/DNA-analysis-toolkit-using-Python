📄 Pages: Analysis & Visualization

This directory contains the core analytical modules of the project. It focuses on biological sequence processing and exploratory analysis of raw datasets prior to downstream workflows.

📂 Folder Structure
File	Description
fasta_analysis.py	Performs parsing and analysis of FASTA files (DNA/protein sequences). Includes sequence extraction, length computation, and gc content detection logic.
raw_analysis.py	Conducts initial data preprocessing and exploratory data analysis (EDA) on raw datasets for gc content analysis and ORF reads finding between strands.

🚀 Getting Started
Prerequisites

Ensure the required dependencies are installed:

pip install streamlit pandas matplotlib
Usage
1. FASTA Sequence Analysis

Run sequence-level analysis on a FASTA file:

python pages/fasta_analysis.py --input data/sequences.fasta

2. Raw Data Analysis

Perform baseline preprocessing and exploratory analysis:

python pages/raw_analysis.py 

📊 Methodology
FASTA Processing
Utilizes the python library for efficient and reliable parsing of FASTA files
Extracts sequence identifiers and biological sequences
Computes sequence length and enables orf/pattern detection
Raw Data Pipeline
Identifies and handles missing or inconsistent values
Generates descriptive statistics for feature distributions
Supports early-stage insight generation prior to advanced analysis
Streamlit for web development using python language offering various kinds of functionalities
Future task, to fix some minor bugs and deploy using streamlit
⚙️ Notes
Ensure input files follow correct formatting standards (FASTA for sequence data, txt or fasta format for raw datasets)
Scripts are modular and can be integrated into larger bioinformatics pipelines or visualization layers (e.g., Streamlit dashboards)
