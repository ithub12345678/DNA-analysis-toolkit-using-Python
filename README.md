# 🧬 DNA Analysis Toolkit (GC Content + ORF Finder)

## Overview

This project is a **bioinformatics pipeline built in Python** that performs fundamental DNA sequence analysis, including:

* GC content analysis (global + sliding window)
* ORF (Open Reading Frame) detection
* Sequence parsing and preprocessing
* Basic biological insights through statistics and visualization

It is designed as a **learning-oriented yet practical tool**, combining computational techniques with biological interpretation.

---

## Features

### 1. GC Content Analysis

* Calculates:

  * GC percentage
  * AT percentage
  * Nucleotide distribution
* Supports **sliding window analysis** for local GC variation
* Visualization graphs specifying nucleotide count and GC% sliding window real-time metric.
* Useful for identifying:

  * Gene-rich regions
  * GC bias patterns

---

### 2. ORF Finder

* Detects Open Reading Frames using:

  * Start codon: `AUG` (Methionine)
  * Stop codons: `UAA`, `UAG`, `UGA`
* Supports:

  * Frame-wise analysis (via translation approach)
  * Extraction of:

    * ORF sequences
    * Start and end positions
* Outputs:

  * All detected ORFs
  * Longest ORF (optional extension)
  * Detailed analysis table with helpful ORF related insights
  * Count comparison between all the 6 reading frames

---

### 3. Sequence Handling

* FASTA parsing
* Cleaning and validation of sequences
* Multi-sequence support

---

## Methodology

### ORF Detection Strategy

Instead of scanning nucleotides directly, the pipeline:

1. Generates **3 forward reading frames**
2. Translates each into a **protein sequence**
3. Identifies ORFs using:

   * `M` → Start
   * `_` → Stop
4. Maps amino acid positions back to nucleotide coordinates

This ensures:

* Frame correctness
* Simpler logic
* Improved computational clarity

---

### GC Content Calculation

* Global GC%:
  [
  GC% = \frac{G + C}{A + T + G + C} \times 100
  ]

* Sliding window analysis:

  * User-defined window size
  * Step-wise GC variation across sequence

---

## Output

### ORF Results

* ORF sequences (protein)
* ORF positions (start, end)
* Frame information (optional extension)

### GC Analysis

* GC%
* AT%
* Nucleotide counts

### Visualizations (if enabled)

* GC content plot (sliding window)
* ORF distribution (future extension)
* Length distribution (future extension)

---

## Project Structure

```
├── Home.py
├── parser.py
├── Analysis.py
├── utils/
├── Pages/
└── README.md
```

---

## Installation

```bash
git clone https://github.com/your-username/dna-analysis-toolkit.git
cd dna-analysis-toolkit
pip install -r requirements.txt
```

---

## Usage

### Run the tool

```bash
python Home.py
```

### Input

* FASTA file OR raw DNA sequence

### Output

* Console-based results
* (Optional) Streamlit interface for visualization

---

## Future Improvements

* File input for ORF detection
* Interactive visualization (Streamlit dashboard)
* Codon usage analysis
* Protein property analysis (molecular weight, pI)
* Basic motif detection feature 
* Export results (CSV / JSON / FASTA)

---

## Learning Objectives

This project demonstrates:

* Sequence parsing and manipulation
* Reading frame logic in genomics
* Efficient pattern detection
* Integration of biology with programming

---

## Limitations

* Simplified ORF detection (no introns/exons handling)
* Uses standard genetic code only

---

## Contributing

Contributions are welcome. You can:

* Improve performance
* Add new biological analyses
* Enhance visualization

---

## License

This project is open-source and available under the MIT License.

---

## Author

Developed as part of a **bioinformatics learning pipeline** combining Python programming with molecular biology concepts.
