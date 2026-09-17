# 🧬 Signal Peptide Prediction

## Laboratory of Bioinformatics II — Module 2

A computational study of **signal peptide prediction in eukaryotic proteins**, combining classical bioinformatics methods with machine-learning approaches.

---
## Progress update
- setup VM, VPN, and github 
- Data collection
Positive: 2972 retrieved → 2961 final
Negative: 20975 retrieved → 20975 final

---

## 📌 Project Overview

Signal peptides (SPs) are short N-terminal amino-acid sequences that direct newly synthesized proteins into the **secretory pathway**, enabling their transport through the endoplasmic reticulum (ER), Golgi apparatus, and ultimately to destinations such as the cell membrane or extracellular space.

Predicting signal peptides from protein sequences is an important problem in **protein function prediction, subcellular localization, and genome/proteome annotation**.

This project investigates different computational approaches for identifying signal peptides and their cleavage sites, starting from classical statistical methods and progressing toward machine-learning approaches.

> **Main question:**  
> Can we reliably distinguish proteins containing signal peptides from proteins without signal peptides using sequence information?

---

## 🎯 Objectives

- Collect and curate a high-quality **eukaryotic protein dataset** from UniProtKB.
- Construct positive and negative datasets for signal peptide prediction.
- Analyze sequence and dataset characteristics.
- Extract informative sequence features.
- Implement the **von Heijne weight-matrix approach** for signal peptide/cleavage-site prediction.
- Implement a **feature-based Support Vector Machine (SVM)** approach.
- Evaluate the models using cross-validation and an independent/blind test set.
- Compare the approaches and analyze their strengths and limitations.
- Discuss more recent approaches to signal peptide prediction.

The project therefore covers both:
1. **Signal peptide detection** — does the protein contain a signal peptide?
2. **Cleavage-site prediction** — where is the signal peptide cleaved?

---

## 🧪 Dataset

Protein sequences are collected from **UniProtKB** and restricted to **eukaryotic proteins**.

### Positive Dataset

Positive examples satisfy the project criteria, including:

- Reviewed UniProt entries
- No protein fragments
- Experimental evidence for a signal peptide
- Protein-level evidence for protein existence
- Protein length ≥ 40 amino acids
- Signal peptide length > 13 amino acids
- Known signal peptide cleavage site

### Negative Dataset

Negative examples satisfy:

- Reviewed UniProt entries
- No protein fragments
- Protein-level evidence for protein existence
- Protein length ≥ 40 amino acids
- No annotated signal peptide at any evidence level
- Experimental localization to selected non-secretory compartments

The datasets will be stored in both **TSV** and **FASTA** formats.

---

## 🗂️ Repository Structure

```text
LB2/
│
├── README.md
├── .gitignore
│
├── src/
│   └── collect_uniprot_data.py
│
├── data/
│   ├── raw/
│   └── processed/
│
├── notebooks/
│
├── results/
│
└── docs/
