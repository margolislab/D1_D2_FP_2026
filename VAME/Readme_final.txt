Margolis Lab: SNI Longitudinal Behavior & Photometry Pipeline
This repository contains the data processing, alignment, and statistical modeling pipeline for integrating TDT fiber photometry data with automated behavioral tracking (SimBA and VAME). It is specifically designed to analyze D1 and A2A neuronal dynamics in the dorsolateral striatum during naturalistic behaviors in Spared Nerve Injury (SNI) and Sham rodent models.

🛠 Prerequisites & Dependencies
MATLAB: Requires the Signal Processing Toolbox and Statistics and Machine Learning Toolbox.

TDT SDK: The MATLAB TDT bin2mat SDK must be installed and added to your MATLAB path.

Python: Required if running the VAME behavioral segmentation pipeline locally.

SimBA: Behavioral classifications should be pre-processed and exported as .csv or .xlsx files containing bout Start_time and End_time.

📂 Repository Structure
Plaintext
├── run_pipeline.py                     # Python orchestrator for the VAME workflow
├── step1_photometry_batch.m            # Shared TDT photometry cleaning & Z-scoring
├── step2_vame_photometry_alignment.m   # VAME motif synchronization and PETH extraction
├── step3_vame_auroc_permutation.m      # VAME 10k permutation testing
├── step4_vame_stats_recovery.m         # VAME Hedges' g & KS-test statistics
└── master_simba_photometry.m           # Unified SimBA alignment, PSTH, and stats pipeline
🚀 Workflow Instructions
Phase 1: Neural Data Preprocessing (Shared)
Script: step1_photometry_batch.m

Function: Serves as the universal first step for both VAME and SimBA workflows.

Process: Imports raw TDT streams, synchronizes them to Blackbox video TTL pulses (PC0_) to remove initialization surges, performs linear detrending, and calculates robust Median/MAD-based Z-scores to isolate movement-corrected signals (zGCaMP - zIso).

Output: Centralized .mat files containing cleaned_Z and exact timestamps (export_data.timestamps).

Phase 2: Behavioral Alignment & Statistics
Choose your behavioral tracking pathway below:

Pathway A: SimBA Alignment
Script: master_simba_photometry.m

Process: Prompts for your SimBA summary Excel file and the folder containing the Phase 1 .mat files. It dynamically filters for specified behaviors (e.g., locomotion, grooming, rearing), genotypes, and timepoints (e.g., Baseline 2 vs. Week 12).

Analysis: Extracts bout-level snippets using absolute time (Start_time), performs 10,000-shuffle permutation tests for auROC calculation, and downsamples data (~100Hz) for clean vector plotting.

Outputs:

Illustrator-ready vector PSTHs (.pdf and .svg).

auROC distribution histograms.

GraphPad Prism-ready histogram data tables.

Console printout of manuscript-ready statistics: Kolmogorov-Smirnov (KS) tests, Variance F-tests, and Hedges' g effect sizes (corrected for small sample sizes).

Pathway B: VAME Alignment
Scripts: run_pipeline.py -> step2 through step4 MATLAB scripts.

Process: Run run_pipeline.py to initiate the VAME 15-state HMM segmentation natively in Python. The script will then automatically trigger the downstream MATLAB scripts.

Analysis: Syncs VAME community onset frames with the cleaned photometry Z-scores using a sampling rate ratio conversion (video FPS to TDT Hz). Extracts peri-event time histograms and runs permutation testing on bout-level metrics (Peak Z-score, AUC).

Outputs: PDF/SVG PETH overlays, side-by-side metric boxplots, bout duration correlation plots, and manuscript-ready statistical reports.

📊 Statistical Note
To account for longitudinal experimental designs, the effect size for distribution shifts (e.g., Baseline vs. Week 12) is calculated using Hedges' g, which includes a correction factor for small sample sizes. Significance testing for state-dependent calcium activity relies on two-tailed permutation testing against 10,000 temporally shifted null distributions.