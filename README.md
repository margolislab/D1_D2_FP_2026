# D1_D2_FP_2026
Analysis of fiber photometry for evoked, supervised, and unsupervised learning
# ARM Mechanostimulation: Photometry TTL Alignment Pipeline

This script processes TDT fiber photometry data aligned strictly to hardware-generated TTL pulses from Automated Mechanostimulation (ARM) equipment. 

Unlike the VAME and SimBA pipelines—which synchronize neural traces to visually decoded behavioral state timestamps—this pipeline relies on exact millisecond-precision event markers (`PC1_` epocs) to calculate evoked stimulus responses, sliding-window AUROC discriminability, and calcium decay kinetics.

## 🛠 Prerequisites
* **MATLAB**: Requires the Curve Fitting Toolbox (for exponential decay tracking) and the Signal Processing Toolbox.
* **TDT SDK**: The MATLAB TDT bin2mat SDK must be installed.
* **AUC Functions**: Requires the custom `auc.m` and `auc_bootstrap.m` functions in your MATLAB path for the sliding Hanley AUROC calculations.

## ⚙️ Processing Workflow

1. **Signal Detrending & Cleaning**
   * Bypasses the first 8000 samples to strip initialization noise and warm-up transients.
   * Calculates a 2nd-degree polynomial line of best fit to remove slow baseline drift.
   * Performs robust Median/MAD-based Z-scoring to isolate true GCaMP transients from the Isosbestic control.
   * Replaces statistical artifacts (>4 SD or <-2 SD) with linear interpolation and a 5-point median filter.

2. **Event Alignment**
   * Extracts ARM TTLs (`PC1_` onset).
   * Generates a 20-second per-event window (`-10s` to `+10s`) for the stimulus response, and a separate 10-second spontaneous activity window.
   * Compiles event-aligned rasters and per-animal mean traces (`Animal_Traces.mat`).

3. **Sliding Window AUROC**
   * Slides a 200ms integration window (at 50ms steps) across the post-stimulus timeline.
   * Compares each post-stimulus window against a 1-second pre-stimulus baseline.
   * Extracts maximum discriminability (Peak AUROC), latency to significance (via bootstrapping), and total AUC integral.

4. **Calcium Kinetics & Decay Analysis**
   * Isolates the trace following the peak post-stimulus response.
   * Fits a single-exponential decay curve (`A*exp(-x/tau)+C`) to calculate tau ($\tau$) and half-decay latency.
   * Calculates a highly robust **90%-70% decay slope** (avoiding late-trace plateau noise).

## 📊 Outputs

Running this script generates several files in the selected TDT data folder:

* `Animal_Traces.mat`: Contains the raw matrix of curated trials, the mean trace, and the corresponding time vector.
* `ZAlign_summary.svg`: A 4-panel plot showing event vs. spontaneous rasters and mean traces.
* `AUROC_plot.svg`: A plot of AUROC discriminability over time.
* `Decay_Fit.svg`: A visual validation of the exponential decay curve fit and baseline estimation.
* `AUROC_metrics.csv`: Peak AUROC, peak time, latency to significance, and duration above significance threshold.
* `ZScore_CalciumTrace_Metrics.csv`: Traditional trace metrics including Max Amplitude, Peak Latency, Decay $\tau$, Half-decay, and the robust 90-70% decay slope.
