# Mental Arithmetic EEG Signal Analysis

## Overview

This project analyzes **EEG signals recorded during repeated mental subtraction tasks** to investigate cognitive workload, brain region activation, and functional connectivity.  

Using **MATLAB**, we processed and visualized EEG data from six subjects, employing **power spectral density (PSD)** and **magnitude-squared coherence (MSC)** to assess brain activity in theta (θ), alpha (α), and beta (β) frequency bands.  

## Features

- **Baseline-normalized PSD:** Highlights cognitive task activation relative to rest.  
- **Frequency Bands Analyzed:** δ (1–4 Hz), θ (4–8 Hz), α (8–13 Hz), β1 (13–20 Hz), β2 (20–30 Hz), γ (30–40 Hz).  
- **Hemispherical Analysis:** Observed left hemisphere dominance in task engagement.  
- **Functional Connectivity:** Evaluated using magnitude-squared coherence, revealing strong frontal and fronto-parietal interactions.  
- **Artifact Handling:** Removed excessive noise and detrended signals to ensure stability.  
- **Visualization:** Topographical maps, boxplots, and coherence maps for intuitive understanding of cognitive workload.  

## Key Findings

- Increased theta, alpha, and beta power in left frontal and parieto-occipital regions during mental arithmetic.  
- Left hemisphere dominance during familiar, practiced tasks.  
- Significant functional connectivity observed in frontal regions, suggesting coordinated neural processing.  

## Methodology

1. **Data Collection:** EEG with 10–20 electrode placement; 500 Hz sampling; high-pass (0.5 Hz), low-pass (45 Hz), and notch (50 Hz) filtering applied.  
2. **Preprocessing:** Artifact removal, detrending, and selection of stable 45-second windows.  
3. **Analysis:**  
   - PSD computed via Welch method (Hamming window, 10 s, 0.1 s overlap).  
   - MSC computed for all electrode pairs to assess connectivity.  
   - Band power integrated per frequency band.  
4. **Visualization:** MATLAB vectorial plots for topography, distribution, and coherence.  
