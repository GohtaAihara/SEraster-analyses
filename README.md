# SEraster-analyses

This repository contains code scripts to reproduce analyses and figures
in our *Bioinformatics* paper [Aihara G. et al. (2024), "SEraster: a
rasterization preprocessing framework for scalable spatial omics data
analysis",
*Bioinformatics*](https://doi.org/10.1093/bioinformatics/btae412), which
evaluates the utility of an R package SEraster.

## SEraster

`SEraster` is a rasterization preprocessing framework that aggregates
cellular information into spatial pixels to reduce resource requirements
for spatial omics data analysis. `SEraster` is available on
[GitHub](https://github.com/JEFworks-Lab/SEraster).

## Data

All of the datasets that were used in the manuscript are publicly
available datasets. References to original datasets used in the analyses
are provided in the manuscript.

## Contents

The code in this repository is organized in the
[analyses/](https://github.com/GohtaAihara/SEraster-analyses/tree/main/analyses)
directory, and file names indicate corresponding purposes.

### Preprocessing

All of the preprocessing was done with "preprocessing\_[spatial omics
technology]\_[SpeciesTissue].R" files:

Biological datasets (see Supplementary Information S3 in the manuscript
for more information):

-   [preprocessing_merfish_mouseBrain.R](https://github.com/GohtaAihara/SEraster-analyses/blob/main/analyses/preprocessing_merfish_mouseBrain.R)
-   [preprocessing_codex_humanIntestine.R](https://github.com/GohtaAihara/SEraster-analyses/blob/main/analyses/preprocessing_codex_humanIntestine.R)
-   [preprocessing_xenium_mousePup.R](https://github.com/GohtaAihara/SEraster-analyses/blob/main/analyses/preprocessing_xenium_mousePup.R)
-   [preprocessing_codex_humanSpleen.R](https://github.com/GohtaAihara/SEraster-analyses/blob/main/analyses/preprocessing_codex_humanSpleen.R)
    (only used in our [*bioRxiv*
    preprint](https://doi.org/10.1101/2024.02.01.578436))

Simulated datasets (see Supplementary Information S5 in the manuscript
for more information):

-   [preprocessing_nnsvg_simulations.R](https://github.com/GohtaAihara/SEraster-analyses/blob/main/analyses/preprocessing_nnsvg_simulations.R)
-   [preprocessing_crawdad_simulations.R](https://github.com/GohtaAihara/SEraster-analyses/blob/main/analyses/preprocessing_crawdad_simulations.R)

### Analyses

All of the analyses was done with "run\_[method]*[spatial omics
technology]*[SpeciesTissue].R" files:

-   [run_nnSVG_merfish_mouseBrain.R](https://github.com/GohtaAihara/SEraster-analyses/blob/main/analyses/run_nnSVG_merfish_mouseBrain.R):
    This file produced Figure 2A-G; Supplementary Figure 1;
    Supplementary Figure 1.
-   [run_nnSVG_nnSVG_simulations.R](https://github.com/GohtaAihara/SEraster-analyses/blob/main/analyses/run_nnSVG_nnSVG_simulations.R):
    This file produced Figure 2H, I.
-   [run_nnSVG_xenium_mousePup.R](https://github.com/GohtaAihara/SEraster-analyses/blob/main/analyses/run_nnSVG_xenium_mousePup.R):
    This file produced Figure 1; Figure 4B-C.
-   [run_CooccurrenceAffinity_codex_humanIntestine.R](https://github.com/GohtaAihara/SEraster-analyses/blob/main/analyses/run_CooccurrenceAffinity_codex_humanIntestine.R):
    This file produced Figure 3A-F; Supplementary Figure 3.
-   [run_CooccurrenceAffinity_xenium_mousePup.R](https://github.com/GohtaAihara/SEraster-analyses/blob/main/analyses/run_CooccurrenceAffinity_xenium_mousePup.R):
    This file produced Figure 4A, D-F; Supplementary Figure 4.
-   [run_CooccurrenceAffinity_crawdad_simulation.R](https://github.com/GohtaAihara/SEraster-analyses/blob/main/analyses/run_CooccurrenceAffinity_crawdad_simulation.R):
    This file produced Figure 3G, H.
-   [run_CooccurrenceAffinity_codex_humanSpleen.R](https://github.com/GohtaAihara/SEraster-analyses/blob/main/analyses/run_CooccurrenceAffinity_codex_humanSpleen.R):
    This file produced Figure 3A-B in our *bioRxiv* preprint.
-   [run_SOMDE_merfish_mouseBrain.ipynb](https://github.com/GohtaAihara/SEraster-analyses/blob/main/analyses/run_SOMDE_merfish_mouseBrain.ipynb):
    This file applied `SOMDE` (Python package) to downsample the
    preprocessed MERFISH mouse brain dataset. The downsampled datasets
    were used in
    [run_nnSVG_merfish_mouseBrain.R](https://github.com/GohtaAihara/SEraster-analyses/blob/main/analyses/run_nnSVG_merfish_mouseBrain.R)
    to produce Figure 2E.

Others:

-   [functions.R](https://github.com/GohtaAihara/SEraster-analyses/blob/main/analyses/functions.R):
    This file contains utility functions that were used in other files
    for convenience.
    -   `calculateDensity`: computes the density of a matrix.
    -   `rotateAroundCenter`: rotates x,y coordinates at a given angle
        around a midrange point. Later implemented in `SEraster` as
        `permutateByRotation`.
    -   `calculatePerformanceMetrics`: computes performance metrics
        (e.g. TPR, PPV, TNR) based on given ground truth and predicted
        labels.
    -   `gg_color_hue`: splits the ggplot2 standard rainbow color
        palette to a given integer.
