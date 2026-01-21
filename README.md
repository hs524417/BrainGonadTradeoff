# BrainGonadTradeoff

# Overview

This repository contains the data and analysis code for the manuscript "Male reproductive tactics drive correlated selection for large gonads with small brains." This study examines the relationship between brain and gonad size in _Poecilia parae_ and related species to determine whether negative correlations arise from energetic trade-offs or correlated selection. The findings suggest that reproductive strategies shape brain evolution through correlated selection rather than direct resource allocation constraints, challenging traditional assumptions about sexual selection’s impact on brain size.

# Repository Contents

TissueInvestment_R_Code.R: R script for performing statistical analyses on the dataset.

TissueInvestmentData.csv: Contains raw data used in the statsical analysis, including:
  Body size (length and weight), Brain weight, Male Gonad weight, Neuron/glia ratio, Liver Weight, Intestine Length, and Brain Cell Count 

Neuron_Glia_Ratio.xlsx: Excel file used to compute the neuron/glia ratio. The Neuron/glia ratio is calculated using the number of cells stained with GFP and DAPI divided by the number of cells that only stained DAPI.

Brain_Soup_Cell_Counts.xlsx: Excel file used to calculate total number of brain cells. This file includes:
  Raw cell counts (Number of DAPI stained cells across 10 aliquots), Total number of cells (unfiltered), Filtered Total Cells (CV<0.15, and CV<0.10) (CV= coefficient of vairation)
  

# Citation

If you use this data in your research, please cite : Stec, H., Zhang, G.Y., Sandkam, B.A. (2025). Male reproductive tactics drive correlated selection for large gonads with small brains. _biorxiv_. https://doi.org/10.1101/2025.03.22.644732

