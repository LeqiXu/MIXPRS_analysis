# MIXPRS_analysis  
This repository contains the codes used for **simulation studies and real-data analyses** in the **MIXPRS** paper.  

## Overview  
**MIXPRS** is a data fission-based multi-population PRS integration framework designed to effectively combine PRS derived from multiple populations and methods. The MIXPRS pipeline requires **only GWAS summary statistics and LD reference panels**.

For the **command-line tool** of **MIXPRS**, please visit: [https://github.com/LeqiXu/MIXPRS](https://github.com/LeqiXu/MIXPRS)

## Evaluation and Comparison  
The efficacy of MIXPRS is evaluated through **extensive simulations** and **real-data applications to 22 quantitative traits and four binary traits** across **five continental populations** (EUR, EAS, AFR, SAS, AMR), using data from **UK Biobank (UKBB) and All of Us (AoU)**.  

MIXPRS is asses under **three data scenarios**:  
- No tuning data available
- Tuning and testing data from the same cohort
- Tuning and testing data from different cohorts  

MIXPRS is compared with **seven state-of-the-art PRS methods**:  
- JointPRS
- SDPRX  
- XPASS  
- PRS-CSx  
- MUSSEL  
- PROSPER  
- BridgePRS

## Directory Structure  
- `Simulation/` – Contains scripts for **simulation studies** evaluating MIXPRS and seven exsiting methods under different settings.  
- `Real_data/` – Contains scripts for **real data analysis** evaluating MIXPRS and seven exsiting methods with 22 quantitative traits and 4 binary traits in UKBB and AoU.

## Support
Please direct any problems or questions to Leqi Xu (leqi.xu@yale.edu).
