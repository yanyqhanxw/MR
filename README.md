# MR
This repository contains code and data to conduct MR analysis.
# Pipeline
- Select instruments
    - Significantly associated with exposure: P < 5 × 10E-8
    - PLINK clump: r2 > 0.01, window < 10kb
    - Exclude SNPs directly associated with outcome: P < 5 × 10E-8
    - Harmonize allele
    - RadialMR: exclude outlier pleiotropic SNPs
    - calculate F statistics and exclude SNPs with F < 10
- MR analysis
    - Inverse variance weighted (IVW)
    - Weighted median (WM)
    - MR Egger
    - MRPRESSO
    - reverse causal effect test: MR Steiger
- sensitivity analysis
    - pleiotropy
    - heterogeneity
# Requirements
- R_4.1.1
    - RadialMR_1.0
    - TwoSampleMR_0.5.6
    - MRPRESSO_1.0
