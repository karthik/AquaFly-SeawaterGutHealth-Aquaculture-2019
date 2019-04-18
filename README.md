# AquaFly-SeawaterGutHealth-Aquaculture-2018
This repository contains the data and code to reproduce the results presented in the paper **Gut function and health in post-smolt Atlantic salmon (*Salmo salar*) fed larvae meal from black soldier fly (*Hermetia illucens*)**.

To run the analyses, download this repository as a zipped file. After decompression, open the R project file (AqFl2_GutHealth.Rproj) in the RStudio and run the R codes directly from the project root directory. Versions of R packages for each anaylsis can be found in the session information (*_sessionInfo.txt) under the same folder of R codes.

Below is the file organization in this repository.
```
├── analysis
│   ├── code
│   │   ├── histology.R
│   │   ├── histology_sessionInfo.txt
│   │   ├── organosomatic_index.R
│   │   ├── organosomatic_index_sessionInfo.txt
│   │   ├── qPCR_ref_DI.R
│   │   ├── qPCR_ref_DI_sessionInfo.txt
│   │   ├── qPCR_ref_PI.R
│   │   ├── qPCR_ref_PI_sessionInfo.txt
│   │   ├── qPCR_target_DI.R
│   │   ├── qPCR_target_DI_sessionInfo.txt
│   │   ├── qPCR_target_PI.R
│   │   └── qPCR_target_PI_sessionInfo.txt
│   ├── exploratory_analysis
│   │   ├── histology_histogram.pdf
│   │   ├── OSI_boxPlot.pdf
│   │   ├── OSI_violin.pdf
│   │   ├── qPCR_ref_DI_barPlot.pdf
│   │   ├── qPCR_ref_DI_boxPlot.pdf
│   │   ├── qPCR_ref_DI_pointDiagram.pdf
│   │   ├── qPCR_ref_PI_barPlot.pdf
│   │   ├── qPCR_ref_PI_boxPlot.pdf
│   │   ├── qPCR_ref_PI_pointDiagram.pdf
│   │   ├── qPCR_target_DI_boxPlot.pdf
│   │   ├── qPCR_target_DI_heatmap.pdf
│   │   ├── qPCR_target_DI_violin.pdf
│   │   ├── qPCR_target_PI_boxPlot.pdf
│   │   ├── qPCR_target_PI_heatmap.pdf
│   │   └── qPCR_target_PI_violin.pdf
│   └── model_diagnostics
│       ├── OSI_LMMs_residual.pdf
│       ├── OSI_welch-t_qqplot.pdf
│       ├── qPCR_DI_LMMs_residual.pdf
│       ├── qPCR_DI_welch-t_qqplot.pdf
│       ├── qPCR_PI_LMMs_residual.pdf
│       └── qPCR_PI_welch-t_qqplot.pdf
├── AqFl2_GutHealth.Rproj
├── data
│   ├── clean_data
│   │   ├── AqFl2_qPCR_ref_DI.csv
│   │   ├── AqFl2_qPCR_ref_PI.csv
│   │   ├── AqFl2_qPCR_Target_DI.csv
│   │   ├── AqFl2_qPCR_Target_PI.csv
│   │   └── Interplate_calibration.xlsx
│   ├── raw_data
│   │   ├── AqFl2_histology.csv
│   │   ├── AqFl2_LightCycler_files
│   │   │   ├── DI
│   │   │   │   ├── ref_gene
│   │   │   │   │   ├── AqFl2_DI_actb_ann temp 60.lc96p
│   │   │   │   │   ├── AqFl2_DI_gapdh_ann temp 60.lc96p
│   │   │   │   │   ├── AqFl2_DI_hprt1_ann temp 60.lc96p
│   │   │   │   │   └── AqFl2_DI_rnapo2_ann temp 60.lc96p
│   │   │   │   └── target_gene
│   │   │   │       ├── AqFl2_DI_apoa1_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_apoa4_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_apob_ann temp 63.lc96p
│   │   │   │       ├── AqFl2_DI_aqp8ab_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_casp6_ann temp 60_.lc96p
│   │   │   │       ├── AqFl2_DI_cat_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_cd3gd_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_cd8b_ann temp 59.lc96p
│   │   │   │       ├── AqFl2_DI_chk_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_cldn15_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_cldn25b_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_cyp1a1_ann temp 55.lc96p
│   │   │   │       ├── AqFl2_DI_Ecad_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_fabp2b_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_foxp3_ann temp 62.lc96p
│   │   │   │       ├── AqFl2_DI_hsp70_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_ifng_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_il10_ann temp 62.lc96p
│   │   │   │       ├── AqFl2_DI_il17a_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_il1b_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_il4_ann temp 62.lc96p
│   │   │   │       ├── AqFl2_DI_il6_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_il8_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_mhcI_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_mmp13_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_mta_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_mtp_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_muc2_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_myd88_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_pcna_ann temp 55.lc96p
│   │   │   │       ├── AqFl2_DI_pcyt1a_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_plin2_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_sod1_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_tgfb1_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_DI_tnfa_ann temp 60.lc96p
│   │   │   │       └── AqFl2_DI_zo1_ann temp 60.lc96p
│   │   │   ├── PI
│   │   │   │   ├── ref_gene
│   │   │   │   │   ├── AqFl2_PI60_RefGenes_ann temp 60.lc96p
│   │   │   │   │   ├── AqFl2_PI_actb_ann temp 60.lc96p
│   │   │   │   │   ├── AqFl2_PI_gapdh_ann temp 60.lc96p
│   │   │   │   │   ├── AqFl2_PI_hprt1_ann temp 60.lc96p
│   │   │   │   │   └── AqFl2_PI_rnapo2_ann temp 60.lc96p
│   │   │   │   └── target_gene
│   │   │   │       ├── AqFl2_PI_apoa1_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_apoa4_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_apob_ann temp 63.lc96p
│   │   │   │       ├── AqFl2_PI_aqp8ab_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_casp6_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_cat_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_cd3gd_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_cd8b_ann temp 59.lc96p
│   │   │   │       ├── AqFl2_PI_cdh1_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_chk_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_cldn15_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_cldn25b_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_cyp1a1_ann temp 55.lc96p
│   │   │   │       ├── AqFl2_PI_fabp2b_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_foxp3_ann temp 62.lc96p
│   │   │   │       ├── AqFl2_PI_hsp70_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_ifng_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_il10_ann temp 62.lc96p
│   │   │   │       ├── AqFl2_PI_il17a_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_il1b_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_il4_ann temp 62.lc96p
│   │   │   │       ├── AqFl2_PI_il6_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_il8_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_mhcI_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_mmp13_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_mta_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_mtp_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_muc2_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_myd88_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_pcna_ann temp 55.lc96p
│   │   │   │       ├── AqFl2_PI_pcyt1a_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_plin2_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_sod1_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_tgfb1_ann temp 60.lc96p
│   │   │   │       ├── AqFl2_PI_tnfa_ann temp 60.lc96p
│   │   │   │       └── AqFl2_PI_zo1_ann temp 60.lc96p
│   │   │   └── README.md
│   │   └── AqFl2_organosomatic_index.csv
│   └── README.md
├── LICENSE
├── README.md
└── results
    ├── figures
    │   ├── Figure 1.tiff
    │   ├── Figure 2.tiff
    │   ├── Figure 3.tiff
    │   └── Figure 4.tiff
    └── reference_gene_ranks
        ├── ref_gene_rank_DI.pdf
        └── ref_gene_rank_PI.pdf
```
