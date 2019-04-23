---
title: "README"
author: "Victor Yuan"
date: "March 25, 2019"
output:
  html_document:
    keep_md: true
    toc: true
    toc_depth: 4
    toc_float:
      collapsed: false
    theme: spacelab
---

# Oxygen concentration and placental villous explant

For validating Beristain lab's oxygen project using public placental/decidual scRNAseq data.

## Directory structure

```
Beristain - Villous explants oxygenPreeclampsia score
|   README.html
│   README.md                   
|   README.Rmd 
│   2019-03-19 currated scRNA-seq gene list.xlsx  <- Alex/Jenna's list of hits to validate
|
└───data
│   └───main
|           └───Roser Vento-Tormo       <- The original data obtained from repo 
|           └───interim                 <- Intermediate data that has been transformed
│           └───processed               <- The final, canonical data sets for analysis
|
└───R
|   └───00_data                   <- Downloading & loading data into R
|   └───10_qc_preprocess          <- Quality control and preprocessing data
|   └───20_analysis               <- Analyzing processed data
|
└───outs                          <- Figures and tables
└───reports                       <- powerpoints and documents made to show to others
```

## Data

## Analysis

## References

