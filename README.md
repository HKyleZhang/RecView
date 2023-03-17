[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# RecView<br><sub>- _view_ and _locate_ recombinations positions using pedigree data -</sub>

This R package is designed to distribute the _RecView_ ShinyApp which aims at providing a user-friendly GUI for viewing and locating recombination positions on chromosomes using pedigree data.

<br>

![](https://github.com/HKyleZhang/Thesis_Figure_and_Supplementary/blob/main/3.Paper-III_RecView/Figure/Figure_06.png)

<br>

## Installation

`remotes::install_github("HKyleZhang/RecView")`

+ Installation will install these dependencies: 
  * shiny (>= 1.7.1)
  * shinyjs (>= 2.1.0)
  * shinycssloaders
  * foreach
  * doParallel 
  * parallel,
  * tidyverse (>= 1.3.1)
  * scattermore (>= 0.8) 
  * ggpubr, 
  * ggpmisc, DT, 
  * vcfR 
  * colorspace

<br>

## List of functions

Function  | Description
----------|------------
`make_012gt()` | Formats the genotype file for _RecView_.
`make_012gt_from_vcf()` | Formats the genotype file from VCF file for _RecView_.
`run_RecView_App()` | Invokes _RecView_. _RecView_ provides options to save the result figures and tables to your current working directory.

<br>

## More details about the _RecView_ ShinyApp
### Required input files

File          | Description
--------------|-------------
Genotype file | One .csv file having the genotypes in "012" genotype format with column-wise organisation of the individual/sample genotypes of variants positions. It should have the following columns: `A` is paternal father; `B` is paternal mother; `C` is maternal father; `D` is maternal mother; `AB` is father; `CD` is mother, and the column(s) for offspring. This file can be generated by using `make_012gt()` (or `make_012gt_from_vcf()`).
Scaffold file | One .csv file having the order and orientation of the reference genome scaffolds. It should have the following columns (names are case sensitive): `scaffold`, `size`,	`CHR`, `order`, `orientation`. Note: with a chromosome-level assembly, this file can be tweaked so to make `scaffold` and `CHR` identical, but still keep separate columns.

<br>

### Additional settings
* __Choose offspring(s)__: Choose the offspring for the analysis. It supports multiple selection.

* __Choose chromosome(s)__: Choose the chromosome for the analysis. It supports multiple selection.

* __Locate recombination positions?__ Check 'Yes' to locate recombination positions with either of two algorithms (see below).

* __Algorithms (optional):__ 
  + PD: Proportional Difference algorithm proceeds by specifying a window size (the number of informative SNPs of each flanking window), a step value (k) giving the number of SNPs between each calculated position, and a threshold to trigger denser calculations (at every SNP) to detect local maxima.
  + CCS: Cumulative Continuity Score algorithm calculates a CCS for each position along the chromosome, and (ii) finds putative recombination positions by locating regions where long continuously increasing slopes of CCSs of one grandparent-of-origin is replaced by long continuously increasing slopes of CCSs from the other grandparent. 

* __Radius value (PD optional)__: the number of informative SNPs around the examined position for calculating the proportion of informative SNPs from specific grandparents.

* __Step value (PD optional)__: the step size to move along the chromosome. Larger values decrease the number of positions to be examined, while increasing analysis speed.

* __Finer step value (PD optional)__: the step size to move along the chromosome, after the absolute difference of the proportion of grandparent-of-origin reaches above the threshold. Larger value decreases the positions to be examined, while increasing the analysis speed.

* __Threshold (PD optional):__ the condition to initiate a finer step, and later filter the local maxima for _effectively true_ recombination.

* __Threshold (CCS optional):__ the minimal CCS to consider an _effectively true_ recombination. Larger value is more stringent and captures crossovers, while small value captures both crossovers and non-crossovers. However, small values can also capture artefacts of recombination due to wrongly called genotypes.

* __Saving options (optional):__ 
  + GoO Inference: this option will save inferences of grandparent-of-origin for the selected offspring(s) as csv-file(s) separately for each selected chromosome in the current working directory. 
  + Plots: this option will save the result figures for the selected offspring(s) as pdf-file(s) separately for each selected chromosome in the current working directory.
  + Locations: when _Locate recombination positions?_ is checked "Yes", this option will save the table of the putative recombination locations in the selected offspring(s) as csv-file(s) separately for each selected chromosome in the current working directory.

* __Run analysis__ button: start the analysis!

<br>

## Example workflow
For big VCF file, it is recommended to continue with __Workflow A__.

#### __Workflow A__
  1. Use `--extract-FORMAT-info GT` option in VCFtools to extract genotypes into a single file.
  2. Use `make_012gt()` to format the genotype file. 
  3. prepare scaffold file.
  4. In Rstudio, navigate to the working directory to where you want to save the result figures and tables.
  5. In Rstudio, start the _RecView_ ShinyApp by `run_RecView_App()`; continue with settings and run analysis.

#### __Workflow B__
  1.  Use `make_012gt_from_vcf()` to format the genotype file directly from VCF file. 
  2. prepare scaffold file.
  3. In Rstudio, navigate to the working directory to where you want to save the result figures and tables.
  4. In Rstudio, start the _RecView_ ShinyApp by `run_RecView_App()`; continue with settings and run analysis.
  
<br>

## Please cite:

Zhang H, Hansson B. RecView: an interactive R application for viewing and locating recombination positions using pedigree data. BioRxiv 2022:2022.12.21.521365. doi:10.1101/2022.12.21.521365.

<br>

## Version history
+ version 1.0.0: Dec, 2022
