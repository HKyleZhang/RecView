[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# RecView<br><sub>- *view* and *locate* recombinations positions using pedigree data -</sub>

This R package is designed to distribute the *RecView* ShinyApp which aims at providing a user-friendly GUI for viewing and locating recombination positions on chromosomes using pedigree data.

<br>

## **Installation**

- Developmental `devtools::install_github("HKyleZhang/RecView@dev")`
- Latest (v1.1.0): `devtools::install_github("HKyleZhang/RecView@v1.1.0")`

Click [here](https://github.com/HKyleZhang/RecView) to see all the details about the package.

<br>

## **Changelog**

**Nov. 21, 2024**

- to be added

**Nov. 3, 2024**

- Enable inferring grandparent-of-origin when genotypes of some individuals are missing at all or some sites.
- Enable preview when multiple offspring and chromosomes are selected for analysis.
- Show number of informative sites in GoO figure.
- Reduce RAM usage by changing the way of loading input files.
- Reduce running time for the PD algorithm.
