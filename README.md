# thesis-random-polynomials
This is the companion repository for my Master thesis "Statistical Analysis of Roots of Random Polynomials".
The repository contains all the MATLAB files used to perform the simulations described in Chapter 5 of the thesis.

## Content
- The *Resources* folder contains a set of functions providing useful utilities for the simulations. Each function contains a description of its usage. Some functions have been created by the outhor of this repository, some others have been adapted from external sources. In that case the original source is credited in the function description and cited in the thesis.
- Simulation files, named with the format *DegreeN_NoiseXxxxx_SimulYyyyyyy*, where:
 -- *N*: Degree of the polynomial
 -- *Xxxxx*: Circular, (Diagonal,) FullMatrix
 -- *Yyyyyyy*: a short description of the simulation (e.g. BiasVsIter, MseVsSomething, ...)
- Logistic regression analyisis files, named with the format *DegreeN_LogisticRegressionXxxxx* where *Xxxxx* indicates the type of statistical test result (e.g. Henze-Zirkler, Hotelling's T2, ...) that is analyzed with the logistic regression procedure described in the thesis.

## Thesis experiments and repository files
The following table shows which repository files are associated with which Experiment in the thesis. Files not mentioned in the table have been used for a preliminary analysis and are not directly associated to any experiment in the thesis.

| Plugin | README |
| ------ | ------ |
| Dropbox | [plugins/dropbox/README.md][PlDb] |
| GitHub | [plugins/github/README.md][PlGh] |
| Google Drive | [plugins/googledrive/README.md][PlGd] |
| OneDrive | [plugins/onedrive/README.md][PlOd] |
| Medium | [plugins/medium/README.md][PlMe] |
| Google Analytics | [plugins/googleanalytics/README.md][PlGa] |


