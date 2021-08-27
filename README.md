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
The following table shows which repository files are associated with which Experiment in the thesis. Files not mentioned in the table are not directly associated to any experiment in the thesis.

| Thesis Experiment number | Files in the repository |
| ------ | ------ |
| Experiment 5.1 | Degree02\_NoiseCircular\_SimulBiasVsIterOne.m |
| Experiment 5.2 | Degree02\_NoiseCircular\_SimulMseVsSnr.m |
| Experiment 5.3 | Degree02\_NoiseCircular\_DistanceSNR.m |
| Experiment 5.4 | Degree02_NoiseCircular_ProjectionAndOthers_MULTIPLETIMES.m <br /> Degree02_LogisticRegressionGaussHZ.m <br /> Degree02_LogisticRegressionHotT2.m |
| Experiment 5.5 | Degree02_NoiseCircular_SimulManova.m <br /> Degree02_LogisticRegressionManovaGaussHZ.m <br /> Degree02_LogisticRegressionManovaMBox.m <br /> Degree02_LogisticRegressionManovaMean.m |
| Experiment 5.6 | Degree03_NoiseCircular_ProjectionAndOthers_MULTIPLETIMES.m <br /> Degree03_LogisticRegressionGaussHZ.m <br /> Degree03_LogisticRegressionHotT2.m |
| Experiment 5.7 | Degree03_NoiseCircular_SimulManova.m <br /> Degree03_LogisticRegressionManovaGaussHZ.m <br /> Degree03_LogisticRegressionManovaMBox.m <br /> Degree03_LogisticRegressionManovaMean.m |
| Experiment 5.8 | Degree02\_NoiseFullMatrix\_SimulBiasVsIterOne.m |
| Experiment 5.9 | Degree02_NoiseFullMatrix_ProjectionAndOthers_MULTIPLETIMES.m <br /> Degree02_LogisticRegressionGaussHZ.m <br /> Degree02_LogisticRegressionHotT2.m |
| Experiment 5.10 | Degree02_NoiseFullMatrix_SimulManova.m <br /> Degree02_LogisticRegressionManovaGaussHZ.m <br /> Degree02_LogisticRegressionManovaMBox.m <br /> Degree02_LogisticRegressionManovaMean.m |


