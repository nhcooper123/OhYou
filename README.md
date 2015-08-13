# A cautionary note on the use of Ornstein Uhlenbeck models in macroevolutionary studies.

## Natalie Cooper and Gavin H. Thomas. 
### Chris Venditti, Andrew Meade and Rob P. Freckleton.

Code, figures, etc. for the [paper](http://linktopaper) in the Biological Journal of the Linnaean Society special issue on "Radiations and Extinctions in Deep Time".

### Compiling the paper

To compile the paper you can use the Makefile:

```
make -C Manuscript
```
The Supporting Information file is available in PDF format [here](https://github.com/nhcooper123/OhYou/blob/master/Manuscript/ohyou_SuppInfo.pdf).
### Analysis code

Code needed to run the Maximum Likelihood simulations and make the figures in the paper can be found in the Analyses folder. 

* `OUsims.R` provides all R code for simulating trees and data (GHT).
* `OUfigures.R` provides all R code for making the simulations figures (Figure 4 and some removed in the revision and reworked as Tables 2-6) in the paper (GHT).
* `OhYouFigures.R` provides all R code for making the figures (1 and 3) and table for the literature review (NC).
* `plotFig2.R` provides all R code for making Figure 2 (RPF).
* `PagelParamPlots.R` provides all R code for making additional figures that were not used in the final version of the manuscript (GHT).

### Simulated data and literature review data

* `OU_simulations.rda` contains all the data from the OU simulations.
* `papers.txt` contains the number of papers found in the literature review.
* `literature.txt` contains a summary of the papers from the literature review.
* `ProfileDataRootyLam50.csv`, `ProfileDataRootyOU50.csv`, `ProfileDataTippyLam50.csv`, `ProfileDataTippyOU50.csv`, `ProfileDataYuleLam50.csv` and `ProfileDataYuleOU50.csv` contain Likelihood profile data used to construct Figure 2. 

### Useful(?) functions

At one point I (NC) began to rewrite the simulations code as a series of functions using the package [`diversitree`](https://github.com/richfitz/diversitree). These are provided in the R folder in the file `OhYou_simulations.R`. Feel free to play around with them!