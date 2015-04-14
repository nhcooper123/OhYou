# A cautionary note on the use of Ornstein Uhlenbeck and other models in macroevolutionary studies.

## Natalie Cooper and Gavin H. Thomas. 
### plus Chris Venditti Andrew Meade and Rob P. Freckleton.

Code, figures, etc. for the [paper](http://linktopaper) in Biological Journal of the Linnaean Society.

### Compiling the paper

To compile the paper you can use the Makefile:

```
make -C Manuscript
```
### Analysis code

Code needed to run the Maximum Likelihood simulations and make the figures in the paper can be found in the Analyses folder. 

* `OUsims.R` provides all R code for simulating trees and data (GHT).
* `OUfigures.R` provides all R code for making the simulations figures in the paper (GHT).
* `OhYouFigures.R` provides all R code for making the figures and table for the literature review (NC).
* `PagelParamPlots.R` provides all R code for making figures for the supplementary material (GHT).

### Simulated data and literature review data

* `OU_simulations.rda` contains all the data from the OU simulations.
* `papers.txt` contains the number of papers found in the literature review.
* `literature.txt` contains a summary of the papers from the literature review.

### Useful(?) functions

At one point I (NC) began to rewrite the simulations code as a series of functions using the package [`diversitree`](https://github.com/richfitz/diversitree). These are provided in the R folder in the file `OhYou_simulations.R`. Feel free to play around with them!