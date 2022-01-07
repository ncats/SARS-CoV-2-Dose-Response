# SARS-CoV-2 Dose Response in Cynomolgous Macaques

Our previous work demonstrated that [Seroconversion and fever are dose-dependent in a nonhuman primate model of inhalational COVID-19](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1009865), based on a dose response experiment of the WA-1. The dose response experiment was repeated with the Gamma variant of SARS-CoV-2 and the results were compared to those from WA-1. This repository contains the code and rough interpretation of the analysis of the Gamma vs WA-1 dose response results.


## Approach

In a previous publication we performed bivariate modeling generalized linear modeling for the fever and seroconversion endpoints simultaneously for just WA-1.  Here we repeat that analysis for the Gamma variant and do a comparison between the two variants for both endpoints. To provide clarity into the effect that modeling choices have on our results, I estimate equivalent parameters and perform similar tests on the effect of the variant using a standard univariate GLM approach (as opposed to the bivariate model). In addition to seroconversion and fever, I explore the effect of the variant on shedding rate.

See the [RMarkdown document](https://github.com/ncats/SARS-CoV-2-Dose-Response/blob/main/Seroconversion-and-Fever-Dose-Response-for-Gamma-and-WA-1-Variants.md) for the complete analysis.

## Dependencies

### Statistical analysis

* [VGAM](https://cran.r-project.org/web/packages/VGAM/index.html)
* [emdbook](https://cran.r-project.org/web/packages/emdbook/index.html)

### Processing and data viz

* [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html)
* [scales](https://cran.r-project.org/web/packages/scales/index.html)
* [cowplot](https://cran.r-project.org/web/packages/cowplot/index.html)
