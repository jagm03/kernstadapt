# kernstadapt

**kernstadapt** is an R package for adaptive kernel estimation of the intensity of spatio-temporal point processes.

**kernstadapt** implements functionalities to estimate the intensity of a spatio-temporal point pattern by kernel smoothing with adaptive bandwidth methodology when each data point has its own bandwidth associated as a function of the crowdedness of the region (in space and time) in which the point is observed.

The package presents the intensity estimation through a direct estimator and the partitioning algorithm methodology presented in [González and Moraga
(2022)](https://arxiv.org/pdf/2208.12026.pdf).


## Installation

The stable version on CRAN can be installed using:

```{r, eval=FALSE}
install.packages("kernstadapt")
```

The development version can be installed using **devtools**:

```{r, eval=FALSE}
# install.packages("devtools") # if not already installed
devtools::install_github("jagm03/kernstadapt")
library(kernstadapt)
```

## Main functions

Direct adaptive estimation of the intensity

- `dens.direct()` (non-separable)
- `dens.direct.sep()` (separable)


Adaptive intensity estimation using a partition algorithm

- `dens.par()` (non-separable)
- `dens.par.sep()` (separable)

Bandwidths calculation

- `bw.abram.temp()` (temporal)

Separability test

- `separability.test()`

