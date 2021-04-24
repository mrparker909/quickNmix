
<!-- README.md is generated from README.Rmd. Please edit that file -->

# quickNmix

The goal of quickNmix is to aid in the fitting of asymptotic N-mixture
models, which are computed significantly faster than their canonical
counterpart when population sizes are large.

## Installation

You can install the released version of quickNmix from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("quickNmix")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mrparker909/quickNmix")
```

## Example

This is a basic example which shows how to fit a model with site varying
lambda, and time varying pdet:

``` r
library(quickNmix)

tictoc::tic()
nit = anmu[1:2,1:5] # ancient murrelet chick counts
mod = fitNmix(nit=nit,
              K=400, # upper bound on population size
              l_s_c=list(c(0,1)), # lambda site covariate
              p_t_c=list(c(0,1,1,1,1)),
              control=list(reltol=1e-5))
#> Warning: executing %dopar% sequentially: no parallel backend registered
tictoc::toc()
#> 564.72 sec elapsed
```

model AIC value:

``` r
mod$model_results$AIC
#> [1] 77.5027
```

lambda estimates for each site:

``` r
mod$model_results$estimate_matrices$lambda
#>          [,1]
#> [1,]  32.3202
#> [2,] 135.3340
```

gamma estimates for each site and time:

``` r
mod$model_results$estimate_matrices$gamma
#>          [,1]     [,2]     [,3]     [,4]     [,5]
#> [1,] 25.84412 25.84412 25.84412 25.84412 25.84412
#> [2,] 25.84412 25.84412 25.84412 25.84412 25.84412
```

omega estimates for each site and time:

``` r
mod$model_results$estimate_matrices$omega
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 0.9572106 0.9572106 0.9572106 0.9572106 0.9572106
#> [2,] 0.9572106 0.9572106 0.9572106 0.9572106 0.9572106
```

pdet estimates for each site and time:

``` r
mod$model_results$estimate_matrices$pdet
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 0.9901256 0.4421962 0.4421962 0.4421962 0.4421962
#> [2,] 0.9901256 0.4421962 0.4421962 0.4421962 0.4421962
```
