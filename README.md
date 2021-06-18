
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

## Example 1

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

## Example 2

This is a basic example which shows how you can use multiple cores to
compute:

``` r
library(quickNmix)

# uses library doParallel, here we use 2 cores
doParallel::registerDoParallel(cores = 2) 

tictoc::tic()
nit  = anmu[c(2,5),1:12] # ancient murrelet chick counts
mod2 = fitNmix(nit=nit, K=400)
tictoc::toc()
#> 2908.69 sec elapsed

mod2
#> $optim_results
#> $optim_results$par
#>     B_l_0     B_g_0     B_o_0     B_p_0 
#>  4.995500  4.279050 -2.047628  3.029386 
#> 
#> $optim_results$value
#> [1] 107.5519
#> 
#> $optim_results$counts
#> function gradient 
#>      231      100 
#> 
#> $optim_results$convergence
#> [1] 1
#> 
#> $optim_results$message
#> NULL
#> 
#> 
#> $model_results
#> $model_results$NLL
#> [1] 107.5519
#> 
#> $model_results$AIC
#> [1] 223.1037
#> 
#> $model_results$estimate_matrices
#> $model_results$estimate_matrices$lambda
#>          [,1]
#> [1,] 147.7468
#> [2,] 147.7468
#> 
#> $model_results$estimate_matrices$gamma
#>          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
#> [1,] 72.17185 72.17185 72.17185 72.17185 72.17185 72.17185 72.17185 72.17185
#> [2,] 72.17185 72.17185 72.17185 72.17185 72.17185 72.17185 72.17185 72.17185
#>          [,9]    [,10]    [,11]    [,12]
#> [1,] 72.17185 72.17185 72.17185 72.17185
#> [2,] 72.17185 72.17185 72.17185 72.17185
#> 
#> $model_results$estimate_matrices$omega
#>           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#> [1,] 0.1142923 0.1142923 0.1142923 0.1142923 0.1142923 0.1142923 0.1142923
#> [2,] 0.1142923 0.1142923 0.1142923 0.1142923 0.1142923 0.1142923 0.1142923
#>           [,8]      [,9]     [,10]     [,11]     [,12]
#> [1,] 0.1142923 0.1142923 0.1142923 0.1142923 0.1142923
#> [2,] 0.1142923 0.1142923 0.1142923 0.1142923 0.1142923
#> 
#> $model_results$estimate_matrices$pdet
#>           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#> [1,] 0.9538842 0.9538842 0.9538842 0.9538842 0.9538842 0.9538842 0.9538842
#> [2,] 0.9538842 0.9538842 0.9538842 0.9538842 0.9538842 0.9538842 0.9538842
#>           [,8]      [,9]     [,10]     [,11]     [,12]
#> [1,] 0.9538842 0.9538842 0.9538842 0.9538842 0.9538842
#> [2,] 0.9538842 0.9538842 0.9538842 0.9538842 0.9538842
#> 
#> 
#> 
#> $model_data
#> $model_data$K
#> [1] 400
#> 
#> $model_data$nit
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
#> [1,]  134   61   93   84   92   65  105   86   69    84    87    68
#> [2,]  148   93   90  103   94   58   65   65   53    63    67    63
#> 
#> $model_data$l_s_c
#> NULL
#> 
#> $model_data$g_s_c
#> NULL
#> 
#> $model_data$g_t_c
#> NULL
#> 
#> $model_data$o_s_c
#> NULL
#> 
#> $model_data$o_t_c
#> NULL
#> 
#> $model_data$p_s_c
#> NULL
#> 
#> $model_data$p_t_c
#> NULL
#> 
#> $model_data$SMALL_a_CORRECTION
#> [1] FALSE
```

Multi-threading is implemented using R packages doParallel and foreach.
Note that multi-threading is used to split the computation of the
transition probability matrix by rows. This may not be efficient on some
architectures, and will not be efficient for small K: efficiency
increases with increasing K. Alternative choices for multi-core
processing are being considered.
