
# AttBART

<!-- badges: start -->
<!-- badges: end -->

NOTE: Currently produces constant predictions. Not debugged.


The goal of AttBART is to implement Attention Weighted Bayesian Additive Regression Trees by combining Attention Weighted Boosted Trees (Konstantinov et al. 2022) with Bayesian Additive Regression Trees (Chipman et al. 2010).

[Chipman, H. A., George, E. I., & McCulloch, R. E. (2010). _BART: Bayesian additive regression trees_. The Annals of Applied Statistics, 4(1), 266-298.](https://doi.org/10.1214/09-AOAS285)

[Konstantinov, A., Utkin, L., & Kirpichenko, S. (2022, April). _AGBoost: Attention-based Modification of Gradient Boosting Machine_. In 2022 31st Conference of Open Innovations Association (FRUCT) (pp. 96-101). IEEE.](https://arxiv.org/abs/2207.05724)

## Installation

You can install the development version of AttBART like so:

``` r
library(devtools)
install_github("EoghanONeill/AttBART")
```

## Example

``` r
library(AttBART)
## basic example code


# Simulate a Friedman data set
friedman_data = function(n, num_cov, sd_error){
  x = matrix(runif(n*num_cov),n,num_cov)
  y = 10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5] + rnorm(n, sd=sd_error)
  return(list(y = y,
              x = x))
}
# Training data
data = friedman_data(200, 10, 1)
y = data$y
x = data$x

# Test data
data_test = friedman_data(100, 10, 1)
y.test = data_test$y
x.test = data_test$x


fit.attbart <- attBart_no_w(x,y, feature_weighting = TRUE)

y.test.hat_abart <- predict_attbart_test(fit.attbart, x.test, type = "mean")

plot(y.test, y.test.hat_abart); abline(0, 1)
cor(y.test, y.test.hat_abart)

```

