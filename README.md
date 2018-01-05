# Transfer-Entropy
Code to implement transfer entropy (Shannon)

# Usage
``` r
library(RTransferEntropy)
set.seed(42)

dt <- data.table(series = rnorm(10000))
dt[, sample := code_sample(series)]
dt
#>            series sample
#>     1:  1.3709584      2
#>     2: -0.5646982      2
#>     3:  0.3631284      2
#>     4:  0.6328626      2
#>     5:  0.4042683      2
#>    ---                  
#>  9996: -0.5653259      2
#>  9997: -1.4819016      2
#>  9998: -0.2989931      2
#>  9999:  0.2146498      2
#> 10000:  0.7732340      2
```
# Example using simulated data 
Simulate a simple model to obtain two time series that are not independent (see simulation study in Dimpfl and Peter (2013)),
i.e. one time series is lag of the other plus noise. In this case, one expects significant information flow from x to y 
and none from y to x.

``` r
n<- 100000
x <- rep(0, n + 1)
y <- rep(0, n + 1)

for (i in seq(n)) {
  x[i + 1] <- 0.2 * x[i] + rnorm(1, 0, 2)
  y[i + 1] <- x[i] + rnorm(1, 0, 2)
}

x <- x[-1]
y <- y[-1]
