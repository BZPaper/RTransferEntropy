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
#>     1:  1.3709584      1
#>     2: -0.5646982      1
#>     3:  0.3631284      1
#>     4:  0.6328626      1
#>     5:  0.4042683      1
#>    ---                  
#>  9996: -0.5653259      3
#>  9997: -1.4819016      3
#>  9998: -0.2989931      3
#>  9999:  0.2146498      3
#> 10000:  0.7732340      3
```
