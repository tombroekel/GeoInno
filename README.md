
# regcomplex

The main function in this package is *structural_diversity()* calculates
the measure of **structural diversity** as defined by Broekel (2019)
from patent data. It also features a number of support functions and the
function *complexity_frontier()* used to transform an index of
technological complexity to the so-called **complexity frontier**, which
Mewes & Broekel (2022) use to aggregate the technology-specific index to
the level of a spatial unit (regions or countries).

The package requires the packages igraph, tidyverse, widyr, netdist to
be installed.

## References

Broekel, T., 2019. Using structural diversity to measure the complexity
of technologies. PLoS One 14, 1–27.
<https://doi.org/10.1371/journal.pone.0216856>

Mewes, L., Broekel, T., 2022. Technological complexity and economic
growth of regions. Res. Policy 51, 104156.
<https://doi.org/10.1016/j.respol.2020.104156>

## Installation

You can install the development version of regcomplex like this:

``` r
install.packages("devtools")
library(devtools) 
devtools::install_github("tombroekel/regcomplex", force = T)
```

## Example

This is a basic example which shows you how to use the main function
*structural_diversity()*

``` r
library(regcomplex)

# Create a sample data.frame mimicking common patent data, i.e., a data.frame with four columns: appln_id, technology code, cpc code, year
set.seed(123)
appln_id <- sample(200, size=1000, replace = T)
cpc.1 <- sample(LETTERS[1:10], size=1000, replace = T)
cpc.2 <- sample(LETTERS[1:25], size=1000, replace = T)
cpc <- paste0(cpc.1,cpc.2)

df <- data.frame(appln_id=appln_id, cpc=cpc) %>% arrange(appln_id)
df <- df %>% mutate(tech = substr(cpc,1,1))
year <- sample(c(1999:2005), size=length(unique(appln_id)), replace = T)
df <- df %>% mutate(year = year[appln_id])

#### Calculate structural diversity #### 

complexity <- structural_diversity(df)
```
