
# GeoInno

The package’s aim is to provide a set of tools for working with patent
data that are frequently employed in the context of research related to
the Geography of Innovation literature. It will be continously expanded.
At this stage, its main function is *structural_diversity()*, which
calculates the measure of **structural diversity** as defined by Broekel
(2019) from patent data. It also features a number of support functions
and the function *complexity_frontier()* used to transform an index of
technological complexity to the so-called **complexity frontier**, which
Mewes & Broekel (2022) use to aggregate the technology-specific index to
the level of a spatial unit (regions or countries).

The package requires the packages **igraph**, **tidyverse**,
**widyr**,**future**, **future.apply**, **data.table**, **netdist** to
be installed. The package **netdist** can be found here:
<https://github.com/alan-turing-institute/network-comparison/>

## References

Broekel, T., 2019. Using structural diversity to measure the complexity
of technologies. PLoS One 14, 1–27.
<https://doi.org/10.1371/journal.pone.0216856>

Mewes, L., Broekel, T., 2022. Technological complexity and economic
growth of regions. Res. Policy 51, 104156.
<https://doi.org/10.1016/j.respol.2020.104156>

## Installation

You can install the development version of GeoInno like this:

``` r
install.packages("devtools")
library(devtools) 
devtools::install_github("tombroekel/GeoInno", force = T)
```

## Example

This is a basic example which shows you how to use the main function
*structural_diversity()*

``` r
library(GeoInno)

#Calculate structural diversity using artificial data set pat.df, which is created by create_sample_data()

complexity <- structural_diversity(pat.df)
```
