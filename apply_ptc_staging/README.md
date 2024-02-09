# Applying PTC staging

The file `ptc_staging.R` contains some functions to help with applying the staging model defined in the manuscript to new data.  The major functions are as follows:

- `ptc.uptake()`: Calculate the tau uptake in each Pattern of Tau Covariance (PTC).
- `ptc.wscores()`: Convert tau uptakes in each PTC into W-scores, using the models described in the paper that are estimated from ADNI-CU or OASIS-CU.
- `ptc.staging()`: Convert binary assessments of tau pathology in each PTC to staging assignments.
- `ptc.staging.pipeline():` Strings the previous three functions into one larger function.

## Requirements

Software-wise,  `ptc_staging.R` requires the R packages [dplyr](https://dplyr.tidyverse.org/) and [this.path](https://cran.r-project.org/web/packages/this.path/index.html).

You must also provide the input tau data, in the form of a dataframe with regional uptakes in 68 cortical gray matter regions.  To see the order/naming of these regions, you can call `tau.roi.names()` or view the "Region" column in `ptcs/ptcs.csv`.

## Examples

The following will show some examples of using the functions in `ptc_staging.R`.  The functions can be loaded with R:

```

```
