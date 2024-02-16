# Applying PTC staging

The file `ptc_staging.R` contains some functions to help with applying the staging model defined in the manuscript to new data.  The major functions are as follows:

- `ptc.uptake()`: Calculate the tau uptake in each Pattern of Tau Covariance (PTC).
- `ptc.wscores()`: Convert tau uptakes in each PTC into W-scores, using the models described in the paper that are estimated from ADNI-CU or OASIS-CU.
- `ptc.staging()`: Convert binary assessments of tau pathology in each PTC to staging assignments.
- `ptc.staging.pipeline():` Strings the previous three functions into one larger function.

## Requirements

Software-wise,  `ptc_staging.R` requires the R packages [dplyr](https://dplyr.tidyverse.org/) and [this.path](https://cran.r-project.org/web/packages/this.path/index.html).

You must also provide the input tau data, in the form of a dataframe with regional uptakes in 68 cortical gray matter regions.  To see the order/naming of these regions, you can call `tau.roi.names()` or view the "Region" column in `ptcs/ptcs.csv`.  Age and sex information is also required for using W-score models.

## Examples

The following will show some examples of using the functions in `ptc_staging.R`.  The functions can be loaded with R:

```R
source(ptc_staging.R)
```

These examples will use the ADNI dataset which is used throughout this repository.  It can be generated using `run_build_datasets.sh`, and will be saved in `derivatives/adni`.

```R
path <- '../derivatives/adni/main_data.csv'
df <- read.csv(path)
df <- df[df$Group == 'TrainingBaseline', ]
```

### Calculating uptakes in PTCs

`ptc.uptake()` can be used to calculate uptakes in each PTC.  This is done by simply multiplying the PTC matrix (normalized W from the NMF outputs) with a matrix containing the regional tau uptake in each FreeSurfer gray matter region.

This step may require some data wrangling to convert your FreeSurfer regional uptakes into the same order as the order used in the paper.

```R
# select the tau ROI columns
# the order expected by ptc.uptake() is the same order used in ADNI
all.columns <- colnames(df)
rois <- all.columns[grepl('CTX_.*_SUVR', all.columns)]
tau.uptakes <- ptc.uptake(df, tau.roi.columns = rois)
```

The last line gives you a dataframe with the tau uptake (SUVR) in each PTC:

```R
> colnames(tau.uptakes)
[1] "MedialTemporal"        "RightParietalTemporal" "LeftParietalTemporal" 
[4] "Precuneus"             "Occipital"             "LateralFrontal"     
[7] "Sensorimotor"          "Orbitofrontal"

> dim(tau.uptakes)
[1] 418   8
```

This function could also be used for calculating the intensity of other biomarkers in each PTC, such as amyloid uptake or gray matter volume!  The units of the output are dependent on the units of the input.

### Converting tau uptakes to W-scores

`ptc.wscore` allows you to use the W-score models that are applied in the manuscript.  This enables a normative approach for converting continuous tau uptakes into binary assessments of tau positivity.  Use the `model` argument to select the ADNI-CU or OASIS-CU model (or use "both" to group the two of them together).

Your data MUST include columns for age and sex, as these are the covariates used in the W-score models trained in the main text.  The sex column should be a binary variable, with 0 to indicate female and 1 to indicate male (the example below does this conversion for the ADNI data).

The `cutoff` argument allows you to optionally binarize the W-scores at a given threshold (2.5 was used in the manuscript):

```R
# adding age & sex to the uptakes table
tau.uptakes$Age <- df$Age
tau.uptakes$Sex <- ifelse(df$Gender == 'Male', 1, 0)

tau.wscores <- ptc.wscores(tau.uptakes,
                          model = 'adni',
                          cutoff = 2.5,
                          age.column = 'Age',
                          sex.column = 'Sex')
```

The above call generates binary assessments of pathology in each PTC.  This is the same data shown in Figure 2A:

```R
> dim(tau.wscores)
[1] 418   8

> colSums(tau.wscores)
       MedialTemporal  LeftParietalTemporal RightParietalTemporal 
                  191                   158                   160 
            Precuneus             Occipital        LateralFrontal 
                  152                   118                   112 
         Sensorimotor         Orbitofrontal 
                   88                    87 
```

### Staging based on binary uptakes in each PTC

The function `ptc.staging` can use binary tau measures in each PTC to generate a staging assignment for each subject: 0, 1, 2, 3, 4, or NS.  This is based on the system shown in Figure 3A & B.

```R
> tau.stages <- ptc.staging(tau.wscores)

> length(tau.stages)
[1] 418

> table(tau.stages)
tau.stages
  0   1   2   3   4  NS 
190  24  29  35  98  42 
```

If you don't want to use the W-score models for binarization (or if you have a different preferred method for converting tau uptakes to binary measures of pathology), you can supply your own binary table:

```R
# simple approach based on an SUVR cutoff
> tau.uptakes <- ptc.uptake(df, tau.roi.columns = rois)
> my.binary <- as.data.frame(ifelse(tau.uptakes >= 1.20, 1, 0))
> my.stages <- ptc.staging(my.binary)
> table(my.stages)
my.stages
  0   1   2   3   4  NS 
123  25  60  37 108  65 
```

### Pipeline function

All the above steps have been stitched into one function, which outputs a list containing the uptakes, W-scores, binary tau measures, and stages:

```R
> output <- ptc.staging.pipeline(df, tau.roi.columns = rois)
```
