
**Group Iterative Multiple Model Estimation (GIMME)**
=====================================================

The GIMME algorithm is a continually maintained R package, gimme.

For up-to-date tutorials for both new and old options, [check out our online *gimme* tutorials](https://tarheels.live/gimme/tutorials/) and [CRAN documentation](https://cran.r-project.org/package=gimme) . 

Developers or those who fixed a bug are invited to submit changes here at the GitHub repository.

**The Basics**
==============

-   GIMME can be used to estimate the unified SEM (uSEM; Kim et al., 2007; Gates et al., 2010).

-   Missing data is not a problem.

-   Heterogeneous data is not a problem:

    -   No "group" or "common" structure will be forced unless it truly describes the majority.

    -   Individual-level nuances will surface after a group or common structure is fit (provided one exists).

    -   If desired, subgroups of individuals with similar patterns of effects will be generated to aid the researcher in finding similar patterns among the varied individual models.

-   Works well with as little as 3 or as many as 20 variables.

-   Can be freely downloaded by installing the package "gimme" in R.

-   Requires at least T = 30 time points per person / per variable, with T = 60 and above recommended. 

**Running GIMME**
=================

**1. Create two new folders (i.e., directories)**

-   Create a source folder for your time series. This can be anywhere that you have permission to read and write.

-   Nothing can be in the source folder other than the time series data.

-   Create an output folder for your results. This must be different from the above folder.

**2. Extract the time series for your variables**

-   Have each variable be a column, with the rows being the observation (e.g., scan in fMRI or a day in daily diary studies).

-   Substitute NA for missing values.

-   Have a separate file for each individual/session.

-   Put each file in the source folder you created in step 1. Do not put anything else in this folder.

-   Files must be either comma, space, or tab delimited.

**3. Installing gimme with R**

-   Open an R script and enter into the console: `install.packages("gimme")`

-   Once gimme has been installed, you will need to load the package by entering: `library(gimme)`

**4. Running gimme**

The *gimme* (or equivalently, *gimmeSEM*) function requires that you input:

-   The path to the directory containing your data

-   How data are separated (e.g., comma-separated values)

-   Whether the data files contain a header row

All other fields are optional and will go to defaults if no user input is provided. If no output directory is indicated, all information is stored as R objects (see tutorial linked above for details).

``` r
fit <- gimme(         # can use "gimme" or "gimmeSEM"
  data = '',          # source directory where your data are 
  out = '',           # output directory where you'd like your output to go
  sep = "",           # how data are separated. "" for space; "," for comma, "\t" for tab-delimited
  header = ,          # TRUE or FALSE, is there a header
  ar = TRUE,          # TRUE (default) or FALSE, start with autoregressive paths open
  plot = TRUE,        # TRUE (default) or FALSE, generate plots
  subgroup = FALSE,   # TRUE or FALSE (default), cluster individuals based on similarities in effects
  paths = NULL,       # option to list paths that will be group-level (semi-confirmatory)
  groupcutoff = .75,  # the proportion that is considered the majority at the group level
  subcutoff = .5      # the proportion that is considered the majority at the subgroup level
)        
```

While *gimme* is running you will see information iterate in the command window. The algorithm will tell you when it is finished.

**Output**
==========

-   The output directory will contain:

    -   **indivPathEstimates**: Contains estimate, standard error, p-value, and z-value for each path and each individual

    -   **summaryFit**: Contains model fit information for individual-level models. If subgroups are requested, this file also indicates the subgroup membership for each individual.

    -   **summaryPathCountMatrix**: Contains counts of total number of paths, both contemporaneous and lagged, estimated for the sample. The row variable is the outcome and the column variable is the predictor variable.

    -   **summaryPathCounts**: Contains summary count information for paths identified at the group-, subgroup (if subgroup = TRUE), and individual-level.

    -   **summaryPathPlots**: Produced if plot = TRUE. Contains figure with group, subgroup (if subgroup = TRUE), and individual-level paths for the sample. Black paths are group-level, green paths are subgroup-level, and grey paths are individual-level, where the thickness of the line represents the count.

-   The subgroup output directory (if subgroup = TRUE) will contain:

    -   **subgroup*k*PathCounts**: Contains counts of relations among lagged and contemporaneous variables for the **k**th subgroup

    -   **subgroup*k*Plot**: Contains plot of group, subgroup, and individual level paths for the **k**th subgroup. Black represents group-level paths, grey represents individual-level paths, and green represents subgroup-level paths.

    -   *Note: if a subgroup of size n = 1 is discovered, subgroup-level output is not produced. Subgroups of size one can be considered outlier cases*

-   In individual output directory (*where id represents the original file name for each individual*):

    -   ***id*Betas**: Contains individual-level estimates of each path for each individual.

    -   ***id*StdErrors**: Contains individual-level standard errors for each path for each individual.

    -   ***id*Plot**: Contains individual-level plots. Red paths represent positive weights and blue paths represent negative weights.

**FAQ**
=======

**How many time points do I need?** This is a difficult question since it will be related to the number of variables you are using. Rules of thumb for any analysis can generally be used: the more the better! Having at least 100 time points is recommended, but adequate results have been obtained in simulation studies with only T = 60.

**Do all individuals have to have the same number of observations (T)?** No.

**How many people do I need in my sample?** For regular *gimme*, reliable results are obtained with as few as 10 participants. Remember that in this context, power to detect effects is determined by the number of time points rather than the number of individuals. Still, having at least 10 individuals helps *gimme* to detect signal from noise by looking for effects that consistently occur.

**What do I do if I obtain an error?** Do some initial trouble-shooting. 1. Ensure that all of your individuals have the same number of variables (columns) in their data sets. 2. Ensure that all variables have variability (i.e., are not constant). *gimme* will let you know if this is the case. 3. Ensure your path directories are correct. 4. Ensure that the columns are variables and the rows contain the observations across time. 5. If all of this is correct, please email the error you received, code used to run *gimme*, and the data (we promise not to use it or share it) to: <gimme@unc.edu>.
