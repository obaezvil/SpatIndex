# SpatIndex <img src="./inst/logos/SpatIndex.png" align="right" width="150" />

Spatial indices to monitor global change 

Simple way of computing hydroclimatological indices in R. The users are able to compute state-of-the-art indices by providing the file path to daily tiff files.

The indices included in this package to date are:

- **Rx1day**: Maximum 1-day precipitation

- **Rx5day**: Maximum consecutive 5-day precipitation

- **SDII**: Simple precipitation intensity index

- **R10mm**: Annual count of days when PRCP ≥ 10mm

- **R20mm**: Annual count of days when PRCP ≥ 20mm

- **Rnnmm**: Annual count of days when PRCP ≥ nn mm, where nn is a user-defined threshold

- **CDD**: Maximum length of dry spell: maximum number of consecutive days with RR < 1mm

- **CWD**: Maximum length of wet spell: maximum number of consecutive days with RR ≥ 1mm

- **R95p**: Annual total PRCP when RR > 95th percentile

- **R99p**: Annual total PRCP when RR > 99th percentile

- **R95pTOT**: Contribution to total precipitation from very wet days

- **R99pTOT**: Contribution to total precipitation from extremely wet days

- **PRCPTOT**: Annual total precipitation on wet days


### Installation

Install development version from GitHub :

```r
remotes::install_github("obaezvil/SpatInex")
```

The functions have the same name of the indices. For example, if the user want to calculate the maximum 1-day precipitation (Rx1day), the function with the same name should be applied:

```r
Rx1day(rst.path, vct = NULL,  temporal.scale = c("total", "annual"),
  start.date,  end.date,  date.fmt = "%Y-%m-%d",  pattern = NULL)
```

Where:

- *rst.path* File path to daily raster files for the period of analysis. These files should include the date in any format.

- *vct* Vector file of the study area (Optional). It will be used to crop the spatial extent of the raster files if required. This parameter is set to 'NULL'; therefore, the index will be calculated over the entire raster extent by default.

- *temporal.scale* either 'total' to use all the period of record (i.e., all the files in the parent folder 'rst.path'), or 'annual' to compute the index annually (From Jan to Dec).

- *start.date* Position where the dates of the raster files start. For example, for files named 'ProductA_1989-02-28.tif', the *start.date* will indicate the position where the date starts (i.e., 10).

- *end.date* Position where the dates of the raster files end. For example, for files named 'ProductA_1989-02-28.tif', the *end.date* will indicate the position where the date ends (i.e., 19).

- *date.fmt* Format of the dates included in the file names (default = "%Y-%m-%d").

- *pattern* Set to NULL. Is there a specific pattern to list the raster files?

An example of the implementation of **Rx1day** over Sub-Saharan Africa for 1981--2021 can be observed in the following Figure:

<img src="./inst/logos/Sub-Saharan_Africa_Rx1day.png" align="center" width="900" />
