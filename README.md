# SpatIndex <img src="./inst/logos/SpatIndex.png" align="right" width="150" />

Spatial indices to monitor droughts and global change 

This R package provides a simple way of computing spatially-distributed drought and extreme indices. The users are able to compute meteorological, agricultural, and hydrological indices as well as the extreme CLIMDEX just by providing raster files of the required hydroclimatological variables. The drought and extreme indices that are implemented in **SpatIndex** are sumarised in the following figure:

<img src="./inst/logos/Fig_indices.png" align="center" width="1000" />

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
