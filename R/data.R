#' Ancient Murrelet Chick Counts
#'
#' @description Parker et al. (2020) include a data set collected by the 
#' Laskeek Bay Conservation Society on yearly Ancient Murrelet chick counts
#' from the year 1990 to 2006. The data is collected for six sampling sites
#' on East Limestone Island.
#' 
#' @references Parker, M.R.P., Pattison, V. & Cowen, L.L.E. Estimating population abundance using counts from an auxiliary population. Environ Ecol Stat 27, 509–526 (2020). \doi{https://doi.org/10.1007/s10651-020-00455-3}
#'
#' @format A data frame with 6 rows and 17 columns. Each row represents a sampling location, and each column represents a sampling occasion:
#' @source Parker et al. (2020) \doi{https://doi.org/10.1007/s10651-020-00455-3}
"anmu"

#' Golden Eagle Counts Data
#'
#' @description Golden Eagle counts for the years 1993 to 2020, collected by the Rocky Mountain Eagle Research Foundation (RMERF). Counts are made during the spring from April 1st until March 22nd each year. Data is available from Eaglewatch.ca.
#' 
#' @references Rocky Mountain Eagle Research Foundation (RMERF). EagleWatch.ca. (2020).
#'
#' @format A data frame with 28 rows and 11 columns. Each row represents a spring observation period, and each column represents a variable:
#' \describe{
#'   \item{Year}{observation year}
#'   \item{Hours}{hours spent collecting observations}
#'   \item{Eagles}{number of Golden Eagles observed}
#'   \item{TotPrec}{total precipitation measured in mm}
#'   \item{UniqueObservers}{number of principle observers who performed observations}
#'   \item{PeterSherrington_Year}{indicator variable which is 1 for years in which Peter Sherrington was the most prevalent principle observer}
#'   \item{observerPC1}{top 5 principal component scores calculated from the top twelve most prevalent principal observers}
#'   \item{observerPC2}{top 5 principal component scores calculated from the top twelve most prevalent principal observers}
#'   \item{observerPC3}{top 5 principal component scores calculated from the top twelve most prevalent principal observers}
#'   \item{observerPC4}{top 5 principal component scores calculated from the top twelve most prevalent principal observers}
#'   \item{observerPC5}{top 5 principal component scores calculated from the top twelve most prevalent principal observers}
#' }
#' @source RMERF (2020) EagleWatch.ca
"eagles"