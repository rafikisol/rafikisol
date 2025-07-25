#'WoSIS data
#'
#'A small data frame of a subset WoSIS data in Southern Malawi. The format is how
#'it should be if problems are experienced when using the catSpline function.
#'
#'@format ## `Neno`
#' A data frame with 212 rows and 4 columns:
#' \describe{
#'   \item{id}{Profile number or order, should align with a shapefile}
#'   \item{top}{top depth of a horizon in the profile}
#'   \item{botoom}{bottom depth of a horizon in the profile}
#'   \item{Class}{USDA soil texture classes}
#' }
#'@source <https://www.isric.org/explore/isric-soil-data-hub>

'Neno'

#'WoSIS harmonised soil texture classes
#'
#'Global data of harmonised USDA soil texture classes of 101,610 profiles with their
#'uncertainties in a Geopackage. Hormonised to GlobalSoilMap standardised depths.
#'
#'@format ## `global_data`
#' List with two sf layers:
#' \describe{
#'   \item{Harmonised}{The harmonised soil texture classes}
#'   \item{Uncertainties}{The uncertainties associated with each profiles and depth}
#' }
#'@source <https://www.isric.org/explore/isric-soil-data-hub>

'global_data'
