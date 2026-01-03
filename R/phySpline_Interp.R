#'@title Interpolate phySplines
#'
#'@name phySpline_Interp
#'
#'@author Trevan Flynn
#'@description
#'A helper function to interpolate values from the phySpline() function for visualisation purposes from
#'the coefficients and optional classes. Can interpolate categorical or continuous phySplines.
#'
#'@param depths a vector of the depths to interpolate on to.
#'@param coeffs the coefficients for the profile returned from phySpline().
#'@param looup for categorical splines, the lookup table used from phySplines().
#'@param vlow the lowest value to interpolate
#'@param vhigh the highest value to interpolate
#'
#'@return a dataframe with the depths, the value per depth and optionally the class.
#'
#'@export

phySpline_Interp <- function(depths = seq(0, 160, by = 1), coeffs, lookup = NULL, vlow = NULL, vhigh = NULL) {

  latent_vals <- sapply(depths, function(z) {
    # Find which horizon depth z belongs to
    idx <- which(z >= coeffs$top & z < coeffs$bottom)

    if(length(idx) == 0){
      # Extrapolation
      if(z < min(coeffs$top)) {
        # Above first horizon: linear extrapolate using first horizon
        u_j <- coeffs$top[1]
        alpha_j <- coeffs$alpha[1]
        b_j <- coeffs$b[1]
        gamma_j <- coeffs$gamma[1]
        return(alpha_j + b_j * (z - u_j) + gamma_j * (z - u_j)^2)
      } else if(z > max(coeffs$bottom)) {
        # Below last horizon: linear extrapolate using last horizon
        u_j <- coeffs$top[nrow(coeffs)]
        alpha_j <- coeffs$alpha[nrow(coeffs)]
        b_j <- coeffs$b[nrow(coeffs)]
        gamma_j <- coeffs$gamma[nrow(coeffs)]
        return(alpha_j + b_j * (z - u_j) + gamma_j * (z - u_j)^2)
      } else {
        return(NA)
      }
    }

    # Inside horizon
    u_j     <- coeffs$top[idx]
    alpha_j <- coeffs$alpha[idx]
    b_j     <- coeffs$b[idx]
    gamma_j <- coeffs$gamma[idx]

    alpha_j + b_j * (z - u_j) + gamma_j * (z - u_j)^2
  })

  # Apply vlow / vhigh bounds if provided
  if(!is.null(vlow)) latent_vals <- pmax(latent_vals, vlow)
  if(!is.null(vhigh)) latent_vals <- pmin(latent_vals, vhigh)

  # Map to classes if lookup is provided
  if(!is.null(lookup)) {
    class_vals <- sapply(latent_vals, function(v) {
      if(is.na(v)) return(NA)
      lookup$class[which.min(abs(v - lookup$numeric))]
    })
  } else {
    class_vals <- rep(NA, length(latent_vals))
  }

  data.frame(
    Depths = depths,
    Value = latent_vals,
    Class  = class_vals
  )
}
