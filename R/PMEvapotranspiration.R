
#' @title PMEvapotranspiration
#' @description FAO-56 Penman-Monteith evapotranspiration
#'
#' @param T_max Maximum temperature in degree Celsius
#' @param T_min Minimum temperature in degree Celsius
#' @param Latitude Latitude in decimal degree
#' @param altitude altitude in meter
#' @param RH_min Minimum Relative humidity in percentage
#' @param RH_max Maximum Relative humidity in percentage
#' @param u2 The average wind speed (m/s) measured at 2 m above the ground level
#' @param n Actual daily sunshine duration
#' @param J number of the day in the year between 1 (1 January) and 365 or 366 (31 December)
#'
#' @return Returns a list containing the calculated reference evapotranspiration (ET0) and other intermediate
#' values used in the calculation. The list includes:
#' \itemize{
#'   \item \code{ET0}: Reference evapotranspiration calculated using FAO-56 Penman-Monteith Equation.
#'   \item \code{T_mean}: Mean temperature.
#'   \item \code{P}: Atmospheric pressure.
#'   \item \code{gamma}: Psychrometric constant.
#'   \item \code{delta}: Slope of vapor pressure curve.
#'   \item \code{e_s}: Saturation vapor pressure.
#'   \item \code{e_a}: Actual vapor pressure.
#'   \item \code{Ra}: Extraterrestrial radiation.
#'   \item \code{Rso}: Clear sky solar radiation.
#'   \item \code{Rs}: Solar radiation.
#'   \item \code{Rns}: Net shortwave radiation.
#'   \item \code{Rnl}: Net longwave radiation.
#'   \item \code{Rn}: Net radiation.
#' }
#' @examples{
#' calculateET0(T_max = 25, T_min = 15, Latitude = 45,
#'  altitude = 100, RH_min = 20,
#'  RH_max = 80, u2 = 2, n = 10, J = 150)
#'}
#' @export
#'
#' @references
#' #' \itemize{
#'\item Córdova, M., Carrillo-Rojas, G., Crespo, P., Wilcox, B., and Célleri, R. (2015). Evaluation of the Penman-Monteith (FAO 56 PM) method for calculating reference evapotranspiration using limited data. Mountain Research and Development, 35(3), 230-239.
#' \item Debnath, S., Adamala, S., and Raghuwanshi, N. S. (2015). Sensitivity analysis of FAO-56 Penman-Monteith method for different agro-ecological regions of India. Environmental Processes, 2, 689-704.
#' }

calculateET0 <- function(T_max, T_min, Latitude, altitude, RH_min, RH_max, u2, n, J) {
  sigma <- 4.903*10^-9 #  the Stefan-Boltzmann constant MJ⋅K^−4⋅^m−2⋅da^y−1
  Gsc <- 0.0820  # Solar constant [MJ m^-2 min^-1]
  T_mean <- (T_min + T_max) / 2
  Z <- altitude
  P <- (101.3*((293 - 0.0065*Z)/293)^5.26)
  gamma <- 0.000665 * P
  delta <- 4098*(0.6108*exp((17.27 * T_mean) / (T_mean + 237.3))) / ((T_mean + 237.3)^2)
  e_Tmin <- 0.6108 * exp((17.27 * T_min) / (T_min + 237.3))
  e_Tmax <- 0.6108 * exp((17.27 * T_max) / (T_max + 237.3))
  e_s <- (e_Tmin + e_Tmax) / 2
  e_a <- (e_Tmin * (RH_max / 100) + e_Tmax * (RH_min / 100)) / 2
  phi <- Latitude * pi / 180
  dr <- 1 + 0.033 * cos(2 * pi * J / 365)
  declination <- 0.409 * sin((2 * pi * J / 365) - 1.39)
  omega_s <- acos(-tan(phi) * tan(declination))
  Ra <- (24 * 60 / pi) * Gsc * dr * (omega_s * sin(phi) * sin(declination) + cos(phi) * cos(declination) * sin(omega_s))
  Rso <- (0.75 + 2 * 10^-5 * Z) * Ra
  N <- (24 / pi) * omega_s
  a <- 0.25
  b <- 0.50
  Rs <- (a + b * (n / N)) * Ra
  Rns <- (1 - 0.23) * Rs
  Rnl <- sigma * (((T_max + 273.16)^4 + (T_min + 273.16)^4) / 2) * (0.34 - 0.14 * sqrt(e_a)) * (1.35 * (Rs / Rso) - 0.35)
  Rn <- Rns - Rnl
  G <- 0 #soil heat flux is often assumed to be negligible for daily time steps
  ET0 <- (0.408 * delta * (Rn - G) + gamma * (900 / (T_mean + 273)) * u2 * (e_s - e_a)) / (delta + gamma * (1 + 0.34 * u2))

  # Return a list containing all the calculated values
  Result<- list(ET0 = ET0, T_mean = T_mean, P = P, gamma = gamma, delta = delta, e_s = e_s, e_a = e_a, Ra = Ra, Rso = Rso, Rs = Rs, Rns = Rns, Rnl = Rnl, Rn = Rn)
return(Result)
}


