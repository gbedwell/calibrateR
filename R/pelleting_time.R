pelleting_time <- function( speed, 
                            speed.units = c( "rpm", "rcf" ),
                            rmin,
                            rmax,
                            S ){
  
  if( rmax < 50 ){
    warning( "Ensure that radius values are given in mm units.",
             call. = FALSE )
  }
  
  speed.units = match.arg( speed.units )
  
  if( speed.units == "rcf" ){
    speed = rpm_rcf_convert( from = "rcf", 
                             to = "rpm", 
                             rmax = rmax, 
                             speed = speed )
  }
  
  k <- 2.53E11 * log( rmax / rmin ) / speed^2
  hours <- k / S
  mins <- hours * 60
  
  return( c( hours = hours, minutes = mins ) )
}


