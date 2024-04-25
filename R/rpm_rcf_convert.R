rpm_rcf_convert <- function( from = c( "rpm", "rcf" ),
                             to = c( "rpm", "rcf" ),
                             rmax,
                             speed ){
  if( rmax < 50 ){
    warning( "Ensure that rmax is given in mm units.",
             call. = FALSE )
  }

  if( from == "rpm" && to == "rcf" ){
    new.speed <- 1.118E-5 * ( rmax / 10 ) * ( speed )^2
  } else{
    if( from == "rcf" && to == "rpm" ){
      new.speed = sqrt( speed / ( 1.118E-5 * ( rmax / 10 ) ) )
    } else{
      stop( "Can only go from RPM to RCF or vice versa.",
            call. = FALSE )
    }
  }
  return( new.speed )
}
