#' K Calculation (Vectorized)
#'
#' Compute the K values of the LP3 distribution used by Bulletin 17B
#'
#' @param AEPs numeric - vector of annual exceedance probabilities
#' @param skew numeric - skewness of distribution
#'
#' @return vector of K values
#'
#' @keywords internal
#'
k_17b <- function(AEPs, skew){
  
  K <- sapply(AEPs, kcalc, skew)
  
  return(K)
  
}


#' K Calculation
#'
#' Compute the K values of the LP3 distribution used by Bulletin 17B
#'
#' @param AEP numeric - annual exceedance probability
#' @param skew numeric - skewness of distribution
#'
#' @return vector of K values
#'
#' @keywords internal
#'
kcalc <- function(AEP, skew){
  
  if(AEP < 0 | AEP > 1){
    stop(paste("Invalid probability:", prob))
  }
  
  prob = 1 - AEP #Convert to nonexceedance probability
  
  K <- .Fortran("s_kfxx",
                    as.double(skew), #Number of observations, censored or uncensored
                    as.double(prob), #Vector of lower interval values
                    as.double(-99) #Computed K value
  )
  
  return(K[[3]])
  
}