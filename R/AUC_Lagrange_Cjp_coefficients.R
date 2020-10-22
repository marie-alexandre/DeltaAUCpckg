#' @title Time Dependent Coefficients Cjp for AUC Lagrange Interpolation Method
#' @description This function calculate the time-dependent coefficients Cjp involved in the calculation of the are under the curve when the Lagrange interpolation method is used.
#' 
#' @param ind_j a numerical scalar indicating the value of the index j.
#' @param ind_p a numerical scalar indicating the value of the index p
#' @param t a numerical vector of time points (x-axis coordinates) to consider for the AUC calculation.
#' @return a numerical scalar corresponding to the coefficient Cjp evaluated for j = \code{ind_j} and p = \code{ind_p}.
#' @details The coefficients Cij are defined as (METTRE l'equation sous format Latex).
#' @rdname AUC_Lagrange_Cjp_coefficients
#' @export 
#' 
AUC_Lagrange_Cjp_coefficients <- function(ind_j,ind_p,t){
  if(ind_j == 2){
    term1 <- (t[ind_j] - t[ind_j-1])*prod(sapply(seq(0,2),function(l){
      if(l!=ind_p){
        return(t[1+l])
      }else{
        return(1)
      }
    },simplify=TRUE))
    term2 <- (t[ind_j]^2 - t[ind_j-1]^2)/2*sum(sapply(seq(0,2),function(l){
      if(l!=ind_p){
        return(t[1+l])
      }else{
        return(0)
      }
    },simplify=TRUE))
    term3 <- (t[ind_j]^3 - t[ind_j-1]^3)/3
    Cjp <- term1 - term2 + term3
  }else if(ind_j == length(t)){
    term1 <- (t[ind_j] - t[ind_j-1])*prod(sapply(seq(0,2),function(l){
      if(l!=ind_p){
        return(t[ind_j-2+l])
      }else{
        return(1)
      }
    },simplify=TRUE))
    term2 <- (t[ind_j]^2 - t[ind_j-1]^2)/2*sum(sapply(seq(0,2),function(l){
      if(l!=ind_p){
        return(t[ind_j-2+l])
      }else{
        return(0)
      }
    },simplify=TRUE))
    term3 <- (t[ind_j]^3 - t[ind_j-1]^3)/3
    Cjp <- term1 - term2 + term3
  }else{
    term1 <- (t[ind_j] - t[ind_j-1])*prod(sapply(seq(0,3),function(l){
      if(l!=ind_p){
        return(t[ind_j-2+l])
      }else{
        return(1)
      }
    },simplify=TRUE))
    term2 <- (t[ind_j]^2 - t[ind_j-1]^2)/2*sum(sapply(seq(0,2),function(l1){
      if(l1!=ind_p){
        res <- 0
        for(l2 in seq(l1+1,3)){
          if(l2!=ind_p){
            res <- res + t[ind_j-2+l1]*t[ind_j-2+l2]
          }
        }
        return(res)
      }else{
        return(0)
      }
    },simplify=TRUE))
    term3 <- (t[ind_j]^3 - t[ind_j-1]^3)/3*sum(sapply(seq(0,3),function(l){
      if(l!=ind_p){
        return(t[ind_j-2+l])
      }else{
        return(0)
      }
    },simplify=TRUE))
    term4 <- (t[ind_j]^4 - t[ind_j-1]^4)/4
    Cjp <- -term1 + term2 - term3 + term4
  }
  return(Cjp)
}
