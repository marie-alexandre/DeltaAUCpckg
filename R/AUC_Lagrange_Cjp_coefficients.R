#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param ind_j PARAM_DESCRIPTION
#' @param ind_p PARAM_DESCRIPTION
#' @param t PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' @rdname AUC_Lagrange_Cjp_coefficients
#' @export 
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
