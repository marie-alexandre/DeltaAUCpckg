#' @title Optimal Research B-Spline Internal Knots
#' @description This function uses AIC criterion to estimate optimal position and number of knots by fitting free-knot splines to data with one independent and one dependent variable. It is assumed that knots are estimated for least-squares splines with no penalty using genetic algorithm. 
#' @param data a data frame containing the independent variable (\emph{y}) and the dependent one (\emph{x}). 
#' @param degree an integer scalar indicating the degree of the spline it. By default, this value is fixed at 3
#' @param minknot an integer scalar indicating the minimum number of knots to consider. By default, this variable is fixed at 2
#' @param maxknot an integer scalar indicating the maximum number of knots to consider. By default, this variable is fixed at 2.
#' @param ... Further arguments to be passed (see \link[freeknotsplines]{freeknotfit} for more details).
#' 
#' @return A numerical vector of optimal knots whose number can varied from \code{minknot} to \code{maxknot}
#' 
#' @seealso 
#'  \code{\link[freeknotsplines]{freeknotfit}} for more details about the method. 
#'  
#' @rdname Optimal_knot_research
#' @export 
#' @importFrom freeknotsplines freelsgen
#' 
Optimal_knot_research <- function(data,degree=3,minknot=2,maxknot=2,...){
  knot_research <- NULL
  
  tested_numknots <- seq(minknot,maxknot)
  
  tmp_research <- lapply(seq(1,length(tested_numknots)),function(i){
    res <- tryCatch(freeknotsplines::freelsgen(x=data$x,y=data$y,degree=degree,numknot=tested_numknots[i],...),
                    error=function(cond){"error"})
    return(res)
  })
  if(length(unique(format(tmp_research))) == 1 & class(tmp_research[[1]]) != "freekt"){
    stop("Unable to estimate optimal knots with these arguments")
  }else{
    tmp_research <- tmp_research[which(unlist(lapply(tmp_research, function(x) class(x))) == "freekt")]
    AIC <- unlist(lapply(tmp_research, function(x){tryCatch(AIC(x), error=function(cond){return(Inf)})}))
    if(abs(min(AIC)) < Inf){
      knot_research <- tmp_research[[which.min(AIC)]]
    }else{
      stop("Unable to estimate optimate knots of finite value of AIC.")
    }
  }
  
  if(is.null(knot_research)){
    return(NULL)
  }else{
    return(knot_research@optknot)
  }
}
