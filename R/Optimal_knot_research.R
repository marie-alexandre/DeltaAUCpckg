#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data PARAM_DESCRIPTION
#' @param degree PARAM_DESCRIPTION
#' @param minknot PARAM_DESCRIPTION
#' @param maxknot PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @seealso 
#'  \code{\link[freeknotsplines]{freeknotfit}}
#' @rdname Optimal_knot_research
#' @export 
#' @importFrom freeknotsplines freelsgen
Optimal_knot_research <- function(data,degree,minknot,maxknot){
  knot_research <- NULL
  
  tested_numknots <- seq(minknot,maxknot)
  
  tmp_research <- lapply(seq(1,length(tested_numknots)),function(i){
    res <- tryCatch(freeknotsplines::freelsgen(x=data$x,y=data$y,degree=degree,numknot=tested_numknots[i],seed=555,stream=0),
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
