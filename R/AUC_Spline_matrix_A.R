#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param time PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname AUC_Spline_matrix_A
#' @export 
AUC_Spline_matrix_A <- function(time){
  m <- length(time)
  hj <- NULL  # Warning: length(hj)=m-1
  for(j in 2:m){
    hj <- c(hj,(time[j]-time[j-1]))
  }
  Matrix_A <- matrix(0,ncol=m,nrow=m)
  Matrix_A[1,1] <- 1/hj[1] ; Matrix_A[1,2] <- (1/hj[1] + 1/hj[2]) ; Matrix_A[1,3] <- 1/hj[2]
  Matrix_A[m,m-2] <- 1/hj[length(hj)-1] ; Matrix_A[m,m-1] <- (1/hj[length(hj)-1]+1/hj[length(hj)]) ; Matrix_A[m,m] <- 1/hj[length(hj)]
  
  for(j in 2:(m-1)){
    Matrix_A[j,j-1] <- hj[j-1]/6
    Matrix_A[j,j] <- (hj[j-1]+hj[j])/3
    Matrix_A[j,j+1] <- hj[j]/6
  }
  return(Matrix_A)
}
