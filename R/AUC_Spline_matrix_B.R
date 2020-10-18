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
#' @rdname AUC_Spline_matrix_B
#' @export 
AUC_Spline_matrix_B <- function(time){
  m <- length(time)
  hj <- NULL  # Warning: length(hj)=m-1
  for(j in 2:m){
    hj <- c(hj,(time[j]-time[j-1]))
  }
  Matrix_B <- matrix(0,ncol=m,nrow=m)
  
  for(j in 2:(m-1)){
    Matrix_B[j,j-1] <- 1/hj[j-1] 
    Matrix_B[j,j] <- -(1/hj[j-1]+1/hj[j])
    Matrix_B[j,j+1] <- 1/hj[j]
  }
  return(Matrix_B)
}
