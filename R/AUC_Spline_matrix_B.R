#' @title Spline interpolation Method - Matrix of the zero orer derivative coefficients
#' @description In the area under the curve calculation using the spline interpolation method, the vector of the second derivative of the outcome of interest \emph{Y} is expressed as \eqn{A Y^{''} = B Y  + F}. This function calculate calculate the matrix B.
#' @param time a numerical vector of time points (x-axis cooordinates).
#' @return  a tridiagonal matrix corresponding to the weights of the variable of interest in the spline interpolation method. In this version, the matrix is build considering the "not-a-knot" spline boundary conditions.
#' @details METTRE l'équation au format Latex reliant Y'' et Y + Mettre la definition de la matrice.
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
