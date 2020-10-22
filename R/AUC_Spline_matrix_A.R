#' @title Spline Interpolation Method - Matrix of Second Derivative Coefficients
#' @description In the area under the curve calculation using the spline interpolation method, the vector of the second derivative of the outcome of interest \emph{Y} is expressed as \eqn{A Y^{''} = B Y  + F}. This function calculate calculate the matrix A.
#' @param time a numerical vector of time points (x-axis cooordinates).
#' @return a tridiagonal matrix corresponding to the weights of the second derivative of the variable of interest in the spline interpolation method. In this version, the matrix is build considering the "not-a-knot" spline boundary conditions.
#' @details METTRE l'équation au format Latex reliant Y'' et Y + Mettre la définition de la matrice. 
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
