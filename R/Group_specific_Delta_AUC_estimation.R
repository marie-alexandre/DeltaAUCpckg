#' @title Difference of AUC of Two Group-Specific Polynomial Marginal Dynamics
#' @description This function estimates the difference of area under the curve of marginal dynamics from two groups when marginal dynamics are modeled by group-structured polynomials or B-spline curves.

#' @param MEM_Pol_group A list with similar structure than the output provided by the function \link[DeltaAUCpckg]{MEM_Polynomial_Group_structure}. 
#' 
#' A list containing: \tabular{ll}{
#' \code{Model_estimation} \tab either the vector of the marginal parameter estimates for ALL the groups involved in the group-structured polynomial MEM, or a list containing at least this vector labeled 'beta'. \cr
#' \code{Model_features} \tab a list of at least 2 elements: \cr
#' \tab 1. \code{Groups}  -  a vector indicating the names of ALL the groups involved in the group-structured model. \cr
#' \tab 2. \code{Marginal.dyn.feature}  -  a list summarizing the features of the marginal dynamics defined in the model:  \cr
#' \tab \itemize{
#' \item \code{dynamic.type} - a character scalar indicating the chosen type of marginal dynamics. Options are 'polynomial' or 'spline'
#' \item \code{intercept} -  a logical vector summarizing choices about global and group-specific intercepts (Number of groups + 1) elements whose elements are named as ('global.intercept','group.intercept1', ..., 'group.interceptG') if G Groups are involved in the model. For each element of the vector, if TRUE, the considered intercept is considered as included in the model. 
#'  
#' If \code{dynamic.type} is defined as 'polynomial': 
#' \item \code{polynomial.degree} - an integer vector indicating the degree of polynomial functions, one value for each group. 
#' 
#' If \code{dynamic.type} is defined as 'spline':
#' \item \code{spline.degree} - an integer vector indicating the degree of B-spline curves, one for each group. 
#' \item \code{knots} - a list of group-specific internal knots used to build B-spline basis (one numerical vector for each group) (see \link[splines]{bs} for more details).
#' \item \code{df} - a numerical vector of group-specific degrees of freedom used to build B-spline basis, (one for each group). 
#' \item \code{boundary.knots} - a list of group-specific boundary knots used to build B-spline  basis (one vector for each group) (see \link[splines]{bs} for more details).
#' } \cr
#' }
#' 
#' @param Group1 a character scalar indicating the name of the first group whose marginal dynamics must be considered. This group name must belong to the set of group involved in the MEM for which we want to estimate the AUC (see \code{Groups} vector in \code{MEM_Pol_group}). 
#' @param Group2 a character scalar indicating the name of the second group whose marginal dynamics must be considered. This group name must belong to the set of group involved in the MEM for which we want to estimate the AUC (see \code{Groups} vector in \code{MEM_Pol_group}).
#' @param time.G1 a numerical vector of time points (x-axis coordinates) to use for the Group1 AUC calculation.
#' @param time.G2 a numerical vector of time points (x-axis coordinates) to use for the Group2 AUC calculation.
#' @param method a character scalar indicating the interpolation method to use to estimate the AUC. Options are 'trapezoid' (default), 'lagrange' and 'spline'. In this verion, the 'spline' interpolation is implemented with "not-a-knot" spline boundary conditions.
#' @param Averaged a logical scalar. If TRUE, the function return the difference of normalized AUC (nAUC) where nAUC is computated as the AUC divided by the range of time of calculation. If FALSE (dafault), the classic AUC is calculated.
#' @return 
#' A numerical scalar defined as \eqn{\Delta}AUC = AUC2 - AUC1 (or \eqn{\Delta}nAUC = nAUC2 - nAUC1)  with AUC1 (or nAUC1) and  AUC2 (or nAUC) being respectively estimated as the AUC (or nAUC) for the Group1 and for the Group2.
##' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname Group_specific_Delta_AUC_estimation
#' @export 
Group_specific_Delta_AUC_estimation <- function(MEM_Pol_group,Group1,Group2,time.G1,time.G2,method="trapezoid",Averaged=FALSE){
  AUC_Group1 <- Group_specific_AUC_estimation(MEM_Pol_group=MEM_Pol_group,time=time.G1,Groups=Group1,method=method,Averaged=Averaged) 
  AUC_Group2 <- Group_specific_AUC_estimation(MEM_Pol_group=MEM_Pol_group,time=time.G2,Groups=Group2,method=method,Averaged=Averaged) 
  Delta_AUC <- as.numeric(AUC_Group2 - AUC_Group1)
  return(Delta_AUC)
}
