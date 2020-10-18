#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param MEM_Pol_group PARAM_DESCRIPTION
#' @param Group1 PARAM_DESCRIPTION
#' @param Group2 PARAM_DESCRIPTION
#' @param time.G1 PARAM_DESCRIPTION
#' @param time.G2 PARAM_DESCRIPTION
#' @param method PARAM_DESCRIPTION, Default: 'trapezoid'
#' @param Averaged PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
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
