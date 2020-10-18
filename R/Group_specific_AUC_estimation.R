#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param MEM_Pol_group PARAM_DESCRIPTION
#' @param time PARAM_DESCRIPTION
#' @param Groups PARAM_DESCRIPTION, Default: NULL
#' @param method PARAM_DESCRIPTION, Default: 'trapezoid'
#' @param Averaged PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' @seealso 
#'  \code{\link[splines]{bs}}
#' @rdname Group_specific_AUC_estimation
#' @export 
#' @importFrom splines bs
Group_specific_AUC_estimation <- function(MEM_Pol_group,time,Groups=NULL,method="trapezoid",Averaged=FALSE){
  '%notin%' <- Negate('%in%') 

  # Step 1: Verification of the type of arguments
  # ----- #
  Check_argument_Group_specific_AUC(MEM_Pol_group,time,Groups,method,Averaged)
  
  Model_features <- MEM_Pol_group$Model_features
  Marginal_dynamics <- Model_features$Marginal.dyn.feature
  if(is.null(Groups)){
    Groups <- Model_features$Groups
  }
  if(is.numeric(time)){
    time <- lapply(seq(1,length(Groups)),function(g) return(time))
  }
  
  # Extraction of population parameters according to their groups
  Population_params <- MEM_Pol_group$Model_estimation$beta
  MEM_groups <- Model_features$Groups
  
  global_intercept <- Marginal_dynamics$intercept["global.intercept"]
  ind_params <- 0
  Group_parameters <- list()
  for(g in 1:length(MEM_groups)){
    params <- NULL
    if(global_intercept){
      params <- c(params,Population_params[1])
      if(g == 1){
        ind_params <- ind_params + 1
      }
    }
    if(Marginal_dynamics$dynamic.type == "spline"){
      Nb_group_params <- as.numeric(1*Marginal_dynamics$intercept[paste("group.intercept",g,sep="")] + 
                                      length(Marginal_dynamics$knots[[MEM_groups[g]]]) + Marginal_dynamics$spline.degree[g])
    }else if(Marginal_dynamics$dynamic.type == "polynomial"){
      Nb_group_params <- as.numeric(1*Marginal_dynamics$intercept[paste("group.intercept",g,sep="")] + 
                                      Marginal_dynamics$polynomial.degree[g])
    }
    params <- c(params,Population_params[(ind_params+1):(ind_params+Nb_group_params)])
    ind_params <- ind_params + Nb_group_params
    Group_parameters[[MEM_groups[g]]] <- params
  }
  
  # Step 2: Calculation of AUC
  # ----- #
  Estimated_AUC <- NULL
  for(g in 1:length(Groups)){
    time_group <- time[[g]]
    beta_group <- Group_parameters[[Groups[g]]]
    Pop_Covariate <- NULL
    
    if(global_intercept){
      Pop_Covariate <- cbind(Pop_Covariate,rep(1,length(time_group)))
    }
    
    # Extraction of information about model
    if(Marginal_dynamics$dynamic.type == "polynomial"){
      # Creation of covariate matrix
      Covariate_poly_group <- do.call(cbind,lapply(1*isFALSE(Marginal_dynamics$intercept[paste("group.intercept",g,sep="")]):Marginal_dynamics$polynomial.degree[g],function(d) time_group^d))
      Pop_Covariate <- cbind(Pop_Covariate,Covariate_poly_group)
    }else if(Marginal_dynamics$dynamic.type == "spline"){
      # Creation of covariate matrix
      if(is.null(Marginal_dynamics$boundary.knots[[Groups[g]]])){
        Covariate_spline_group <- splines::bs(x=time_group,knots=Marginal_dynamics$knots[[Groups[g]]],df=Marginal_dynamics$df[g],degree=Marginal_dynamics$spline.degree[g])
      }else{
        Covariate_spline_group <- splines::bs(x=time_group,knots=Marginal_dynamics$knots[[Groups[g]]],df=Marginal_dynamics$df[g],degree=Marginal_dynamics$spline.degree[g],Boundary.knots=Marginal_dynamics$boundary.knots[[Groups[g]]])
      }
      if(Marginal_dynamics$intercept[paste("group.intercept",g,sep="")]){
        Covariate_spline_group <- cbind(rep(1,length(time_group)),Covariate_spline_group)
      }
      Pop_Covariate <- cbind(Pop_Covariate,Covariate_spline_group)
    }# End spline covariate
    
    # Estimation of the marginal dynamics
    Group_dynamics <- as.numeric(Pop_Covariate %*% beta_group)
    # Creation of method time weights (W) vector
    time_weights <- AUC_time_weights_estimation(time=time_group,method)
    AUC_group <- as.numeric(Group_dynamics %*% time_weights)
    Estimated_AUC <- c(Estimated_AUC,AUC_group)
  }
  
  names(Estimated_AUC) <- Groups
  if(Averaged){
    Estimated_nAUC <- sapply(seq(1,length(Groups)),function(g) Estimated_AUC[g]/diff(range(time[[g]])),simplify=TRUE)
    Results <- Estimated_nAUC
  }else{
    Results <- Estimated_AUC
  }
  return(Results)
}
