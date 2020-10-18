#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param MEM_Pol_group PARAM_DESCRIPTION
#' @param Group1 PARAM_DESCRIPTION
#' @param Group2 PARAM_DESCRIPTION
#' @param time.G1 PARAM_DESCRIPTION
#' @param time.G2 PARAM_DESCRIPTION
#' @param method PARAM_DESCRIPTION, Default: 'trapezoid'
#' @param Group.dependence PARAM_DESCRIPTION, Default: TRUE
#' @param Averaged PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[ArgumentCheck]{addError}}
#'  \code{\link[splines]{bs}}
#' @rdname Group_specific_Var_Delta_AUC_estimation
#' @export 
#' @importFrom ArgumentCheck newArgCheck addError finishArgCheck
#' @importFrom splines bs
Group_specific_Var_Delta_AUC_estimation <- function(MEM_Pol_group,Group1,Group2,time.G1,time.G2,method="trapezoid",Group.dependence=TRUE,Averaged=FALSE){
  
  # Step 1: Verification of the type of arguments
  # ----- #
  # Verification of 'Group.dependence'
  Check_group.dep <- ArgumentCheck::newArgCheck()
  if(isFALSE(is.logical(Group.dependence))){
    ArgumentCheck::addError(
      msg = "The argument 'Group.dependence' must be defined as boolean variable",
      argcheck = Check_group.dep
    )
  }
  ArgumentCheck::finishArgCheck(Check_group.dep)
  # Verification of Group1 and Group2 not null
  Check_groups_null <-  ArgumentCheck::newArgCheck()
  if(is.null(Group1) || is.null(Group2)){
    ArgumentCheck::addError(
      msg = "One of the two Groups to compared has been assigned to 'NULL' value",
      argcheck = Check_groups_null
    )
  }
  ArgumentCheck::finishArgCheck(Check_groups_null)
  
  Groups <- c(Group1,Group2)
  time <- list(time.G1,time.G2)
  Check_argument_Group_specific_Var_AUC(MEM_Pol_group,time,Groups,method,Averaged)
  
  
  Model_features <- MEM_Pol_group$Model_features
  Marginal_dynamics <- Model_features$Marginal.dyn.feature
  
  # Extraction of population parameters according to their groups
  Population_params <- MEM_Pol_group$Model_estimation$beta
  Population_variance <- MEM_Pol_group$Model_estimation$varFix
  MEM_groups <- Model_features$Groups
  
  global_intercept <- Marginal_dynamics$intercept["global.intercept"]
  ind_params <- 0
  Group_parameters <- list()
  Group_index_params <- list()
  for(g in 1:length(MEM_groups)){
    params <- NULL
    index_params <- NULL
    if(global_intercept){
      params <- c(params,Population_params[1])
      index_params <- c(index_params,1)  
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
    index_params <- c(index_params,seq(ind_params+1,ind_params+Nb_group_params))
    # variance_matrix <- Population_variance[index_params,index_params]
    ind_params <- ind_params + Nb_group_params
    Group_parameters[[MEM_groups[g]]] <- params
    Group_index_params[[MEM_groups[g]]] <- index_params
  }
  Groups_index_params <- unique(sort(as.numeric(unlist(Group_index_params))))
  Groups_Variance <- Population_variance[Groups_index_params,Groups_index_params]
  
  # Step 2: Calculation of the variance of Delta AUC
  # ----- #
  if(Group.dependence){
    Combined_Pop_Covariate <- NULL
    Pop_Covariate <- list()
    for(g in 1:length(Groups)){
      time_group <- time[[g]]
      
      if(Marginal_dynamics$dynamic.type == "polynomial"){
        # Creation of covariate matrix
        Covariate_group <- do.call(cbind,lapply(1*isFALSE(Marginal_dynamics$intercept[paste("group.intercept",g,sep="")]):Marginal_dynamics$polynomial.degree[g],function(d) time_group^d))
      }else if(Marginal_dynamics$dynamic.type == "spline"){
        # Creation of covariate matrix
        if(is.null(Marginal_dynamics$boundary.knots[[Groups[g]]])){
          Covariate_group <- splines::bs(x=time_group,knots=Marginal_dynamics$knots[[Groups[g]]],df=Marginal_dynamics$df[g],degree=Marginal_dynamics$spline.degree[g])
        }else{
          Covariate_group <- splines::bs(x=time_group,knots=Marginal_dynamics$knots[[Groups[g]]],df=Marginal_dynamics$df[g],degree=Marginal_dynamics$spline.degree[g],Boundary.knots=Marginal_dynamics$boundary.knots[[Groups[g]]])
        }
        if(Marginal_dynamics$intercept[paste("group.intercept",g,sep="")]){
          Covariate_group <- cbind(rep(1,length(time_group)),Covariate_group)
        }
      } # End spline covariate
      Pop_Covariate[[Groups[g]]] <- Covariate_group
    }
    Combined_Pop_Covariate <- matrix(0,ncol=ncol(Pop_Covariate[[1]]) + ncol(Pop_Covariate[[2]]),
                                     nrow=nrow(Pop_Covariate[[1]]) + nrow(Pop_Covariate[[2]]))
    Combined_Pop_Covariate[1:nrow(Pop_Covariate[[1]]),1:ncol(Pop_Covariate[[1]])] <- Pop_Covariate[[1]]
    Combined_Pop_Covariate[(nrow(Pop_Covariate[[1]])+1):nrow(Combined_Pop_Covariate),(ncol(Pop_Covariate[[1]])+1):ncol(Combined_Pop_Covariate)] <- Pop_Covariate[[2]]
    if(global_intercept){
      Combined_Pop_Covariate <- cbind(rep(1,nrow(Combined_Pop_Covariate)),Combined_Pop_Covariate)
    }
    
    time_weights.G1 <- AUC_time_weights_estimation(time=time[[1]],method)
    time_weights.G2 <- AUC_time_weights_estimation(time=time[[2]],method)
    if(Averaged){
      Global_time_weights <- c(rep(0,length(time_weights.G1)),time_weights.G2)/(diff(range(time.G2))) -c(time_weights.G1,rep(0,length(time_weights.G2)))/(diff(range(time.G1)))
    }else{
      Global_time_weights <- c(rep(0,length(time_weights.G1)),time_weights.G2) -c(time_weights.G1,rep(0,length(time_weights.G2)))
    }
    Var_Delta_AUC <- as.numeric(Global_time_weights %*% Combined_Pop_Covariate %*% Groups_Variance %*% t(Combined_Pop_Covariate) %*% Global_time_weights)
  }else{
    Var_AUC_Group1 <- Group_specific_Var_AUC_estimation(MEM_Pol_group=MEM_Pol_group,time=time.G1,Groups=Group1,method=method,Averaged=Averaged)
    Var_AUC_Group2 <- Group_specific_Var_AUC_estimation(MEM_Pol_group=MEM_Pol_group,time=time.G2,Groups=Group2,method=method,Averaged=Averaged)
    Var_Delta_AUC <- Var_AUC_Group1 + Var_AUC_Group2
  }
  return(Var_Delta_AUC)
}
