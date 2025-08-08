
#' @title Create an object of the class MTDP (Multiple Threshold for Dynamic Panel models).
#'
#' @param data a tibble of data. The number of rows of \code{data} must be the sample size \code{iT} times individuals number N.
#' @param thre_type the threshold regression type. Options are
#' \describe{
#' \item{\code{"jump"}}{for estimating the threshold regression with jump effects.}
#' \item{\code{"kink"}}{for estimating the kink regression with continuous threshold effects.}
#' } Default is \code{"jump"}.
#' @param var_dep the name of dependent variable.
#' @param var_thre the name of threshold variable.
#' @param var_base the names for the baseline covariates, i.e. \eqn{x_{it}}. Remark: for dynamic model, \code{lag_y} should be included.
#' @param var_rx the names of the regime-dependent variables; the \code{lag_y} and the \code{intercept} are allowed, but need to be specified.
#' @param var_endo the names for all endogenous variables including the regressors and threshold variables.
#' @param var_iv the names for additional instrument variables. Default is NULL.
#' @param iT the size for time.
#'
#' @return An object of the class MTDP for later usage in \code{\link{MTDP}}.
#'
#' @references
#' Wan and Li. (2025). Dynamic panel models with multi-threshold effects and endogeneity. *Statistics*
#'
#'@seealso \code{\link{MTDP}}
#'
#' @examples
#'
#' N <- 100
#' iT <- 10
#' y <- matrix(NA,N,iT)
#' x1 <- x2 <- matrix(NA,N,iT)
#' alpha <- rnorm(N,0,1)
#' rho1 <- runif(N,-0.9,0.9)
#' rho2 <- runif(N,-0.9,0.9)
#' for(t in 1:iT){
#'   xx1 <- if(t>1) x1[,t-1] else 0
#'  x1[,t] <- rho1*xx1 + 0.25*alpha + rnorm(N,0,0.5)
#'  xx2 <- if(t>1) x2[,t-1] else 0
#'  x2[,t] <- rho2*xx2 + 0.5*alpha + rnorm(N,0,0.5)
#'  yy <- if(t>1) y[,t-1] else 0
#'  y[,t] <- x1[,t] + alpha  + (1-2*x1[,t])*(x2[,t]> -0.25) + (-1+ 2* x1[,t])*(x2[,t]>0.25) + 0.2*yy + rnorm(N,0,0.5)
#'}
#'y <- c(t(y))
#'x1 <- c(t(x1))
#'x2 <- c(t(x2))
#'data = cbind(y,x2,x1,1)
#'colnames(data) = c("y","x2","x1","cons")
#'dep = "y"
#'thre = "x2"
#'var_rx = c("cons","x1")
#'var_base = c("lag_y","x1")
#'var_endo = NULL
#'var_iv = NULL
#'use = NewMTDP(data,dep,thre,var_base,var_rx, var_endo,var_iv,iT)
#'
#'
#' @export
#'
NewMTDP <- function(data,thre_type=c("jump"),var_dep,var_thre,var_base,var_rx,
                    var_endo,var_iv=NULL,iT)
{
  #if(!is_tibble(data)) stop(simpleError("data should be a tibble!"))
  ret = list(); class(ret) = "MTDP"
  N <- dim(data)[1]/iT


  if(length(var_dep)>1) stop(simpleError("Only one dependent variable!"))
  y_mat = data[,var_dep]
  y_mat = matrix(y_mat,ncol=iT,byrow=TRUE)

  if(var_dep %in% var_endo) stop(simpleError("The reponse variable cannot be endogeneous!"))


  k_x <- length(var_base)
  if("lag_y" %in% var_base){
    cat("The current estimation is performed for dynamic panel model! \n")
    var_base = setdiff(var_base,"lag_y")
    flag_static = 0
    var_base_name = c("lag_y")
  }else{
    if("lag_y" %in% var_rx) stop(simpleError("The lagged y should not be contained in the regime-dependent variable for static panel !"))
    cat("The current estimation is performed for static panel model! \n")
    if(length(var_base) < 1) stop(simpleError("The baseline variables for static panel should contain some informative covarites!"))
    flag_static = 1
    var_base_name = c()
  }
  var_base_name = c(var_base_name,var_base)


  #// define the matrix for the baseline covariates
  if(any(!var_base %in% colnames(data))) stop("The dependent variables do not belong to the data set.")
  lag_y_mat = cbind(matrix(1,N,1),y_mat[,1:(iT-1)])
  x_mat = lag_y_mat
  data_x <- as.matrix(data[,var_base])
  if(any(apply(data_x,2,function(x) all(x==1)))){
    stop("The intercept in the first regime can not be estimated !")
  }

  k_b <- length(var_base)
  if(k_b >= 1){
    for(jj in 1:(k_b)){
      x_temp = data_x[,jj]
      x_temp = matrix(x_temp,ncol=iT,byrow=TRUE)
      x_mat = rbind(x_mat,x_temp)
    }
  }

  if(length(var_thre)>1) stop(simpleError("Only one threshold variable!"))
  if(var_thre == "lag_y"){ #The special case that threshold variable is the lagged $y$.
    q_mat = lag_y_mat
  }else{
    q_mat = data[,var_thre]
    q_mat = matrix(q_mat,ncol=iT,byrow=TRUE)
  }


  #/*.........regime-dependent variables.........*\

  rx_mat = NULL
  var_rx_name = c()
  #if("cons" %in% var_rx){
  #  rx_mat = matrix(1,N,T)
  #  var_rx = setdiff(var_rx,"cons")
  #  var_rx_name = c(var_rx_name,"cons")
  #}
  if(thre_type == "kink"){
    rx_mat = q_mat
    var_rx_name = var_thre
    cat("Kink threshold regression; (the regime-dependent variable should be the threshold variable) !\n")
  }else if(thre_type == "jump"){
    cat("Jump threshold regression ! \n")
    if("lag_y" %in% var_rx){
      rx_mat = rbind(rx_mat,lag_y_mat)
      var_rx = setdiff(var_rx,"lag_y")
      var_rx_name = c(var_rx_name,"lag_y")
    }
    if(length(var_rx) != 0){
      var_rx_name = c(var_rx_name,var_rx)
      if(any(!var_rx %in% colnames(data))) stop("The regime-dependent variables do not belong to the data set.")
      data_rx <- as.matrix(data[,var_rx])
      for(jj in 1:(ncol(data_rx))){
        rx_temp = data_rx[,jj]
        rx_temp = matrix(rx_temp,ncol=iT,byrow=TRUE)
        rx_mat = rbind(rx_mat,rx_temp)
      }
    }
  }


  k_rx = nrow(rx_mat)/N

  #/*.........instrumental variables.........*\
  k_inst = length(var_iv)
  if(k_inst >= 1){
    if(any(!var_iv %in% colnames(data))) stop("The instrument variables do not belong to the data set.")
    data_IV = as.matrix(data[,var_iv])
    if(k_inst >= 2){
      #when there exist 2 or more additional IV's
      IV_mat = matrix(data_IV[,1],ncol=iT,byrow = TRUE)
      for(jj2 in 2:k_inst){
        #//For-loop constructing (Nk_inst * T) IV_mat
        IV_temp = data_IV[,jj2]
        IV_temp = matrix(IV_temp,ncol=iT,byrow=TRUE)
        IV_mat = rbind(IV_mat,IV_temp)
      }
    }else{
      #// When there exists only 1 IV.
      IV_mat = matrix(data_IV[,1],ncol=iT,byrow=TRUE)
    }
  }else{
    IV_mat = NULL;k_inst =0
  }

  y_mat_fd = FD_con(y_mat,iT) #// This will be (N*(T-1)) matrix
  x_mat_fd = FD_con(x_mat,iT)


  #//#### define the IV i.e. variable Z in the paper ####//
  var_base1 = var_base
  if((!(var_thre %in% var_base)) & (var_thre != "lag_y")){
    var_base1 = c(var_thre,var_base)
  }
  k_endo = length(var_endo)
  if(k_endo > 0){
    var_iv_list = setdiff(var_base1,var_endo)
    var_iv_list <- c(var_iv_list,var_endo)
  }else{
    var_iv_list <- var_base1
  }
  k_exo = 0
  x_iv_mat = lag_y_mat
  data_iv_x <- as.matrix(data[,var_iv_list])
  k_b1 <- length(var_base1)
  if(k_b1 >= 1){
    for(jj in 1:(k_b1)){
      x_temp = data_iv_x[,jj]
      x_temp = matrix(x_temp,ncol=iT,byrow=TRUE)
      x_iv_mat = rbind(x_iv_mat,x_temp)
    }
  }
  k_1 <- nrow(x_iv_mat)/N
  k_endo_n = k_endo + k_exo
  moment_mat = moment_mat_con(iT, k_1, k_endo_n)
  mt_c = ncol(moment_mat); mt_r = nrow(moment_mat);
  x_iv_mat_fd = FD_con(x_iv_mat,iT)
  moment_mat = moment_mat + cbind(matrix(0,mt_r, mt_c - k_exo), matrix(2,mt_r, k_exo))

  Z_mat = Z_mat_con(y_mat, x_iv_mat, x_iv_mat_fd, moment_mat, k_endo_n,flag_static)
  if(k_inst >= 1) {
    ziv = Z_mat_IV_add(Z_mat, moment_mat, IV_mat, k_inst,flag_static)
    moment_mat = ziv$moment_mat; Z_mat = ziv$Z_mat
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #                   The rx regressor process
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Gvc = GMM_var_con(y_mat_fd,x_mat_fd,rx_mat,q_mat,moment_mat,N,iT,flag_static)
  ret$y_mat_rep = Gvc$y_mat_rep; ret$x_temp1 = Gvc$x_temp1;
  ret$rx_temp2 = Gvc$rx_temp2;ret$rx_temp3 = Gvc$rx_temp3
  ret$q_temp1 = Gvc$q_temp1; ret$q_temp2 = Gvc$q_temp2

  ret$iT = iT; ret$N = N; ret$flag_static = flag_static;
  ret$q_mat = q_mat; ret$Z_mat = Z_mat;
  ret$moment_mat=moment_mat;
  ret$var_rx_name = var_rx_name; ret$var_base_name  = var_base_name
  ret$thre_type = thre_type

  return(ret)
}
