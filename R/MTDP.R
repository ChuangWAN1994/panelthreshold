

#' @title: Estimate the Multiple Threshold for Dynamic Panel models.
#'
#' @description  The function needs the return value (an object of the class MTDP) from the \code{\link{NewMTDP}}. It copies the object, reuses its contents to estimate the corresponding MTDP model, and then returns a new object of the class MTDP containing the results from the estimation. The user can choose to save the return value to a new object or simply to overwrite the object returned from \code{NewMTDP}.
#'
#' @param use an object of the class MTDP, created by \code{\link{NewMTDP}} function.
#' @param im the number of thresholds. The default value is 1.
#' @param h_0 the constant of the bandwidth for the kernel function, i.e., \eqn{h=h_0\hat{\sigma}n^{-1/5}}. This is mainly used in estimating the variance-covariance matrix of the estimators. Here \code{h_0} can be a vector, and its range is generally [0.1,2). We can set multiple values, and multiple estimated standard deviations will also be output. Default is 1.
#' @param par_ini the initial values for the parameters used in the "Nelder-Mead" algorithm. If NULL, the quantiles of threshold variable is specified.
#' @param trim_rate the trimming parameter. Default is 0.15.
#' @param avg_num the number of randomizations in constructing the weight matrix at the first step of GMM estimation. Default is 0.
#' @param nboots the number of bootstrap repetitions in model specification test. If 0,  no test will be performed, and the output p-value will be -1. Default is 0.
#' @param grid_num number of grid points in setting the searching domain in the model specification test. It will be used only when \code{nboots} is greater than 0. Default is 30.
#'
#' @return a new object of the class MTDP containing the results from the estimation.
#'
#' The object is a list containing the components:
#'
#' \item{J_hat}{the objective value.}
#' \item{est}{a vector of estimates of all unknown parameters.}
#' \item{se}{a vector of the standard errors of all the estimates.}
#' \item{ep_collect}{the estimated residuals.}
#' \item{L}{the number of instruments}
#' \item{boots_p}{the p value for the model specification test.}
#'
#' @references
#' Chuang Wan and Yi Li. (2025). Dynamic panel models with multi-threshold effects and endogeneity. *Statistics*
#'
#' @seealso \code{\link{NewMTDP}}
#'
#' @export
#'
#' @examples
#' N <- 100
#' iT <- 10
#' y <- matrix(NA,N,iT)
#' x1 <- x2 <- matrix(NA,N,iT)
#' alpha <- rnorm(N,0,1)
#' rho1 <- runif(N,-0.9,0.9)
#' rho2 <- runif(N,-0.9,0.9)
#' for(t in 1:iT){
#'  xx1 <- if(t>1) x1[,t-1] else 0
#'  x1[,t] <- rho1*xx1 + 0.25*alpha + rnorm(N,0,0.5)
#'  xx2 <- if(t>1) x2[,t-1] else 0
#'  x2[,t] <- rho2*xx2 + 0.5*alpha + rnorm(N,0,0.5)
#'  yy <- if(t>1) y[,t-1] else 0
#'  y[,t] <- x1[,t] + alpha  + (1-2*x1[,t])*(x2[,t]> -0.25) + (-1+ 2* x1[,t])*(x2[,t]>0.25) + 0.2*yy + rnorm(N,0,0.5)
#' }
#' y <- c(t(y))
#' x1 <- c(t(x1))
#' x2 <- c(t(x2))
#' data = cbind(y,x2,x1,1)
#' colnames(data) = c("y","x2","x1","cons")
#' dep = "y"
#' thre = "x2"
#' var_rx = c("cons","x1")
#' var_base = c("lag_y","x1")
#' var_endo = NULL
#' var_iv = NULL
#' use = NewMTDP(data,thre_type=c("jump"),dep,thre,var_base,var_rx,var_endo,var_iv,iT)
#' im = 2
#' h_0 = c(0.5)
#' par = NULL
#' trim_rate=0.15;
#' grid_num=100;
#' avg_num=0
#' obj_mtdp <- MTDP(use,im,h_0,par,trim_rate,avg_num,nboots = 0 )
#' obj_mtdp$est
#'
#'
#'
MTDP <- function(use,
                 im=1,
                 h_0=1,
                 par_ini=NULL,
                 trim_rate=0.15,
                 avg_num=0,
                 nboots = 0,
                 grid_num=30)
{
  if(class(use)!="MTDP")
    stop(simpleError("The argument 'use' is not an object of class 'PSTR'"))
  ret = use

  # get the data here
  y_mat_rep = ret$y_mat_rep; x_temp1 = ret$x_temp1;
  rx_temp2 = ret$rx_temp2;rx_temp3 = ret$rx_temp3
  q_temp1 = ret$q_temp1; q_temp2 = ret$q_temp2
  iT = ret$iT; N = ret$N; flag_static = ret$flag_static;
  q_mat = ret$q_mat; Z_mat = ret$Z_mat
  moment_mat = ret$moment_mat
  var_rx = ret$var_rx_name; var_base = ret$var_base_name
  thre_type = ret$thre_type

  if(flag_static == 1) {
    W_n = GMM_W_n_static(Z_mat, moment_mat, iT)
    L = sum(moment_mat) - sum(moment_mat[, 2])
    end_p = ncol(x_temp1)
    x_temp1 = x_temp1[, (N+1):end_p]
    rx_temp2 = rx_temp2
    rx_temp3 = rx_temp3

  }else{
    W_n = GMM_W_n_con(Z_mat, moment_mat, iT)
    L = sum(moment_mat) - sum(moment_mat[1,])
  }

  g_1n_bar = rowSums(Z_mat*y_mat_rep)/N

  if(im == 0){
    #// First step GMM estimation
    k_x <- ncol(x_temp1)/N
    g_2n_bar = matrix(0,L,(k_x))

    for(kk3 in 1:k_x){
      g_2n_bar[,kk3] = rowSums(Z_mat*x_temp1[,((kk3-1)*N+1):(kk3*N)])/N
    }
    library(MASS)
    GMM_est = ginv(t(g_2n_bar)%*% W_n %*% g_2n_bar) %*% (t(g_2n_bar)%*%W_n%*%g_1n_bar)

    #
    ##// calculate the error for linear model
    ep_collect <- matrix(0,L,N)
    ind_vec = seq(0,(k_x-1)*N, length= k_x)

    for(ii2 in 1:L){
      for(ii3 in 1:N){
        ii3_rep = matrix(ii3,k_x,1)
        xreg = c(x_temp1[ii2,(ind_vec+ii3_rep)]) #(x_temp1[ii2,ind_vec+ii3_rep])
        ep_collect[ii2,ii3] = y_mat_rep[ii2,ii3] - t(GMM_est) %*% matrix(xreg,ncol=1)
      }
    }

    W_n_2 = GMM_W_n_2_con(ep_collect, Z_mat)
    GMM_est_2 = ginv(t(g_2n_bar)%*% W_n_2 %*% g_2n_bar) %*% (t(g_2n_bar)%*%W_n_2%*%g_1n_bar)

    ep_collect_2 <- matrix(0,L,N)
    for(ii2 in 1:L){
      for(ii3 in 1:N){
        ii3_rep = matrix(ii3,k_x,1)
        xreg = c(x_temp1[ii2,(ind_vec+ii3_rep)]) #(x_temp1[ii2,ind_vec+ii3_rep])
        ep_collect_2[ii2,ii3] = y_mat_rep[ii2,ii3] - t(GMM_est_2) %*% matrix(xreg,ncol=1)
      }
    }

    g_n_bar = g_1n_bar - g_2n_bar %*% GMM_est_2
    J_result <- t(g_n_bar) %*% W_n_2 %*% g_n_bar  #W_n_2

    eval = list()
    eval$ep_collect = ep_collect_2
    eval$J_hat = drop(J_result);
    eval$est = as.vector(GMM_est_2)
    eval$L = L
    return(eval)

  }

  q_var <- as.vector(q_mat)
  q_var <- unique(q_var)
  grid_lower = quantile(q_var,trim_rate)
  grid_upper = quantile(q_var,1-trim_rate)

  if(is.null(par_ini) | missing(par_ini)){
    par_ini <- as.matrix(quantile(q_var,seq(0,1,length=im+2)[-c(1,im+2)]))
  }else{
    par_ini <- as.matrix(par_ini)
  }

  cr <- consRegion(im,grid_lower,grid_upper)
  ui <- cr$ui; ci <- cr$ci
  if(thre_type == "jump"){
    NM <- try(MMoptim_jump(ui,ci,par_ini,y_mat_rep,g_1n_bar,Z_mat,
                           x_temp1,rx_temp2,rx_temp3,
                           q_temp1,q_temp2,W_n),silent = TRUE)
  }else if(thre_type == "kink"){
    NM <- try(MMoptim_kink(ui,ci,par_ini,y_mat_rep,g_1n_bar,Z_mat,
                           x_temp1,rx_temp2,rx_temp3,
                           q_temp1,q_temp2,W_n),silent=TRUE)
  }

  if(is(NM,"try-error")) stop(simpleError("The two-stage profile estimation fails !"))
  par_est <- as.matrix(NM$par_est); GMM_est <- as.matrix(NM$GMM_est)
  ep_collect <- NM$ep_collect; J_result <- NM$J_result


  if(avg_num > 0){
    set.seed(123)
    par_est_avg = cbind(par_est)

    for(ii in 1:avg_num){
      Weight_1 <- GMM_weight(N,iT,L,moment_mat,flag_static)
      W_n_avg <- GMM_W_n_2_con(Weight_1,Z_mat)
      if(thre_type == "jump"){
        obj_avg <- try(MMoptim_jump(ui,ci,par_ini,y_mat_rep,g_1n_bar,Z_mat,
                                    x_temp1,rx_temp2,rx_temp3,
                                    q_temp1,q_temp2,W_n_avg),silent = TRUE)
      }else if(thre_type == "kink"){
        obj_avg <- try(MMoptim_kink(ui,ci,par_ini,y_mat_rep,g_1n_bar,Z_mat,
                                    x_temp1,rx_temp2,rx_temp3,
                                    q_temp1,q_temp2,W_n_avg),silent=TRUE)
      }

      if(is(obj_avg,"try-error")) next()
      par_est_avg_1 <- as.matrix(obj_avg$par_est)
      par_est_avg = cbind(par_est_avg,par_est_avg_1)
    }
    par_est <- as.matrix(apply(par_est_avg,1,function(x) mean(x,na.rm = TRUE)))

    if(thre_type == "jump"){
      obj <- jump_GMM_est_Fun(y_mat_rep,g_1n_bar,par_est,Z_mat,x_temp1,rx_temp2,rx_temp3,
                              q_temp1,q_temp2,W_n)
    }else if(thre_type == "kink"){
      obj <- kink_GMM_est_Fun(y_mat_rep,g_1n_bar,par_est,Z_mat,x_temp1,rx_temp2,rx_temp3,
                              q_temp1,q_temp2,W_n)
    }

    ep_collect_1 = obj$ep_collect;
    W_n_2 = GMM_W_n_2_con(ep_collect_1, Z_mat)
    if(thre_type == "jump"){
      obj2 <- jump_GMM_est_Fun(y_mat_rep,g_1n_bar,par_est,Z_mat,x_temp1,rx_temp2,rx_temp3,
                               q_temp1,q_temp2,W_n_2)
    }else if(thre_type == "kink"){
      obj2 <- kink_GMM_est_Fun(y_mat_rep,g_1n_bar,par_est,Z_mat,x_temp1,rx_temp2,rx_temp3,
                               q_temp1,q_temp2,W_n_2)
    }
    ep_collect = obj2$ep_collect; GMM_est = obj2$GMM_est; J_result = obj2$J_result

  }


  est <-  c(GMM_est , par_est)
  ke <- length(est)
  omega_hat = GMM_W_n_2_con(ep_collect, Z_mat)


  if(thre_type == "jump"){
    k_h <- length(h_0)
    sds <- matrix(NA,k_h,ke)
    for(ik in 1:k_h){
      cov_mat = cov_mat_con_Multithreshold(GMM_est,par_est,q_mat,Z_mat,x_temp1,rx_temp2,
                                           rx_temp3, q_temp1,q_temp2,h_0[ik],omega_hat)

      sds[ik,] = c((diag(cov_mat))^(1/2))
    }
  }else if(thre_type == "kink"){
    sds <- matrix(NA,1,ke)
    cov_mat = cov_mat_con_Multikink(GMM_est,par_est,q_mat,Z_mat,x_temp1,rx_temp2,
                                    rx_temp3, q_temp1,q_temp2,omega_hat)

    sds[1,] = c((diag(cov_mat))^(1/2))
  }


  coef_names <- c(var_base)
  for(kki in 1:im){
    coef_names <- c(coef_names,paste0("R",kki,kki+1,"-",var_rx))
  }
  coef_names <- c(coef_names, paste0("thre-",1:im))
  colnames(sds) = names(est) = coef_names
  if(thre_type == "jump"){
    rownames(sds) = paste0("jump-h_0-",h_0)
  }else if(thre_type == "kink"){
    rownames(sds) = paste0("kink")
  }


  eval = list()
  eval$J_hat = drop(J_result);
  eval$est = est
  eval$se = sds
  eval$ep_collect = ep_collect
  eval$L = L

  #//...........................Test..............................\\

  if(nboots > 0){


    grid_th = quantile(q_var,seq(trim_rate,1-trim_rate,length=grid_num))

    Gmat <- set_grid(im,grid_num,grid_th,trim_rate)
    if(thre_type == "jump"){
      ob <- sup_wald_cal_multithreshold(y_mat_rep,Z_mat,x_temp1,rx_temp2,rx_temp3,
                                        q_temp1,q_temp2,im,Gmat,N,W_n)
      sup_wald_ori <- ob$sup_wald
      cov_wald_rec <- ob$cov_wald_rec
      wald_rec = rep(NA,nboots) #To save sup-wald's
      for(bb2 in 1:nboots){

        sup_wald_b = sup_wald_fb_multithreshold(y_mat_rep,Z_mat,x_temp1,rx_temp2,rx_temp3,
                                                q_temp1,q_temp2,ep_collect,omega_hat,cov_wald_rec,
                                                Gmat,N)
        wald_rec[bb2] = sup_wald_b
      }

    }else if(thre_type == "kink"){
      ob <- sup_wald_cal_multikink(y_mat_rep,Z_mat,x_temp1,rx_temp2,
                                   rx_temp3,q_temp1,q_temp2,im,Gmat,
                                   N,W_n)
      sup_wald_ori <- ob$sup_wald
      cov_wald_rec <- ob$cov_wald_rec
      wald_rec = rep(NA,nboots) #To save sup-wald's
      for(bb2 in 1:nboots){

        sup_wald_b = sup_wald_fb_multikink(y_mat_rep,Z_mat,x_temp1,rx_temp2,rx_temp3,
                                           q_temp1,q_temp2,ep_collect,omega_hat,cov_wald_rec,
                                           Gmat,N)
        wald_rec[bb2] = sup_wald_b
      }

    }

    boots_p = mean(wald_rec>sup_wald_ori,na.rm=TRUE)

  }else{
    boots_p = -1
  }

  eval$boots_p = boots_p

  return(eval)
}
