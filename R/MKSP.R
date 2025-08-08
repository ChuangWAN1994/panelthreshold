
#' @title Estimate the Multiple Kink effects for Static Panel (MKSP) models.
#'
#' @description This function can be used to estimate the multi-kink effects for static panel model using a profile two-stage method develop by Wan et al. (2023). The estimator is also an useful initial value for the minimization of the GMM estimation in MTDP with kink effects.
#' @param data data data a tibble of data. The number of rows of \code{data} must be the sample size \code{iT} times individuals number N.
#' @param iT the size for time.
#' @param var_dep the name of dependent variable.
#' @param var_thre the name of threshold variable.
#' @param var_base the names for the baseline covariates, i.e. \eqn{x_{it}}. Remark: for static panel model, \code{lag_y} can not be included.
#' @param im the number of kink effects.
#' @param trim_rate the trimming parameter. Default is 0.15.
#'
#' @return a list holding
#' \item{beta}{the estimator of the regression coefficients.}
#' \item{gamma}{the estimator of the threshold parameters.}
#' \item{beta.se}{the estimated standard error for \code{beta}.}
#' \item{gamma.se}{the estimated standard error for \code{gamma}}
#'
#'
#' @references
#' Wan et al. (2023). Multikink quantile regression for longitudinal data with application to progesterone data analysis. *Biometrics*
#' Zhang et al. (2017). Panel kink regression with an unknown threshold. *Economics Letter*
#' Wan and Li. (2025). Dynamic panel models with multi-threshold effects and endogeneity. *Statistics*
#'
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
#'   xx1 <- if(t>1) x1[,t-1] else 0
#'   x1[,t] <- rho1*xx1 + 0.25*alpha + rnorm(N,0,0.5)
#'   xx2 <- if(t>1) x2[,t-1] else 0
#'   x2[,t] <- rho2*xx2 + 0.5*alpha + rnorm(N,0,0.5)
#'   yy <- if(t>1) y[,t-1] else 0
#'   y[,t] <- x1[,t] + x2[,t] + alpha  - 2*(x2[,t] + 0.25)*(x2[,t]> -0.25) + 2*(x2[,t] - 0.25)*(x2[,t]>0.25) + rnorm(N,0,0.5)
#' }
#' y <- c(t(y))
#' x1 <- c(t(x1))
#' x2 <- c(t(x2))
#' data = cbind(y,x2,x1)
#' colnames(data) = c("y","x2","x1")
#' dep = "y"
#' thre = "x2"
#' var_base = c("x1","x2")
#' MKSP(data,iT,dep,thre,var_base,im=2,trim_rate=0.15)
#'
#' @export
MKSP <- function(data,iT,var_dep,var_thre,var_base,im,trim_rate=0.15)
{
  data <- as.matrix(data)
  N <- dim(data)[1]/iT
  n <- N*iT

  if(im <= 0){
    stop("number of threshold must be positive integer !")
  }

  if(length(var_dep)>1) stop(simpleError("Only one dependent variable!"))
  y = as.matrix(data[,var_dep])
  yt <- tr(y,N,iT)

  xx = as.matrix(data[,var_base])
  kx <- ncol(xx)

  if(any(apply(xx,2,function(x) all(x==1)))){
    stop("The regressor in the first regime can not contain the intercept !")
  }

  if(length(var_thre)>1) stop(simpleError("Only one threshold variable!"))
  kink.thre <- (data[,var_thre])
  kink.mat <- matrix(rep(kink.thre,im),ncol=im,byrow=FALSE)

  dd <- unique(kink.thre)
  cr <- consRegion(im,quantile(dd,trim_rate),quantile(dd,1-trim_rate))
  ui <- cr$ui; ci <- cr$ci

  regress <-function(y,x){
    if (qr(x)$rank == ncol(x)) beta <- qr.solve(x,y)
    if (qr(x)$rank < ncol(x)) beta <- (qr(t(x) %*% x)$qr) %*% t(x) %*% y
    beta
  }


  objctive_fun <- function(par){
    par.rep <- matrix(rep(par,n),ncol=im,byrow=TRUE)
    xxk <- cbind(xx,pmax((kink.mat-par.rep),0))
    xxkt <- tr(xxk, N, iT);
    e <- yt - xxkt %*% regress(yt,xxkt)
    out <- t(e)%*%e
    out
  }

  par.ini <- as.matrix(quantile(dd,seq(0,1,length=im+2)[-c(1,im+2)]))

  obj <- constrOptim(par.ini,f=objctive_fun,
                     grad=NULL, ui,
                     ci,method = "Nelder-Mead")

  rhat <- obj$par
  rhat.rep <- matrix(rep(rhat,n),ncol=im,byrow=TRUE)
  xxk <- cbind(xx,pmax((kink.mat-rhat.rep),0))
  xxkt <- tr(xxk, N, iT);
  beta <- regress(yt,xxkt)
  res <- as.matrix(yt - xxkt %*% beta)
  #res.mat <-

  ##standard error
  delta <- beta[(kx+1):(kx+im)]
  delta.mat <- matrix(rep(delta,n),nrow=n,ncol=im,byrow=TRUE)
  idx.mat <- as.matrix(kink.mat > rhat.rep)
  idx.mat.t <- as.matrix(tr(idx.mat,N,iT))
  gdev <- cbind(xxkt,-delta.mat * idx.mat.t)
  Dn <- matrix(0,kx+2*im,kx+2*im)
  dvec1 <- apply(idx.mat.t * (res %*% matrix(1,nrow=1,ncol=im)),2,sum)
  dvec <- c(rep(0,kx+im),dvec1)
  Dn <- diag(dvec)

  Hn <- t(gdev)%*%gdev + Dn
  Hn.inv <- solve(Hn)
  Gn <- gdev*(res %*% matrix(1,nrow=1,ncol=ncol(gdev)))
  Sig <- t(Gn) %*% Gn
  covmat <- Hn.inv %*% Sig %*% Hn.inv

  coef_names <- c(paste0("base-",var_base))
  for(kki in 1:im){
    coef_names <- c(coef_names,paste0("kink_slop-",kki,kki+1))
  }

  beta = as.vector(beta)
  theta.se <- sqrt(diag(covmat))
  beta.se <- theta.se[1:(kx+im)]
  gamma.se <- theta.se[(kx+im+1):(kx+2*im)]
  names(beta.se) = names(beta) = coef_names

  gamma <- drop(rhat)
  kink_names <- c()
  for(kki in 1:im){
    kink_names <- c(kink_names,paste0("kink_slop-",kki,kki+1))
  }

  names(gamma) = names(gamma.se) = kink_names

  eval = list()
  eval$beta = beta;
  eval$gamma = gamma
  eval$beta.se = beta.se
  eval$gamma.se = gamma.se
  return(eval)


}

