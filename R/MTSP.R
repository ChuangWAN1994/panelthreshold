
#' @title Estimate the Multiple Threshold for Static Panel (MTSP) models.
#'
#' @description This function mainly implements the sequential estimation procedure proposed by Hansen (1999) designed for static panel model with threshold effects. The estimator obtained by the static panel model is also a useful initial value for the minimization of the GMM estimation in MTDP.
#'
#' @param data data a tibble of data. The number of rows of \code{data} must be the sample size \code{iT} times individuals number N.
#' @param iT the size for time.
#' @param var_dep the name of dependent variable.
#' @param var_thre the name of threshold variable.
#' @param var_base the names for the baseline covariates, i.e. \eqn{x_{it}}. Remark: for static panel model, \code{lag_y} can not be included.
#' @param var_rx the names of the regime-dependent variables. For static panel models, \code{lag_y} is not allowed.
#' @param im the number of thresholds  (maximum 4). Default value is 1.
#' @param trim_rate the trimming parameter. Default is 0.15.
#' @param grid_num number of grid points for the grid-search method. Default is 100.
#'
#' @return a list holding
#' \item{beta}{the estimator of the regression coefficients.}
#' \item{gamma}{the estimator of the threshold parameters.}
#' \item{beta.homo.se}{the estimated standard error for \code{beta} with homogeneous error assumption.}
#' \item{beta.het.se}{the estimated standard error for \code{beta} with heterogeneous error assumption.}
#'
#' @references
#' Hansen. (1999). Threshold effects in non-dynamic panels: Estimation, testing, and inference. *Journal of Econometrics*
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
#'   y[,t] <- x1[,t] + alpha  + (1-2*x1[,t])*(x2[,t]> -0.25) + (-1+ 2* x1[,t])*(x2[,t]>0.25) + rnorm(N,0,0.5)
#' }
#' y <- c(t(y))
#' x1 <- c(t(x1))
#' x2 <- c(t(x2))
#' data = cbind(y,x2,x1,1)
#' colnames(data) = c("y","x2","x1","cons")
#' dep = "y"
#' thre = "x2"
#' var_rx = c("cons","x1")
#' var_base = c("x1")
#' obj_mtsp <- MTSP(data,iT,dep,thre,var_base,var_rx,im=2)
#' obj_mtsp$beta
#' obj_mtsp$gamma
#'
#'
#' @export
MTSP <- function(data,iT,var_dep,var_thre,var_base,var_rx,im=1,
                 trim_rate=0.15,grid_num=100)
{
  data <- as.matrix(data)
  N <- dim(data)[1]/iT
  #------------------------------------------------------------------#
  # The number of thresholds should be less 3                        #
  #------------------------------------------------------------------#

  if(im <= 0 || im >3){
    stop("number of threshold must be positive integer and less than 3!")
  }

  if(length(var_dep)>1) stop(simpleError("Only one dependent variable!"))
  y = as.matrix(data[,var_dep])
  yt <- tr(y,N,iT)

  regx = as.matrix(data[,var_rx])
  xt = as.matrix(data[,var_base])

  if(any(apply(xt,2,function(x) all(x==1)))){
    stop("The regressor in the first regime can not contain the intercept !")
  }


  if(length(var_thre)>1) stop(simpleError("Only one threshold variable!"))
  thresh <- as.matrix(data[,var_thre])
  ty <- N*(iT-1)
  dd <- unique(thresh)
  dd <- as.matrix(sort(dd))
  qnt1 <- grid_num*trim_rate
  sq <- as.matrix(seq(trim_rate,trim_rate+(1/grid_num)*(grid_num-2*qnt1),by=1/grid_num))
  qq1 <- as.matrix(dd[floor(sq*nrow(dd))])
  qn1 <- nrow(qq1)

  #------------------------------------------------------------------#
  # The function for calculating the sum of squared errors           #
  #------------------------------------------------------------------#

  sse_calc <- function(y,x){
    e <- y-x%*%qr.solve(x,y)
    out <- t(e)%*%e
    out
  }


  #sse <- thr_sse(yt,qq,rr)
  thr_sse <- function(y,q,r){
    nq <- nrow(q)
    sse <- matrix(NA,nq,1)
    for (qi in 1:nq){
      if (r[1]==0) {
        rr <- q[qi]
      }else{
        rr <- rbind(r,q[qi])
      }

      rr <- as.matrix(sort(rr))
      nr = nrow(rr)
      #xx <- cbind(xt)
      dd <- matrix(NA,nrow=nrow(thresh),ncol=nr)
      xx <- NULL

      for (j in 1:nr){
        dd[,j] <- (thresh < rr[j])
        d <- dd[,j]
        if (j>1) d <- d - dd[,(j-1)]
        if (j == 1){
          xx <- cbind(xx,tr(xt*d,N,iT))
        }else{
          xx <- cbind(xx,tr(regx*d,N,iT))
        }
      }
      d <- 1-dd[,nr]
      xx <- cbind(xx,tr(regx*d,N,iT))
      sse[qi] <- sse_calc(y,xx)
    }
    sse
  }

  r_est <- function(y,r,trim_rate){
    if (max(r)==0){
      qq <- qq1;
      rr <- 0;
    }else{
      rr <- as.matrix(sort(r))
      i <- as.matrix(seq(1,qn1,by=1))
      nn <- colSums(qq1%*%matrix(1,1,nrow(rr))< matrix(1,nrow(qq1),1)%*%t(rr))
      nn <- as.matrix(nn)
      qnt <- grid_num*trim_rate
      ii1 <- (i%*%matrix(1,1,nrow(nn)))<(matrix(1,qn1,1)%*%t(nn+qnt))
      ii2 <- (i%*%matrix(1,1,nrow(nn)))<(matrix(1,qn1,1)%*%t(nn-qnt))
      ii <- (ii1-ii2)%*%matrix(1,nrow(rr),1)
      qq <- as.matrix(qq1[ii!=1])
    }

    sse <- thr_sse(y,qq,rr)
    rihat <- which.min(sse)
    list(sse_b=sse[rihat],rhat_b=qq[rihat])
  }

  #r = rbind(rhat21,rhat22)
  model <- function(r,trim_rate){
    if (max(r)==0){
      qq <- qq1;
      rr <- 0;
    }else{
      rr <- as.matrix(sort(r))
      i <- as.matrix(seq(1,qn1,by=1))
      nn <- colSums(qq1%*%matrix(1,1,nrow(rr))< matrix(1,nrow(qq1),1)%*%t(rr))
      nn <- as.matrix(nn)
      qnt <- grid_num*trim_rate
      ii1 <- (i%*%matrix(1,1,nrow(nn)))<(matrix(1,qn1,1)%*%t(nn+qnt))
      ii2 <- (i%*%matrix(1,1,nrow(nn)))<(matrix(1,qn1,1)%*%t(nn-qnt))
      ii <- (ii1-ii2)%*%matrix(1,nrow(rr),1)
      qq <- as.matrix(qq1[ii!=1])
    }

    sse <- thr_sse(yt,qq,rr)
    rihat <- which.min(sse)
    rhat <- qq[rihat]
    sse1 <- sse[rihat]

    if (max(r) != 0){
      rrr <- sort(rbind(rr,rhat))
    }else{ rrr <- rhat }

    rrr <- as.matrix(rrr)
    nr <- nrow(rrr)
    xx <- NULL
    dd <- matrix(NA,nrow=nrow(thresh),ncol=nr)
    for (j in 1:nr){
      dd[,j] <- (thresh < rrr[j])
      d <- dd[,j]
      if (j>1) d <- d - dd[,(j-1)]
      if (j == 1){
        xx <- cbind(xx,tr(xt*d,N,iT))
      }else{
        xx <- cbind(xx,tr(regx*d,N,iT))
      }
    }
    d <- 1-dd[,nr]
    xx <- cbind(xx,tr(regx*d,N,iT))
    xxi <- ginv(t(xx)%*%xx)
    beta <- xxi%*%(t(xx)%*%yt)
    e <- yt - xx%*%beta
    xxe <- xx*(e%*%matrix((1),nrow=1,ncol=ncol(xx)))
    xxe <- t(xxe)%*%xxe
    het.covmat <- xxi%*%xxe%*%xxi;
    homo.covmat <- xxi*as.vector((t(e)%*%e))/(ty-N-ncol(xx))
    #sehet <- as.matrix(sqrt(diag(xxi%*%xxe%*%xxi)))
    #sehomo <- as.matrix(sqrt(diag(xxi*as.vector((t(e)%*%e))/(ty-N-ncol(xx)))))
    #beta <- cbind(beta,sehomo,sehet)
    list(rhat=rhat,beta=beta,
         het.covmat = het.covmat,homo.covmat=homo.covmat)
  }


  if(im == 1){

    obj1 = model(r=0,trim_rate)
    beta = obj1$beta
    rhat = obj1$rhat
    het.covmat <- obj1$het.covmat
    homo.covmat <- obj1$homo.covmat

  }else if(im ==2){
    obj1 = model(r=0,trim_rate)
    rhat1 = obj1$rhat
    obj21 <- model(rhat1,trim_rate)
    rhat21 <- obj21$rhat
    obj22 <- model(rhat21,trim_rate)
    rhat22 <- obj22$rhat
    rhat = sort(c(rhat21,rhat22))
    beta = obj22$beta
    het.covmat <- obj22$het.covmat
    homo.covmat <- obj22$homo.covmat

  }else if(im == 3){
    obj1 = model(r=0,trim_rate)
    rhat1 = obj1$rhat
    #trim_2 = trimb[2]
    obj21 <- model(rhat1,trim_rate)
    rhat21 <- obj21$rhat
    obj22 <- model(rhat21,trim_rate)
    rhat22 <- obj22$rhat
    #trim_3 = trim[3]
    obj3 <- model(rbind(rhat21,rhat22),trim_rate)
    rhat31 <- obj3$rhat
    rhat = sort(c(rhat21,rhat22,rhat31))
    beta = obj3$beta
    het.covmat <- obj3$het.covmat
    homo.covmat <- obj3$homo.covmat



  }else if (thnum >=4){
    stop("The number of threshold should be correctly specified.")
  }


  ####

  coef_names <- c(paste0("R1-",var_base))
  for(kki in 1:im){
    coef_names <- c(coef_names,paste0("R",kki+1,"-",var_rx))
  }
  beta = as.vector(beta)
  beta.het.se <- sqrt(diag(het.covmat))
  beta.homo.se <- sqrt(diag(homo.covmat))
  names(beta.het.se) = names(beta.homo.se) = names(beta) = coef_names
  gamma <- rhat

  eval = list()
  eval$beta = beta;
  eval$gamma = gamma
  eval$beta.homo.se = beta.homo.se
  eval$beta.het.se = beta.het.se

  return(eval)


}

