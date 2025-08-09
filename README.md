# panelthreshold

`panelthreshold` is an `R` package to estimate the dynamic multiple threshold regression with threshold 
    effects and endogeneity. To handle endogeneity, we apply the rst-differenced GMM framework
    to estimate the unknown parameters for a given number of thresholds. When the number of thresholds 
    is unknown, a modified information criterion is employed to determine the model complexity.  Moreover, 
    we develop a sup-Wald test statistic to facilitate inference on threshold effects. This R package is based on the paper:

Wan C, Li Y. Dynamic panel models with multi-threshold effects and endogeneity. *Statistics*, 2025: 1-21.
    
# How to install
You can install the development version from GitHub
```r
devtools::install_github("ChuangWAN1994/panelthreshold")
```

# Case 1: Dynamic Panel Model with Two Thresholds and Exogenous Regressors
```r
N <- 100
T <- 10
y <- matrix(NA,N,T)
x1 <- x2 <- matrix(NA,N,T)
alpha <- rnorm(N,0,1)
rho1 <- runif(N,-0.9,0.9)
rho2 <- runif(N,-0.9,0.9)
for(t in 1:T){
  xx1 <- if(t>1) x1[,t-1] else 0
  x1[,t] <- rho1*xx1 + 0.25*alpha + rnorm(N,0,0.5)
  xx2 <- if(t>1) x2[,t-1] else 0
  x2[,t] <- rho2*xx2 + 0.5*alpha + rnorm(N,0,0.5)
  yy <- if(t>1) y[,t-1] else 0
  y[,t] <- x1[,t] + alpha  + (1-2*x1[,t])*(x2[,t]> -0.25) + (-1+ 2* x1[,t])*(x2[,t]>0.25) + 0.2*yy + rnorm(N,0,0.5)
}
y <- c(t(y))
x1 <- c(t(x1))
x2 <- c(t(x2))
data = cbind(y,x2,x1,1)
colnames(data) = c("y","x2","x1","cons")

dep = "y"
thre = "x2"
var_rx = c("cons","x1")
var_base = c("lag_y","x1")
var_endo = NULL
var_iv = NULL
use = NewMTDP(data,thre_type=c("jump"),dep,thre,var_base,var_rx, var_endo,var_iv,T)

im = 2
h_0 = c(0.5)
par = NULL
trim_rate=0.15;
grid_num=0;
avg_num=0
im = 2
h_0 = c(0.5)
MTDP(use,im,h_0,par,trim_rate,grid_num,avg_num,nboots = 0 )
```

# Case 2: Dynamic Panel Model with Single Threshold Where the Threshold Variable is the Lagged Dependent Variable
```r
N <- 200
T <- 10
y <- matrix(NA,N,T)
for(t in 1:T){
  yy <- if(t>1) y[,t-1] else 0
  y[,t] <- (-0.5*yy)*(yy <= 0) + (-1.8+0.7*yy)*(yy>0) + rnorm(N,0,0.5)
}
y <- c(t(y))
data = cbind(1,y)
colnames(data) = c("cons","y")

dep = "y"
thre = "lag_y"
var_rx = c("cons","lag_y")
var_base = c("lag_y")
var_endo = NULL
var_iv = NULL
use = NewMTDP(data,thre_type=c("jump"),dep,thre,var_base,var_rx, var_endo,var_iv,T)

im = 1
h_0 = c(0.5)
par = NULL 
trim_rate=0.15;
grid_num=0;
avg_num=0
MTDP(use,im,h_0,par,trim_rate,grid_num,avg_num,nboots = 0 )
```

# Case 3: Dynamic panel with two thresholds, where both threshold variables and explanatory variables exhibit endogeneity
```r
N <- 100
T <- 10
y <- matrix(NA,N,T)
x1 <- x2 <- matrix(NA,N,T)
alpha <- rnorm(N,0,1)
rho1 <- runif(N,-0.9,0.9)
rho2 <- runif(N,-0.9,0.9)
for(t in 1:T){
  err1 <- rnorm(N,0,0.5)
  err2 <- rnorm(N,0,0.5)
  err <- rnorm(N,0,1)
  xx1 <- if(t>1) x1[,t-1] else 0
  x1[,t] <- rho1*xx1 + 0.25*alpha + err1
  xx2 <- if(t>1) x2[,t-1] else 0 
  x2[,t] <- rho2*xx2 + 0.5*alpha + err2
  yy <- if(t>1) y[,t-1] else 0
  y[,t] <- x1[,t] + alpha  + (1-2*x1[,t])*(x2[,t]> -0.25) + 
    (-1+ 2*x1[,t])*(x2[,t]>0.25) + 0.2*yy + (0.25*err + 0.5*err1 + 0.5*err2)
}
y <- c(t(y))
x1 <- c(t(x1))
x2 <- c(t(x2))
data = cbind(y,x2,x1,1)
colnames(data) = c("y","x2","x1","cons")

dep = "y"
thre = "x2"
var_rx = c("cons","x1")
var_base = c("lag_y","x1")
var_endo = c("x1","x2") # set the endogeneous variable
var_iv = NULL
use = NewMTDP(data,thre_type=c("jump"),dep,thre,var_base,var_rx,var_endo,var_iv,T)

im = 2
h_0 = c(0.5)
par = NULL
trim_rate=0.15;
grid_num=0;
avg_num=0
MTDP(use,im,h_0,par,trim_rate,grid_num,avg_num,nboots = 0 )
```

# Case 4: Dynamic panel model with two thresholds, featuring endogenous threshold variables and explanatory variables, but with available exogenous instrumental variables.

```r


N <- 200
T <- 10
y <- matrix(NA,N,T)
x1 <- x2 <- matrix(NA,N,T)
iv1 <- iv2 <- matrix(NA,N,T)
alpha <- rnorm(N,0,1)
rho1 <- runif(N,-0.9,0.9)
rho2 <- runif(N,-0.9,0.9)
for(t in 1:T){
  rho <- 0.75
  Sig <- matrix(rho,3,3)
  diag(Sig) <- 1
  Err <- mvtnorm::rmvnorm(N,mean=rep(0,3),sigma = Sig) 
  ivv1 <- if(t>1) iv1[,t-1] else 0
  iv1[,t] <- rho1*ivv1 + 0.25*alpha + rnorm(N,0,0.5)
  x1[,t] <- iv1[,t] + Err[,1]
  ivv2 <- if(t>1) iv2[,t-1] else 0
  iv2[,t] <- rho2*ivv2 + 0.5*alpha + rnorm(N,0,0.5)
  x2[,t] <- iv2[,t] + Err[,2]
  yy <- if(t>1) y[,t-1] else 0
  y[,t] <- x1[,t] + alpha  + (-2*x1[,t])*(x2[,t]> -0.25) + ( 2* x1[,t])*(x2[,t]>0.25)  + Err[,3] + 0.2*yy
}
y <- c(t(y))
x1 <- c(t(x1))
x2 <- c(t(x2))
iv1 <- c(t(iv1))
iv2 <- c(t(iv2))
data = cbind(y,x2,x1,1,iv1,iv2)
colnames(data) = c("y","x2","x1","cons","iv1","iv2")

dep = "y"
thre = "x2"
var_rx = c("cons","x1")
var_base = c("lag_y","x1")
var_endo = c("x1","x2") # set the endogeneous variable
var_iv = c("iv1","iv2") # set the additional instrument variable
use = NewMTDP(data,thre_type=c("jump"),dep,thre,var_base,var_rx,
              var_endo,var_iv,T)
im = 2
h_0 = c(0.5)
par = NULL
trim_rate=0.15;
grid_num=0;
avg_num=0

MTDP(use,im,h_0,par,trim_rate,grid_num,avg_num,nboots = 0 )
```
