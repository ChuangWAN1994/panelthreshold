#

tr <- function(yy, N, TT){

  # Output is the demeaned data with the same dimension as input
  # NT * 1 or NT * p

  yy <- as.matrix(yy)

  if(dim(yy)[1] != N*TT) print("Error! Dimension of
                                 inputs in demean is wrong!")

  p <- dim(yy)[2];

  if( p == 1){
    y.temp <- matrix(yy, nrow = TT);
    m <- apply(y.temp, 2, mean);
    y.temp <- y.temp - matrix( rep(m, times = TT), nrow = TT,
                               ncol = N, byrow = TRUE);
    y <- matrix(y.temp, nrow = N*TT)
    return(y)
  }else{
    y <- matrix(0, N*TT, p);
    for(j in 1:p){
      y.temp <- matrix( yy[,j], nrow = TT);
      m <- apply(y.temp, 2, mean);
      y.temp <- y.temp - matrix( rep(m, times = TT), nrow = TT,
                                 ncol = N, byrow = TRUE);
      y[,j] <- matrix(y.temp, nrow = N*TT);
    }
    return(y)
  }
}



#//....................Some main function in \emph{NewMTDP} function.....

FD_con <- function(y_mat,T){
  y_mat_fd = as.matrix(y_mat[,2:T] - y_mat[,1:(T-1)])
}


moment_mat_con <- function(T,k_1,k_endo)
{
  if(k_endo == 0){
    moment_mat = matrix(1,T-1,k_1+1)
    moment_mat[,2] = 0:(T-2)
    #This matrix contains information about moment conditions.
    #(Here we assume there's no endogenous x_{i, t})
    #The first row is for \delta x_{i2} = x_{i2} - x_{i1}
    #The second row is for \delta x_{i3} = x_{i3} - x_{i2}
    #... and the last row is for \delta x_{iT} = x_{iT} - x_{iT-1}
    #The first column is for use of constant(1's).
    #   The second column is for use of lagged y_{i, t}'s.
    #(For different t's we use different lagged terms of y_{i, t})
    #The other columns are for use of differenced x_{i, t}'s

  }else{
    moment_mat = matrix(1,T-1,k_1+1)
    moment_mat[,2] = 0:(T-2)
    moment_mat[,(k_1+2-k_endo):(k_1+1)] = kronecker(matrix(1,1,k_endo),
                                                    matrix(c(0:(T-2))))
    #matrix(c(0:(T-2)),1,k_endo)

    #Additionally, when there are (# = k_endo) endogenous x_{i, t},
    #  The last columns are for use of endogenous lagged x_{i, t}'s
  }

  return(moment_mat)
}



Z_mat_con <- function(y_mat,x_mat,x_mat_fd,moment_mat,k_endo,flag_static){

  N = nrow(y_mat)
  T = ncol(y_mat)
  flag_z = 1 #//index indicator to fill out Z_mat
  k_1 = nrow(x_mat)/N

  if(flag_static == 1){
    ## For static model
    L = sum(moment_mat) - sum(moment_mat[,2])
    #/Here moment condition starts from t=2
    #However, the 2nd column is excluded since it is for lagged y_{i, t}
    #(The 1st row of moment_mat starts from t=2)


    Z_mat = matrix(0,L,N)
    for(kk1 in 1:(T-1)){
      m_temp = moment_mat[kk1,]
      # //Moment conditions for t=(kk1+1)

      for(ll1 in 1:(k_1+1)){
        if(m_temp[ll1] == 0){
          next() # no IV to stack
        }else if(ll1 == 1){
          Z_mat[flag_z,] = matrix(1,1,N)
          flag_z = flag_z+1 # Stack 1's(constant)
        }else if(ll1 ==2){
          next() #// Do not stack lagged y_{i, t}'s
        }else if(ll1 <= (k_1+1-k_endo)){
          Z_mat[flag_z, ] = t(x_mat_fd[((ll1-2)*N+1):((ll1-1)*N),kk1])
          flag_z = flag_z+1 #// Stack \delta x_{i, t}'s for exogeneous x_{i, t}
        }else{
          Z_mat[flag_z:(flag_z+m_temp[ll1]-1),] =
            t(x_mat[((ll1-2)*N+1):((ll1-1)*N),1:(m_temp[ll1])])
          flag_z = flag_z+m_temp[ll1] #// Stack x_{i, t}'s for endogenous x_{i, t}
        }
      }

    }

  } else{

    L = sum(moment_mat[2:nrow(moment_mat),])
    #//Here moment condition stars from t=3
    # (The 1st row of moment_mat starts from t=2)
    Z_mat = matrix(0,L,N)
    for(kk1 in 2:(T-1)){
      m_temp = moment_mat[kk1,]
      # Moment conditions for t=(kk1+1)
      for(ll1 in 1:(k_1+1)){
        if(m_temp[ll1] == 0){
          next()  # No IV to stack. But this will not happen normally.
        }else if(ll1 == 1){
          Z_mat[flag_z,] = matrix(1,1,N)
          flag_z = flag_z+1 #Stack 1's(constant)
        }else if(ll1 == 2){
          Z_mat[flag_z:(flag_z+m_temp[ll1]-1),] =
            t(y_mat[,(kk1-m_temp[ll1]):(kk1-1)]) # Stack lagged y_{i, t}'s
          flag_z = flag_z + m_temp[ll1]
        }else if(ll1 <= (k_1+1-k_endo)){
          Z_mat[flag_z,] = t(x_mat_fd[((ll1-2)*N+1):((ll1-1)*N),kk1])
          flag_z = flag_z+1 #Stack \delta x_{i, t}'s for exogeneous x_{i, t}

        }else{
          Z_mat[flag_z:(flag_z+m_temp[ll1]-1),] =
            t(x_mat[((ll1-2)*N+1):((ll1-1)*N),1:(m_temp[ll1])])
          flag_z = flag_z + m_temp[ll1] # Stack x_{i, t}'s for endogenous x_{i, t}

        }
      }
    }
  }

  return(Z_mat)
}




Z_mat_IV_add <- function(Z_mat,moment_mat,IV_mat,k_inst,flag_static)
{
  #//This is a function to extend Z_mat and moment_mat,
  #//so that they can respond to additional given IV's.
  #//Activated only if there are additional IV's (optional)

  L = nrow(Z_mat)
  N = ncol(Z_mat)
  T = nrow(moment_mat) + 1

  moment_add = matrix(1,T-1, k_inst)
  #//Indicates IV's to add

  flag_z1 = 1
  flag_z2 = 1
  #//Index for Z_mat and new Z_mat, respectively

  Z_copy = Z_mat #// copy of original Z_mat

  if (flag_static == 1){
    #//For static model(special case)

    Z_mat = matrix(0,L + sum(moment_add), N)
    #// This will be new Z_mat for static model

    for(kk4 in 2:T){
      m_i1 = sum(moment_mat[kk4-1,]) - moment_mat[kk4-1, 2]
      #// number of original IV's at t = kk4

      IV_add = matrix(IV_mat[, kk4], ncol=N,byrow=TRUE)
      #//(k_inst * N) matrix, this will be added

      IV_ori = Z_copy[flag_z1:(flag_z1+m_i1-1), ]
      #// original IV's at t = kk4

      flag_z1 = flag_z1 + m_i1 #// Index updating

      Z_mat[flag_z2:(flag_z2+m_i1-1), ] = IV_ori
      Z_mat[(flag_z2+m_i1):(flag_z2+m_i1+k_inst-1), ] = IV_add

      flag_z2 = flag_z2 + m_i1 + k_inst
    }

  }else{
    #//For dynamic model


    Z_mat = matrix(0,L + sum(moment_add) - k_inst, N)
    #// This will be new Z_mat for dynamic model

    for(kk4 in 3:T){
      m_i1 = sum(moment_mat[kk4-1, ])
      #// number of original IV's at t = kk4

      IV_add = matrix(IV_mat[, kk4], ncol=N,byrow=TRUE)
      #//(k_inst * N) matrix, this will be added

      IV_ori = Z_copy[flag_z1:(flag_z1+m_i1-1), ]
      #// original IV's at t = kk4

      flag_z1 = flag_z1 + m_i1 #// Index updating

      Z_mat[flag_z2:(flag_z2+m_i1-1), ] = IV_ori
      Z_mat[(flag_z2+m_i1):(flag_z2+m_i1+k_inst-1),] = IV_add

      flag_z2 = flag_z2 + m_i1 + k_inst


    }
  }
  moment_mat = cbind(moment_mat, moment_add)
  #//Add columns for additional IV

  return(list(moment_mat=moment_mat,Z_mat = Z_mat))
}




GMM_var_con <- function(y_mat_fd, x_mat_fd,
                        rx_mat,
                        q_mat,
                        moment_mat,
                        N,T,flag_static)
{
  if(flag_static == 1){
    L = sum(moment_mat) - sum(moment_mat[,2])
    # In case of static model
  }else{
    L = sum(moment_mat[2:(nrow(moment_mat)),])
    #/In case of dynamic model(default)
    # The first row was excluded since it is for t=2
    #//*For t=2 (LHS) is (y_{i, 2} - y_{i, 1}), so there's no proper (RHS)

  }

  #N = nrow(y_mat_fd) #the number of individuals
  k_x = nrow(x_mat_fd)/N #the number of regressors including lagged y_{i,t-1}
  k_rx = nrow(rx_mat)/N

  y_mat_rep = matrix(0,L,N)
  x_temp1 = matrix(0,L,N*k_x)
  rx_temp2 = rx_temp3 = matrix(0,L,N*k_rx)
  q_temp1 = matrix(0,L,N)
  q_temp2 = matrix(0,L,N)
  flag_r1 = 1
  if(flag_static == 1){
    ## For static model
    for(kk in 1:(T-1)){
      m_temp = sum(moment_mat[kk,]) - sum(moment_mat[kk,2])
      # of stacked moment condition for t = kk
      if(m_temp == 0){
        next()
      }else{
        y_mat_rep[flag_r1:(flag_r1+m_temp-1),] = t(kronecker(y_mat_fd[,kk], matrix(1,1,m_temp)))
        ##/Using # (kronecker product), repeat the desired column
        x_temp1[flag_r1:(flag_r1+m_temp-1),] = t(kronecker(x_mat_fd[,kk],matrix(1,1,m_temp)))

        rx_temp2[flag_r1:(flag_r1+m_temp-1),] = t(kronecker(rx_mat[,kk+1],matrix(1,1,m_temp)))
        rx_temp3[flag_r1:(flag_r1+m_temp-1),] = t(kronecker(rx_mat[,kk],matrix(1,1,m_temp)))

        q_temp1[flag_r1:(flag_r1+m_temp-1),] = t(kronecker(q_mat[,kk+1],matrix(1,1,m_temp)))
        q_temp2[flag_r1:(flag_r1+m_temp-1),] = t(kronecker(q_mat[,kk],matrix(1,1,m_temp)))
        flag_r1 = flag_r1+m_temp

      }
    }
  }else{
    ## For dynamic model (Default)
    for(kk in 2:(T-1)){
      m_temp = sum(moment_mat[kk,]) #// # of stacked moment condition for t = kk
      if(m_temp == 0){
        next()
      }else{
        y_mat_rep[flag_r1:(flag_r1+m_temp-1),]=t(kronecker(y_mat_fd[,kk], matrix(1,1,m_temp)))
        #//Using # (kronecker product), repeat the desired column

        x_temp1[flag_r1:(flag_r1+m_temp-1),] = t(kronecker(x_mat_fd[,kk],matrix(1,1,m_temp)))

        rx_temp2[flag_r1:(flag_r1+m_temp-1),] = t(kronecker(rx_mat[,kk+1],matrix(1,1,m_temp)))
        rx_temp3[flag_r1:(flag_r1+m_temp-1),] = t(kronecker(rx_mat[,kk],matrix(1,1,m_temp)))

        q_temp1[flag_r1:(flag_r1+m_temp-1),] = t(kronecker(q_mat[,kk+1],matrix(1,1,m_temp)))
        q_temp2[flag_r1:(flag_r1+m_temp-1),] = t(kronecker(q_mat[,kk],matrix(1,1,m_temp)))

        flag_r1 = flag_r1+m_temp
      }
    }
  }

  return(list(y_mat_rep=y_mat_rep,
              x_temp1 = x_temp1,
              rx_temp2 = rx_temp2,
              rx_temp3 = rx_temp3,
              q_temp1 = q_temp1,
              q_temp2 = q_temp2))
}


#//..............Some main function in \emph{MTDP} function.....


GMM_W_n_static <- function(Z_mat,moment_mat,T)
{
  #//This is a function to construct weight matrix for 1st step
  #//This is for static model (optional)

  L = nrow(Z_mat)
  N = ncol(Z_mat)

  W_n_temp1 = matrix(0,L, L)
  W_n_temp2 = matrix(0,L, L)

  flag_r2 = 1  #//To indicate the next column to fill up
  flag_z1 = 1
  flag_z2 = 1  # //To indicate the index of Z_mat to copy

  moment_vec = rowSums(moment_mat) - moment_mat[, 2]
  #//for static model

  for(kk2 in 2:T){
    m_num = moment_vec[kk2-1]
    #// number of IV's at period t (where LHS is y_{i, t} - y_{i, t-1})
    if(m_num == 0){
      next() #/ When there's no IV. But this will not happen normally.
    }else{
      if(flag_r2 == 1){
        Z_temp1 = matrix(Z_mat[1:m_num,],nrow=m_num)
        #  (m_1 * N) matrix, where m_1 is # of moment condition at t_0

        W_n_temp2[1:m_num,1:m_num] = 1/N*(Z_temp1 %*% t(Z_temp1))

        flag_r2 = flag_r2 + m_num
        flag_z2 = flag_z2 + m_num
      }else{
        m_num_l = moment_vec[kk2-2] # // lagged m_num
        Z_temp1 = Z_mat[flag_z1:(flag_z1+m_num_l-1), ]
        Z_temp2 = Z_mat[flag_z2:(flag_z2+m_num-1), ]
        #//Both are building blocks of (10) in Seo & Shin(2016).
        flag_z1 = flag_z1 + m_num_l
        flag_z2 = flag_z2 + m_num #// updating

        W_n_temp1[(flag_r2-(m_num_l)):(flag_r2-1), flag_r2:(flag_r2+m_num-1)] =
          -(1/N)*(Z_temp1 %*% t(Z_temp2))

        W_n_temp2[flag_r2:(flag_r2+m_num-1), flag_r2:(flag_r2+m_num-1)] =
          (1/N)*(Z_temp2 %*% t(Z_temp2))

        flag_r2 = flag_r2 + m_num

      }
    }
  }
  W_n <- try(solve(W_n_temp1 + t(W_n_temp1) + 2*W_n_temp2),silent=TRUE)
  if(is(W_n,"try-error")){
    W_n = ginv(W_n_temp1 + t(W_n_temp1) + 2*W_n_temp2)
  }


  return(W_n)
}



GMM_W_n_con <- function(Z_mat,moment_mat,T)
{
  #//This is a function to construct weight matrix for 1st step
  #//This is for dynamic model (default)
  L = nrow(Z_mat)
  N = ncol(Z_mat)

  W_n_temp1 = matrix(0,L,L)
  W_n_temp2 = matrix(0,L,L)

  flag_r2 = 1 #//To indicate the next column to fill up
  flag_z1 = 1
  flag_z2 = 1 #//To indicate the index of Z_mat to copy

  moment_vec = rowSums(moment_mat)

  #//for dynamic model(default)
  for(kk2 in 3:T){
    m_num = moment_vec[kk2-1] # the number of IV's at period t (where LHS is y_{i, t} - y_{i, t-1})

    if(m_num == 0){
      #continue
      #// When there's no IV. But this will not happen normally.
      next()
    }else{
      if(flag_r2 == 1){
        Z_temp1 = matrix(Z_mat[1:(m_num),],nrow=m_num)
        #(m_1*N) matrix where m_1 is the number of moment condition at  t_0

        W_n_temp2[1:m_num,1:m_num] = 1/N*(Z_temp1 %*% t(Z_temp1))
        flag_r2 = flag_r2 + m_num
        flag_z2 = flag_z2 + m_num
      }else{
        m_num_l = moment_vec[kk2-2] #lagged m_num
        Z_temp1 = Z_mat[flag_z1:(flag_z1+m_num_l-1),]
        Z_temp2 = Z_mat[flag_z2:(flag_z2+m_num-1),]
        #//Both are building blocks of (10) in Seo & Shin(2016).

        flag_z1 = flag_z1 + m_num_l
        flag_z2 = flag_z2 + m_num # updating

        W_n_temp1[(flag_r2-m_num_l):(flag_r2-1),flag_r2:(flag_r2+m_num-1)]=
          -1/N*(Z_temp1 %*% t(Z_temp2))
        # W_n_temp2[flag_r2..flag_r2+m_num-1, flag_r2..flag_r2+m_num-1] ///
        #= (1/N)*(Z_temp2 * Z_temp2')
        W_n_temp2[flag_r2:(flag_r2+m_num-1),flag_r2:(flag_r2+m_num-1)] = 1/N*(Z_temp2 %*% t(Z_temp2))
        flag_r2 = flag_r2 + m_num
      }
    }
  }
  W_n = try(solve(W_n_temp1+t(W_n_temp1)+2*W_n_temp2), silent = TRUE)
  if(is(W_n,"try-error")){
    W_n = ginv(W_n_temp1+t(W_n_temp1)+2*W_n_temp2)
  }

  return(W_n)
}



GMM_W_n_2_con <- function(ep_collect,Z_mat)
{
  #//This is a function to construct weight matrix for 2nd step
  L = nrow(Z_mat)
  N = ncol(Z_mat)

  g_temp_1 = matrix(0,L, L)
  g_temp_2 = matrix(0,L, 1)
  #/Matrices that are corresponding to the 1st, 2nd term of (11) in the paper.

  for(ii3 in 1:N){
    g_i_hat = (ep_collect[,ii3])*(Z_mat[,ii3])
    g_temp_1 = g_temp_1 + g_i_hat %*% t(g_i_hat)
    g_temp_2 = g_temp_2 + g_i_hat
  }
  W_n_2 = try(solve(1/N*g_temp_1 - (1/(N^2))*(g_temp_2 %*% t(g_temp_2))),silent = TRUE)
  if(is(W_n_2,"try-error")){
    W_n_2 = ginv(1/N*g_temp_1 - (1/(N^2))*(g_temp_2 %*% t(g_temp_2)))
  }

  return(W_n_2)

}


ep_cal_Multithreshold <- function(r_hat, x_temp1, rx_temp2, rx_temp3,q_temp1,q_temp2,y_mat_rep,GMM_est)
{
  #// This is a function to stack estimated errors from the 1st step

  L = nrow(q_temp1)
  N = ncol(q_temp1)
  k_crx = ncol(x_temp1)/N
  k_rx = ncol(rx_temp2)/N

  thre_num = length(r_hat)

  ep_collect = matrix(0,L,N)
  #(L*N) matrix that will contain all collected residuals from the first GMM

  ind_temp2 = ind_temp3 = list()
  for(kki in 1:thre_num){
    r_th = r_hat[kki]
    ind_temp2[[kki]] = (q_temp1 > r_th)
    ind_temp3[[kki]] = (q_temp2 > r_th)
  }

  #ind_temp2 = (q_temp1 > r_hat)
  #ind_temp3 = (q_temp2 > r_hat)

  ind_vec1 = seq(0,(k_crx-1)*N, length= k_crx)
  ind_vec2 = seq(0,(k_rx-1)*N,length=k_rx)
  #// To pick general k_1 variables

  for(ii2 in 1:L){
    for(ii3 in 1:N){
      ii3_rep1 = matrix(ii3,k_crx,1)
      ii3_rep2 = matrix(ii3,k_rx,1)
      xreg = c(x_temp1[ii2,(ind_vec1+ii3_rep1)]) #(x_temp1[ii2,ind_vec+ii3_rep])
      xreg_thre = NULL
      for(kki in 1:thre_num){
        xreg_thre = c(xreg_thre,
                      (ind_temp2[[kki]][ii2,ii3]*rx_temp2[ii2,ind_vec2+ii3_rep2]-ind_temp3[[kki]][ii2,ii3]*
                         rx_temp3[ii2,ind_vec2+ii3_rep2]))
      }
      ep_collect[ii2,ii3] = y_mat_rep[ii2,ii3] - t(GMM_est) %*%
        matrix(c(xreg,xreg_thre),ncol=1)
    }
  }

  return(ep_collect)


}

ep_cal_multikink <- function(r_hat,x_temp1,rx_temp2,rx_temp3,
                             q_temp1,q_temp2,y_mat_rep,GMM_est)
{
  #// This is a function to stack estimated errors from the 1st step
  #// Especially, for the kink model case - i.e. only the variable
  #// (q_{i, t} - r) * 1_{q_{i, t} > r} changes according to threshold
  L = nrow(q_temp1)
  N = ncol(q_temp1)
  k_crx = ncol(x_temp1)/N

  kink_num = length(r_hat)
  ep_collect = matrix(0,L,N)
  #//(L*N) matrix that will contain all collected residuals from the first GMM

  ind_temp2 = ind_temp3 = list()
  for(kki in 1:kink_num){
    r_th = r_hat[kki]
    ind_temp2[[kki]] = (q_temp1 > r_th)
    ind_temp3[[kki]] = (q_temp2 > r_th)
  }

  ind_vec = seq(0,(k_crx-1)*N,length=k_crx)
  # // To pick general k_1 variables

  for(ii2 in 1:L){
    for(ii3 in 1:N){
      ii3_rep = matrix(ii3,k_crx,1)
      xreg = c(x_temp1[ii2,(ind_vec+ii3_rep)]) #(x_temp1[ii2,ind_vec+ii3_rep])
      xreg_kink = NULL
      for(kki in 1:kink_num){
        xreg_kink = c(xreg_kink,(ind_temp2[[kki]][ii2,ii3] * (q_temp1[ii2,ii3]-r_hat[kki]) -
                                   ind_temp3[[kki]][ii2,ii3] * (q_temp2[ii2,ii3]-r_hat[kki])))
      }
      ep_collect[ii2,ii3] = y_mat_rep[ii2,ii3] - t(GMM_est) %*%
        matrix(c(xreg,xreg_kink),ncol=1)
    }
  }

  return(ep_collect)

}



#//.......................Some function for the Nelder-Meader algorithm


consRegion <- function(im,grid_lower,grid_upper)
{

  ui1 <- diag(1,im); ui2 <- diag(-1,im)
  ci1 <- rep(grid_lower,im); ci2 <- rep(-grid_upper,im)
  if(im==1){
    ui <- rbind(ui1,ui2); ci <- c(ci1,ci2)
  }else if(im >= 2){
    ui3 <- matrix(0,im-1,im)
    for(l in 1:(im-1)){
      ui3[l,(l:(l+1))] <- c(-1,1)
    }
    ci3 <- rep(0,im-1)
    ui <- rbind(ui1,ui2,ui3)
    ci <- c(ci1,ci2,ci3)
  }
  return(list(ui=ui,ci=ci))

}



###jump threshold function
MMoptim_jump <- function(ui,ci,par,y_mat_rep,g_1n_bar,Z_mat,
                         x_temp1,rx_temp2,rx_temp3,
                         q_temp1,q_temp2,W_n)
{
  nr <- ncol(par)
  N <- ncol(y_mat_rep)
  J_vec <- c()
  par_est_mat <- c()
  GMM_est_mat <- c()
  ep_collect_mat <- c()

  for(rr in 1:nr){
    par_cur <- par[,rr]
    mio <- try(constrOptim(par_cur,f=function(r) jump_objective_function(g_1n_bar,r,Z_mat,x_temp1,rx_temp2,rx_temp3,
                                                                         q_temp1,q_temp2,W_n),
                           grad=NULL, ui,
                           ci,method = "Nelder-Mead"),silent=TRUE)
    if(is(mio,"try-error")) next()
    par_est_cur <- as.matrix(mio$par)
    obj <- jump_GMM_est_Fun(y_mat_rep,g_1n_bar,par_est_cur,Z_mat,x_temp1,rx_temp2,rx_temp3,
                            q_temp1,q_temp2,W_n)
    ep_collect_cur = obj$ep_collect;
    W_n_2 = GMM_W_n_2_con(ep_collect_cur, Z_mat)

    mio2 <- try(constrOptim(par_cur,f=function(r) jump_objective_function(g_1n_bar,r,Z_mat,x_temp1,rx_temp2,rx_temp3,
                                                                          q_temp1,q_temp2,W_n_2),
                            grad=NULL, ui,
                            ci,method = "Nelder-Mead"),silent=TRUE)
    if(is(mio2,"try-error")) next()
    par_est_cur_2 <- as.matrix(mio2$par)

    obj2 <- jump_GMM_est_Fun(y_mat_rep,g_1n_bar,par_est_cur_2,Z_mat,x_temp1,rx_temp2,rx_temp3,
                             q_temp1,q_temp2,W_n_2)

    J_vec <- c(J_vec,obj2$J_result)
    par_est_mat <- cbind(par_est_mat,par_est_cur_2)
    GMM_est_mat <- cbind(GMM_est_mat,obj2$GMM_est)
    ep_collect_mat <- cbind(ep_collect_mat,obj2$ep_collect)
  }
  idx <- which.min(J_vec)
  J_result <- J_vec[idx]
  par_est <- par_est_mat[,idx]
  GMM_est <- GMM_est_mat[,idx]
  ep_collect <- ep_collect_mat[,((idx-1)*N+1):(idx*N)]

  return(list(par_est=par_est,GMM_est=GMM_est,
              ep_collect=ep_collect,J_result=J_result))
}


jump_objective_function <- function(g_1n_bar,r,Z_mat,x_temp1,rx_temp2,rx_temp3,
                                    q_temp1,q_temp2,W_n){

  #/This is a function to calculate GMM estimator for given threshold r_th
  L = nrow(q_temp1)
  N = ncol(q_temp1)
  k_crx = ncol(x_temp1)/N
  k_rx = ncol(rx_temp2)/N
  r <- sort(drop(r))
  thre_num <- length(r)
  g_2n_bar = matrix(0,L,(k_crx+k_rx*thre_num))

  for(kk3 in 1:k_crx){
    g_2n_bar[,kk3] = rowSums(Z_mat*x_temp1[,((kk3-1)*N+1):(kk3*N)])/N
  }
  for(kki in 1:thre_num){
    r_th <- r[kki]
    ind_temp2 = (q_temp1 > r_th)
    ind_temp3 = (q_temp2 > r_th)
    for(kk3 in 1:k_rx){
      g_2n_bar[,k_crx+(kki-1)*k_rx+kk3] = rowSums(Z_mat * (rx_temp2[,((kk3-1)*N+1):(kk3*N)]*ind_temp2
                                                           -rx_temp3[,((kk3-1)*N+1):(kk3*N)]*ind_temp3))/N
    }
  }
  GMM_est = try(solve(t(g_2n_bar)%*% W_n %*% g_2n_bar) %*% (t(g_2n_bar)%*%W_n%*%g_1n_bar), silent = TRUE )
  if(is(GMM_est,"try-error")){
    GMM_est = ginv(t(g_2n_bar)%*% W_n %*% g_2n_bar) %*% (t(g_2n_bar)%*%W_n%*%g_1n_bar)
  }
  g_n_bar = g_1n_bar - g_2n_bar %*% GMM_est
  J_result <- t(g_n_bar) %*% W_n %*% g_n_bar

  return(J_result)
}


jump_GMM_est_Fun <- function(y_mat_rep,g_1n_bar,r,Z_mat,x_temp1,rx_temp2,rx_temp3,
                             q_temp1,q_temp2,W_n,flag = 1)
{
  L = nrow(q_temp1)
  N = ncol(q_temp1)
  k_crx = ncol(x_temp1)/N
  k_rx = ncol(rx_temp2)/N
  r <- sort(drop(r))
  thre_num <- length(r)
  g_2n_bar = matrix(0,L,(k_crx+k_rx*thre_num))

  for(kk3 in 1:k_crx){
    g_2n_bar[,kk3] = rowSums(Z_mat*x_temp1[,((kk3-1)*N+1):(kk3*N)])/N
  }
  for(kki in 1:thre_num){
    r_th <- r[kki]
    ind_temp2 = (q_temp1 > r_th)
    ind_temp3 = (q_temp2 > r_th)
    for(kk3 in 1:k_rx){
      g_2n_bar[,k_crx+(kki-1)*k_rx+kk3] = rowSums(Z_mat * (rx_temp2[,((kk3-1)*N+1):(kk3*N)]*ind_temp2
                                                           -rx_temp3[,((kk3-1)*N+1):(kk3*N)]*ind_temp3))/N
    }
  }
  GMM_est = try(solve(t(g_2n_bar)%*% W_n %*% g_2n_bar) %*% (t(g_2n_bar)%*%W_n%*%g_1n_bar), silent = TRUE )
  if(is(GMM_est,"try-error")){
    GMM_est = ginv(t(g_2n_bar)%*% W_n %*% g_2n_bar) %*% (t(g_2n_bar)%*%W_n%*%g_1n_bar)
  }
  if(flag == 1){
    ep_collect = ep_cal_Multithreshold(r, x_temp1, rx_temp2, rx_temp3, q_temp1, q_temp2,
                                       y_mat_rep, GMM_est)
    g_n_bar = g_1n_bar - g_2n_bar %*% GMM_est
    J_result <- t(g_n_bar) %*% W_n %*% g_n_bar
    return(list(ep_collect=ep_collect,GMM_est=GMM_est,
                J_result=J_result))
  }else{
    return(list(GMM_est=GMM_est))
  }


}


###kink threshold function


MMoptim_kink <- function(ui,ci,par,y_mat_rep,g_1n_bar,Z_mat,
                         x_temp1,rx_temp2,rx_temp3,
                         q_temp1,q_temp2,W_n)
{
  nr <- ncol(par)
  N <- ncol(y_mat_rep)
  J_vec <- c()
  par_est_mat <- c()
  GMM_est_mat <- c()
  ep_collect_mat <- c()

  for(rr in 1:nr){
    par_cur <- par[,rr]
    mio <- try(constrOptim(par_cur,f=function(r) kink_objective_function(g_1n_bar,r,Z_mat,x_temp1,rx_temp2,rx_temp3,
                                                                         q_temp1,q_temp2,W_n),
                           grad=NULL, ui,
                           ci,method = "Nelder-Mead"),silent=TRUE)
    if(is(mio,"try-error")) next()
    par_est_cur <- as.matrix(mio$par)
    obj <- kink_GMM_est_Fun(y_mat_rep,g_1n_bar,par_est_cur,Z_mat,x_temp1,rx_temp2,rx_temp3,
                            q_temp1,q_temp2,W_n)
    ep_collect_cur = obj$ep_collect;
    W_n_2 = GMM_W_n_2_con(ep_collect_cur, Z_mat)

    mio2 <- try(constrOptim(par_cur,f=function(r) kink_objective_function(g_1n_bar,r,Z_mat,x_temp1,rx_temp2,rx_temp3,
                                                                          q_temp1,q_temp2,W_n_2),
                            grad=NULL, ui,
                            ci,method = "Nelder-Mead"),silent=TRUE)
    if(is(mio2,"try-error")) next()
    par_est_cur_2 <- as.matrix(mio2$par)

    obj2 <- kink_GMM_est_Fun(y_mat_rep,g_1n_bar,par_est_cur_2,Z_mat,x_temp1,rx_temp2,rx_temp3,
                             q_temp1,q_temp2,W_n_2)

    J_vec <- c(J_vec,obj2$J_result)
    par_est_mat <- cbind(par_est_mat,par_est_cur_2)
    GMM_est_mat <- cbind(GMM_est_mat,obj2$GMM_est)
    ep_collect_mat <- cbind(ep_collect_mat,obj2$ep_collect)
  }
  idx <- which.min(J_vec)
  J_result <- J_vec[idx]
  par_est <- par_est_mat[,idx]
  GMM_est <- GMM_est_mat[,idx]
  ep_collect <- ep_collect_mat[,((idx-1)*N+1):(idx*N)]

  return(list(par_est=par_est,GMM_est=GMM_est,
              ep_collect=ep_collect,J_result=J_result))
}


kink_objective_function <- function(g_1n_bar,r,Z_mat,x_temp1,rx_temp2,rx_temp3,
                                    q_temp1,q_temp2,W_n){

  #/This is a function to calculate GMM estimator for given threshold r_th
  L = nrow(q_temp1)
  N = ncol(q_temp1)
  k_crx = ncol(x_temp1)/N
  r <- sort(drop(r))
  thre_num <- length(r)
  g_2n_bar = matrix(0,L,k_crx+thre_num)

  for(kk3 in 1:k_crx){
    g_2n_bar[,kk3] = rowSums(Z_mat*x_temp1[,((kk3-1)*N+1):(kk3*N)])/N
  }

  for(kki in 1:thre_num){
    r_th <- r[kki]
    ind_temp2 = (q_temp1 > r_th)
    ind_temp3 = (q_temp2 > r_th)
    kink_var1 = ((q_temp1 - r_th)*ind_temp2 - (q_temp2 - r_th)*ind_temp3)
    g_2n_bar[, k_crx + kki] = rowSums(Z_mat*kink_var1)/N
  }
  GMM_est = try(solve(t(g_2n_bar)%*% W_n %*% g_2n_bar) %*% (t(g_2n_bar)%*%W_n%*%g_1n_bar), silent = TRUE )
  if(is(GMM_est,"try-error")){
    GMM_est = ginv(t(g_2n_bar)%*% W_n %*% g_2n_bar) %*% (t(g_2n_bar)%*%W_n%*%g_1n_bar)
  }
  g_n_bar = g_1n_bar - g_2n_bar %*% GMM_est
  J_result <- t(g_n_bar) %*% W_n %*% g_n_bar

  return(J_result)
}


kink_GMM_est_Fun <- function(y_mat_rep,g_1n_bar,r,Z_mat,x_temp1,rx_temp2,rx_temp3,
                             q_temp1,q_temp2,W_n,flag = 1)
{
  L = nrow(q_temp1)
  N = ncol(q_temp1)
  k_crx = ncol(x_temp1)/N
  r <- sort(drop(r))
  thre_num <- length(r)
  g_2n_bar = matrix(0,L,k_crx+thre_num)

  for(kk3 in 1:k_crx){
    g_2n_bar[,kk3] = rowSums(Z_mat*x_temp1[,((kk3-1)*N+1):(kk3*N)])/N
  }

  for(kki in 1:thre_num){
    r_th <- r[kki]
    ind_temp2 = (q_temp1 > r_th)
    ind_temp3 = (q_temp2 > r_th)
    kink_var1 = ((q_temp1 - r_th)*ind_temp2 - (q_temp2 - r_th)*ind_temp3)
    g_2n_bar[, k_crx + kki] = rowSums(Z_mat*kink_var1)/N
  }
  GMM_est = try(solve(t(g_2n_bar)%*% W_n %*% g_2n_bar) %*% (t(g_2n_bar)%*%W_n%*%g_1n_bar), silent = TRUE )
  if(is(GMM_est,"try-error")){
    GMM_est = ginv(t(g_2n_bar)%*% W_n %*% g_2n_bar) %*% (t(g_2n_bar)%*%W_n%*%g_1n_bar)
  }
  g_n_bar = g_1n_bar - g_2n_bar %*% GMM_est

  if(flag == 1){
    ep_collect = ep_cal_multikink(r, x_temp1, rx_temp2, rx_temp3, q_temp1, q_temp2,
                                  y_mat_rep, GMM_est)
    g_n_bar = g_1n_bar - g_2n_bar %*% GMM_est
    J_result <- t(g_n_bar) %*% W_n %*% g_n_bar
    return(list(ep_collect=ep_collect,GMM_est=GMM_est,
                J_result=J_result))
  }else{
    return(list(GMM_est=GMM_est))
  }


}




#//.......................Some function for the grid-search algorithm


set_grid <- function(im,grid_num,grid_th,trim_rate)
{
  #qnt1 <- floor(grid_num * trim_rate)+1
  qnt1 <- floor(trim_rate/(1-2*trim_rate)*grid_num)+1
  RI <- combn(grid_num,im)
  if(im > 1) RI <- RI[,apply(RI,2,function(x) any(abs(diff(x)) >= qnt1))]
  Gmat <- apply(RI,c(1,2),function(x) grid_th[x])
  return(Gmat)
}


##//.......................variance and covariance matrix estimation



cov_mat_con_Multithreshold <- function(GMM_est_2,r_hat_2,q_mat,Z_mat,x_temp1,rx_temp2,
                                       rx_temp3, q_temp1,q_temp2,h_0,omega_hat)
{
  #//This function estimates covariance matrix. First, we construct ingredients.

  #//This function estimates covariance matrix. First, we construct ingredients.

  L = nrow(Z_mat)
  N = ncol(Z_mat)
  h_band = h_0 * 1.06 * N^(-0.2) * sd((matrix(q_mat, ncol=1,byrow=TRUE)))
  k_crx = ncol(x_temp1)/N
  k_rx <- ncol(rx_temp2)/N
  thre_num = length(r_hat_2)

  G_b = matrix(0,L,k_crx)
  G_d = matrix(0,L,k_rx*thre_num)
  G_r = matrix(0,L,thre_num)
  #//Three matrices are \hat{G_b}, \hat{G_d}, \hat{G_r} on Seo & Shin(2016), 173p.

  g_2n_bar = matrix(0,L,k_rx*thre_num+k_crx) #/This is identical to g_2n_bar in Seo & Shin(2016)


  for(ii6 in 1:k_crx){
    temp_crx = x_temp1[,((ii6-1)*N+1):(ii6*N)] #Extract (L*N) for one variable
    g_2n_bar[,ii6] = rowSums(Z_mat*temp_crx)/N
    G_b[,ii6] = -rowSums(Z_mat*x_temp1[,((ii6-1)*N+1):(ii6*N)])/N
  }

  for(kki in 1:thre_num){
    ind_temp2 = (q_temp1 > r_hat_2[kki])
    ind_temp3 = (q_temp2 > r_hat_2[kki])

    for(ii6 in 1:k_rx){
      g_2n_bar[,(k_crx+(kki-1)*(k_rx)+ii6)] = rowSums(Z_mat*(rx_temp2[,((ii6-1)*N+1):(ii6*N)]*ind_temp2-
                                                               rx_temp3[,((ii6-1)*N+1):(ii6*N)]*ind_temp3))/N ##
    }
  }



  for(kki in 1:thre_num){
    G_d[,((kki-1)*k_rx+1):(kki*k_rx)] = -g_2n_bar[,(k_crx+(kki-1)*k_rx+1):(kki*k_rx+k_crx)]
  }


  for(kki in 1:thre_num){
    for(ii7 in 1:N){
      ind_vec = N*(0:(k_rx-1))+ii7

      rx_temp_rc1 = cbind(rx_temp3[,ind_vec]) # no intercept, so matrix(1,L,1) should be removed
      rx_temp_rc2 = kronecker(dnorm((r_hat_2[kki]-q_temp2[,ii7])/h_band),matrix(1,1,k_rx))
      rx_temp_rc3 = cbind(rx_temp2[,ind_vec])
      rx_temp_rc4 = kronecker(dnorm((r_hat_2[kki]-q_temp1[,ii7])/h_band),matrix(1,1,k_rx))

      temp_z = (rx_temp_rc1*rx_temp_rc2 - rx_temp_rc3*rx_temp_rc4) %*%
        GMM_est_2[(k_crx+(kki-1)*k_rx+1):(k_crx+kki*k_rx)]
      G_r[,kki] = G_r[,kki] + Z_mat[, ii7]*temp_z
    }
    G_r[,kki] = G_r[,kki]/(N*h_band)
  }


  G = cbind(G_b, G_d, G_r)
  cov_mat = try(solve(t(G)%*%omega_hat%*%G)/N,silent = TRUE)
  if(is(cov_mat,"try-error")){
    cov_mat = ginv(t(G)%*%omega_hat%*%G)/N
  }

  return(cov_mat)
}

cov_mat_con_Multikink <- function(GMM_est_2,r_hat_2,q_mat,Z_mat,x_temp1,rx_temp2,
                                  rx_temp3,q_temp1,q_temp2,omega_hat)
{
  #//This function estimates covariance matrix. First, we construct ingredients.
  L = nrow(Z_mat)
  N = ncol(Z_mat)
  kink_num = length(r_hat_2)
  #h_band = h_0 * 1.06 * N^(-0.2) * sd((matrix(q_mat, ncol=1,byrow=TRUE)))
  k_1 = ncol(x_temp1)/N

  G_b = matrix(0,L, k_1)
  G_d = matrix(0,L, kink_num)
  G_r = matrix(0,L, kink_num)
  #//Three matrices are \hat{G_b}, \hat{G_d}, \hat{G_r} on Seo & Shin(2016), 173p.

  for(kki in 1:kink_num){
    ind_temp2 = (q_temp1 > r_hat_2[kki])
    ind_temp3 = (q_temp2 > r_hat_2[kki])
    #// This two (L*N) 0-1 matrices contain information whether q_{it} or q_{it-1} is
    #// larger than the threshold r_th
    G_d[, kki] = rowSums(Z_mat*((q_temp1-r_hat_2[kki])*ind_temp2 -
                                  (q_temp2-r_hat_2[kki])*ind_temp3))/N

    d_0 = GMM_est_2[k_1 + kki]
    G_r[,kki] = -d_0*rowSums(Z_mat*(ind_temp2 - ind_temp3))/N

  }

  #// Adjusted for kink model
  for(ii6 in 1:k_1){
    temp_x = x_temp1[,((ii6-1)*N+1):(ii6*N)] #//Extract (L*N) for one variable
    G_b[,ii6] = -rowSums(Z_mat*temp_x)/N
  }

  G = cbind(G_b, G_d, G_r)
  cov_mat = solve(t(G)%*%omega_hat%*%G)/N

  return(cov_mat)
}



##......................................linearity test.................................//


sup_wald_cal_multithreshold <- function(y_mat_rep,Z_mat,x_temp1,rx_temp2,rx_temp3,
                                        q_temp1,q_temp2,im,
                                        Gmat,N,W_n){



  k_1 = ncol(x_temp1)/N
  k_rx <- ncol(rx_temp2)/N
  dlt_len <- k_rx * im
  grid_num <- ncol(Gmat)

  cov_wald_rec = matrix(NA,grid_num*dlt_len,dlt_len)
  #// To record sigma_hat for each $delta$

  #// 1st-step GMM estimate


  g_1n_bar = rowSums(Z_mat*y_mat_rep)/N
  #//
  wald_vec = matrix(NA,grid_num,1)

  for(bb in 1:grid_num){
    #// Specially,omit the indentation for this loop

    gamma_b <- Gmat[,bb]
    oGc <- jump_GMM_est_Fun(y_mat_rep,g_1n_bar,gamma_b,Z_mat,x_temp1,rx_temp2,rx_temp3,
                            q_temp1,q_temp2,W_n)
    ep_collect_b = oGc$ep_collect; #GMM_est = oGc$GMM_est
    W_n_2_b = GMM_W_n_2_con(ep_collect_b, Z_mat)
    oGc2 <- jump_GMM_est_Fun(y_mat_rep,g_1n_bar,gamma_b,Z_mat,x_temp1,rx_temp2,rx_temp3,
                             q_temp1,q_temp2,W_n_2_b)
    ep_collect_2_b = oGc2$ep_collect; GMM_est_2_b = oGc2$GMM_est
    omega_hat_b = GMM_W_n_2_con(ep_collect_2_b, Z_mat)
    V_s = cov_con_Multithreshold_b(gamma_b,Z_mat,x_temp1,rx_temp2,rx_temp3,
                                   q_temp1,q_temp2,omega_hat_b)

    size_c = ncol(V_s) #// size of (V_s'* V_s)
    dlt_hat = as.matrix(GMM_est_2_b[(k_1+1):size_c])

    sel_m = matrix(0,dlt_len,size_c)
    sel_m[,(k_1+1):(size_c)] <- diag(dlt_len)

    cov_wald = try(sel_m %*% solve(t(V_s) %*% V_s) %*% t(sel_m),silent=TRUE)
    if(is(cov_wald,"try-error")){
      cov_wald = sel_m %*% ginv(t(V_s) %*% V_s) %*% t(sel_m)
    }

    wald_b = try(N*t(dlt_hat) %*% solve(cov_wald) %*% dlt_hat,silent = TRUE)
    if(is(wald_b,"try-error")){
      wald_b = N*t(dlt_hat) %*% ginv(cov_wald) %*% dlt_hat
    }

    wald_vec[bb] = wald_b

    cov_wald_rec[((bb-1)*dlt_len+1):(bb*dlt_len),] = cov_wald

  }

  sup_wald = max(wald_vec,na.rm=TRUE)

  return(list(sup_wald=sup_wald,cov_wald_rec=cov_wald_rec))

}



sup_wald_cal_multikink <- function(y_mat_rep,Z_mat,x_temp1,rx_temp2,
                                   rx_temp3,q_temp1,q_temp2,im,Gmat,
                                   N,W_n){

  k_1 = ncol(x_temp1)/N
  grid_num <- ncol(Gmat)

  cov_wald_rec = matrix(NA,grid_num*im,im)
  #// To record sigma_hat for each $r$

  #// 1st-step GMM estimate
  g_1n_bar = rowSums(Z_mat*y_mat_rep)/N
  #//

  wald_vec = rep(NA,grid_num)

  for(ii in 1:grid_num){
    gamma <- Gmat[,ii]

    oGc = try(kink_GMM_est_Fun(y_mat_rep,g_1n_bar,gamma,Z_mat,x_temp1,rx_temp2,rx_temp3,
                               q_temp1,q_temp2,W_n),silent=TRUE)
    if(is(oGc,"try-error")) next()
    GMM_est = oGc$GMM_est; g_2n_bar = oGc$g_2n_bar; g_n_bar= oGc$g_n_bar
    ep_collect = ep_cal_multikink(gamma, x_temp1, rx_temp2, rx_temp3, q_temp1, q_temp2,
                                  y_mat_rep, GMM_est)

    W_n_2 = GMM_W_n_2_con(ep_collect, Z_mat)

    oGc1 = try(kink_GMM_est_Fun(y_mat_rep,g_1n_bar,gamma,Z_mat,x_temp1,rx_temp2,rx_temp3,
                                q_temp1,q_temp2,W_n_2),silent=TRUE)
    if(is(oGc1,"try-error"))  next()
    GMM_est_2 = oGc1$GMM_est;g_2n_bar = oGc1$g_2n_bar; g_n_bar= oGc1$g_n_bar
    ep_collect_2 = ep_cal_multikink(gamma, x_temp1, rx_temp2, rx_temp3, q_temp1, q_temp2,
                                    y_mat_rep, GMM_est_2)
    omega_hat = GMM_W_n_2_con(ep_collect_2, Z_mat)


    ##
    V_s = cov_con_multikink_b(gamma,Z_mat,x_temp1,rx_temp2, rx_temp3,
                              q_temp1, q_temp2,omega_hat)

    size_c = ncol(V_s) #// size of (V_s'* V_s)
    #//verify size_c = k_1+kink_num
    dlt_hat = as.matrix(GMM_est_2[(k_1+1):size_c])

    sel_m = matrix(0,im,size_c)
    sel_m[,(size_c-im+1):size_c] <- diag(im)

    cov_wald = try(sel_m %*% solve(t(V_s) %*% V_s) %*% t(sel_m),silent=TRUE)
    if(is(cov_wald,"try-error")){
      cov_wald = sel_m %*% ginv(t(V_s) %*% V_s) %*% t(sel_m)
    }

    wald_b = try(N*t(dlt_hat) %*% solve(cov_wald) %*% dlt_hat,silent = TRUE)
    if(is(wald_b,"try-error")){
      wald_b = N*t(dlt_hat) %*% ginv(cov_wald) %*% dlt_hat
    }

    wald_vec[ii] = wald_b

    cov_wald_rec[((ii-1)*im+1):(ii*im),] = cov_wald


  }

  sup_wald = max(wald_vec,na.rm=TRUE)

  return(list(sup_wald=sup_wald,cov_wald_rec=cov_wald_rec))

}




cov_con_Multithreshold_b <- function(gamma,Z_mat,x_temp1,rx_temp2,rx_temp3,
                                     q_temp1,q_temp2,omega_hat)
{
  L = nrow(Z_mat)
  N = ncol(Z_mat)
  k_1 = ncol(x_temp1)/N
  k_rx = ncol(rx_temp2)/N
  im <- length(gamma)
  ind_temp2_list <- ind_temp3_list <- list()
  for(kk in 1:im){
    ind_temp2_list[[kk]] <- (q_temp1 > gamma[kk])
    ind_temp3_list[[kk]] <- (q_temp2 > gamma[kk])
  }

  G_b = matrix(0,L,k_1)
  G_d = matrix(0,L,k_rx*im)

  ##covariance matrix
  for(ii6 in 1:k_1){
    temp_x = x_temp1[,((ii6-1)*N+1):(ii6*N)]
    G_b[,ii6] = -rowSums(Z_mat*temp_x)/N
  }

  for(kk in 1:im){
    for(iix in 1:k_rx){
      G_d[,((kk-1)*k_rx+iix)] <- rowSums(Z_mat*(rx_temp2[,((iix-1)*N+1):(iix*N)]*ind_temp2_list[[kk]]-
                                                  rx_temp3[,((iix-1)*N+1):(iix*N)]*ind_temp3_list[[kk]]))/N
    }
  }


  V_s = (chol(omega_hat)) %*% cbind(G_b,G_d)


  return(V_s)

}

cov_con_multikink_b <- function(gamma,Z_mat,
                                x_temp1,rx_temp2, rx_temp3,
                                q_temp1, q_temp2,omega_hat){

  L = nrow(Z_mat)
  N = ncol(Z_mat)
  k_1 = ncol(x_temp1)/N
  kink_num <- length(gamma)
  ind_temp2_list <- ind_temp3_list <- list()
  for(kk in 1:kink_num){
    ind_temp2_list[[kk]] <- (q_temp1 > gamma[kk])
    ind_temp3_list[[kk]] <- (q_temp2 > gamma[kk])
  }

  G_b = matrix(0,L, k_1)
  G_d = matrix(0,L, kink_num)
  for(ii6 in 1:k_1){
    temp_x = x_temp1[,((ii6-1)*N+1):(ii6*N)] #//Extract (L*N) for one variable
    G_b[,ii6] = -rowSums(Z_mat*temp_x)/N
  }
  for(kk in 1:kink_num){
    G_d[,kk] <- rowSums(Z_mat*((q_temp1-gamma[kk])*ind_temp2_list[[kk]]-
                                 (q_temp2-gamma[kk])*ind_temp3_list[[kk]]))/N
  }
  V_s = (chol(omega_hat)) %*% cbind(G_b,G_d)
  return(V_s)
}






sup_wald_fb_multithreshold <- function(y_mat_rep,Z_mat,x_temp1,rx_temp2,rx_temp3,
                                       q_temp1,q_temp2,ep_collect_2,W_n_2,cov_wald_rec,
                                       Gmat,N)
{
  k_1 = ncol(x_temp1)/N
  k_rx <- ncol(rx_temp2)/N
  grid_num <- ncol(Gmat)
  im <- nrow(Gmat)
  dlt_len <- im*k_rx
  L = nrow(W_n_2)

  norm_vec = matrix(qnorm(runif(N,0,1)),nrow=1)
  b_eta <- kronecker(matrix(1,L,1),norm_vec)  #kronecker(matrix(1,1,k_endo), matrix(c(0:(T-2))))

  g_1n_bar_b = rowSums(Z_mat*ep_collect_2*b_eta)/N
  #// This part makes difference for fast bootstrap

  wald_vec = rep(NA,grid_num)

  for(bb in 1:grid_num){
    gamma_b <- Gmat[,bb]
    oGc <- try(jump_GMM_est_Fun(y_mat_rep,g_1n_bar_b,gamma_b,Z_mat,x_temp1,rx_temp2,rx_temp3,
                                q_temp1,q_temp2,W_n_2,flag =0),silent=TRUE)
    if(is(oGc,"try-error")) next
    GMM_est_b = as.matrix(oGc$GMM_est)
    size_c = nrow(GMM_est_b)
    dlt_hat = as.matrix(GMM_est_b[(k_1+1):size_c,])
    cov_wald = cov_wald_rec[((bb-1)*dlt_len+1):(bb*dlt_len),]

    wald_b = try(N*t(dlt_hat) %*% solve(cov_wald) %*% dlt_hat,silent = TRUE)
    if(is(wald_b,"try-error")){
      wald_b = try(N*t(dlt_hat) %*% ginv(cov_wald) %*% dlt_hat,silent = TRUE)
    }
    wald_vec[bb] = wald_b

  }

  sup_wald_b = max(wald_vec,na.rm=TRUE) #

  return(sup_wald_b)


}



sup_wald_fb_multikink <- function(y_mat_rep,Z_mat,x_temp1,rx_temp2,rx_temp3,
                                  q_temp1,q_temp2,ep_collect_2,W_n_2,cov_wald_rec,
                                  Gmat,N)
{
  k_1 = ncol(x_temp1)/N
  L = nrow(W_n_2)
  grid_num <- ncol(Gmat)
  im <- nrow(Gmat)
  norm_vec = matrix(qnorm(runif(N,0,1)),nrow=1)
  b_eta <- kronecker(matrix(1,L,1),norm_vec)  #kronecker(matrix(1,1,k_endo), matrix(c(0:(T-2))))

  g_1n_bar = rowSums(Z_mat*ep_collect_2*b_eta)/N
  #// This part makes difference for fast bootstrap

  wald_vec = rep(NA,grid_num)

  for(ii in 1:grid_num){
    gamma <- Gmat[,ii]

    oGc <- try(kink_GMM_est_Fun(y_mat_rep,g_1n_bar,gamma,Z_mat,x_temp1,rx_temp2,rx_temp3,
                                q_temp1,q_temp2,W_n_2,flag =0),silent=TRUE)
    if(is(oGc,"try-error")) next
    GMM_est = oGc$GMM_est; g_2n_bar = oGc$g_2n_bar; g_n_bar= oGc$g_n_bar

    size_c = nrow(GMM_est)
    dlt_hat = as.matrix(GMM_est[(k_1+1):size_c])

    cov_wald = cov_wald_rec[((ii-1)*im+1):(ii*im),]
    #// This is the same as the cov_wald for original wald test
    wald_b = N*t(dlt_hat) %*% solve(cov_wald) %*% dlt_hat
    wald_vec[ii] = wald_b

  }


  sup_wald_b = max(wald_vec,na.rm=TRUE)

  return(sup_wald_b)


}






#//.............The following remains as before........................\\



GMM_weight <- function(N,T,L,moment_mat,flag_static)
{
  #if(flag_static == 1){
  #  L = sum(moment_mat) - sum(moment_mat[,2])
  # In case of static model
  #}else{
  #  L = sum(moment_mat[2:(nrow(moment_mat)),])
  #}

  eps_mat <- matrix(rnorm(N*T,0,1),nrow=N,ncol=T)
  eps_mat_fd = FD_con(eps_mat,T)

  flag_r1 = 1
  weight_mat_rep = matrix(0,L,N)
  if(flag_static == 1){
    for(kk in 1:(T-1)){
      m_temp = sum(moment_mat[kk,]) - sum(moment_mat[kk,2])
      # of stacked moment condition for t = kk
      if(m_temp == 0){
        next()
      }else{
        weight_mat_rep[flag_r1:(flag_r1+m_temp-1),] = t(kronecker(eps_mat_fd[,kk], matrix(1,1,m_temp)))
        flag_r1 = flag_r1+m_temp
      }
    }
  }else{
    ## For dynamic model (Default)
    for(kk in 2:(T-1)){
      m_temp = sum(moment_mat[kk,]) #// # of stacked moment condition for t = kk
      if(m_temp == 0){
        next()
      }else{
        weight_mat_rep[flag_r1:(flag_r1+m_temp-1),]=t(kronecker(eps_mat_fd[,kk], matrix(1,1,m_temp)))
        flag_r1 = flag_r1+m_temp
      }
    }
  }

  return(weight_mat_rep)


}



