library(MASS)
library(splines)
library(truncnorm)
library(survival)



#### data generation ##########################################################
generate_data <- function(n,cen,model){
  if(model=='OC'){  
    x1 <- rbinom(n,1,0.5)
    x2 <- runif(n,-1,1)
    ## parameter for h(X)
    beta0 <- c(0,0.2,0.2)
    X <- cbind(1,x1,x2)
    theta0 <- c(-0.15,0.3,0.942)
    ## exposure
    #gamma0 <- c(0.25,-0.5,-0.25)
    gamma0 <- c(0.1,-0.2,-0.2)
    A <- rtruncnorm(n,a=-1,b=1,mean=X %*% gamma0,sd=0.5)
    #A <- runif(n,-1,1)
    tau <- 1.2
    phiA <- tau/(1+exp(-5*A))-tau/2
    tmp <- runif(n)
    #t_time = -log(tmp)/exp(X %*% beta0 + (A * tau)*I(X %*% theta0 >= 0) )
    t_time = -log(tmp)/exp(X %*% beta0 + (phiA)*I(X %*% theta0 >= 0) )
    c_time <- runif(n,0,cen)
    obs_time = pmin(c_time,t_time)
    delta = as.numeric(t_time<=c_time)
    data <- data.frame(obs_time,delta,X,A)
    colnames(data) <- c('obs_time','delta','x0','x1','x2','A')
    data <- data[order(data$A),]
  }
  if(model=='OF'){
    x1 <- rbinom(n,1,0.5)
    x2 <- runif(n,-1,1)
    ## parameter for h(X)
    beta0 <- c(0.1,0,-0.3)
    X <- cbind(1,x1,x2)
    Xtmp <- cbind(x1,x2,x2^2)
    theta0 <- c(-0.15,0.3,0.942)
    ## exposure
    gamma0 <- c(0.1,-0.2,-0.2)
    A <- rtruncnorm(n,a=-1,b=1,mean=X %*% gamma0,sd=0.5)
    #A <- runif(n,-1,1)
    tau <- 1.2
    phiA <- tau/(1+exp(-5*A))-tau/2
    tmp <- runif(n)
    #t_time = -log(tmp)/exp(Xtmp %*% beta0 + (A * tau)*I(X %*% theta0 >= 0) )
    t_time = -log(tmp)/exp(Xtmp %*% beta0 + (phiA)*I(X %*% theta0 >= 0) )
    c_time <- runif(n,0,cen)
    obs_time = pmin(c_time,t_time)
    delta = as.numeric(t_time<=c_time)
    data <- data.frame(obs_time,delta,X,A)
    colnames(data) <- c('obs_time','delta','x0','x1','x2','A')
    data <- data[order(data$A),]
  }
  return(data)
}

## estimate loglikelihood parameter
cumfun1 <- function(x,sigma_ini){
  A_new <- seq(-1,1,length.out = 2000)
  int_tmp <- dnorm(A_new,x,sigma_ini)* (A_new-x) 
  return(sum((int_tmp[-length(A_new)]+int_tmp[-1])/2 * diff(A_new)))
}

cumfun2 <- function(x,sigma_ini){
  A_new <- seq(-1,1,length.out = 2000)
  int_tmp <- dnorm(A_new,x,sigma_ini)* (A_new-x)^2 
  return(sum((int_tmp[-length(A_new)]+int_tmp[-1])/2 * diff(A_new)))
}

cumfun3 <- function(x,sigma_ini){
  A_new <- seq(-1,1,length.out = 2000)
  int_tmp <- dnorm(A_new,x,sigma_ini)* (A_new-x)^3 
  return(sum((int_tmp[-length(A_new)]+int_tmp[-1])/2 * diff(A_new)))
}

cumfun4 <- function(x,sigma_ini){
  A_new <- seq(-1,1,length.out = 2000)
  int_tmp <- dnorm(A_new,x,sigma_ini)* (A_new-x)^4
  return(sum((int_tmp[-length(A_new)]+int_tmp[-1])/2 * diff(A_new)))
}



newton <- function(gamma_ini,sigma_ini,iter_max=1000,tol=0.0001){
  for(j in 1:iter_max){
    mu_hat <- X_gps %*% gamma_ini
    diff_A <- data$A-mu_hat
    denom <- pnorm(1,mu_hat,sigma_ini)-pnorm(-1,mu_hat,sigma_ini)
    cum1 <- apply(matrix(mu_hat,ncol=1),1,cumfun1,sigma_ini=sigma_ini)
    cum2 <- apply(matrix(mu_hat,ncol=1),1,cumfun2,sigma_ini=sigma_ini)
    cum4 <- apply(matrix(mu_hat,ncol=1),1,cumfun4,sigma_ini=sigma_ini)
    
    ## score for gamma
    g_score <- colMeans(X_gps * c(diff_A/sigma_ini^2 - cum1/sigma_ini^2 /denom))
    
    ## score for sigma^2
    s1 <- cum2/sigma_ini^4-denom/sigma_ini^2
    s_score <- mean(diff_A^2/sigma_ini^4 - 1/sigma_ini^2 - s1/denom ) * 0.5
    
    ## Hessian for gamma
    gH1 <- -1/sigma_ini^2
    gH2 <- -(s1*denom-(cum1/sigma_ini^2)^2) / denom^2
    g_H <- crossprod(X_gps,X_gps * c(gH1+gH2)) / n
    
    ## Hessian for sigma^2
    sH1 <- -diff_A^2/sigma_ini^6
    sH2 <- 1/(2*sigma_ini^4)
    sH3_tmp1 <- (cum4/(4*sigma_ini^8)-cum2/(5*sigma_ini^6)+denom/(2*sigma_ini^4))
    sH3 <- -(sH3_tmp1*denom - 0.25*cum2/sigma_ini^4 * s1)/denom^2
    s_H <- mean(sH1+sH2+sH3)
    
    gamma_new <- gamma_ini - solve(g_H)%*%g_score
    sigma_new <- sigma_ini - 1/s_H * s_score
    if(prod(abs(gamma_new-gamma_ini)<tol)==T){break}
    gamma_ini <- gamma_new
    sigma_ini <- sigma_new
    #print(c(gamma_new,sigma_new))
  }
  return(drop(c(gamma_new,sigma_new)))
}

## calculate expectation in Randomized trial
expect_fun <- function(A_new,sp_new,den){
  tmp <- sp_new * den
  tmp_med <- (tmp[-1,]+tmp[-length(A_new),])/2
  return(colSums(tmp_med * diff(A_new)))
}

## calculate expectation in Observational study
expect_fun_ob <- function(mu,sigma_ini){
  den <- dtruncnorm(A_new,-1,1,mean=mu,sd=sigma_ini)
  tmp <- sp_new * den
  tmp_med <- (tmp[-1,]+tmp[-length(A_new),])/2
  return(colSums(tmp_med * diff(A_new)))
}

expect_fun_ob1 <- function(mu,sigma_ini){
  den <- dtruncnorm(A_new,-1,1,mean=mu,sd=sigma_ini) * (A_new-mu)/sigma_ini^2
  tmp <- sp_new * den
  tmp_med <- (tmp[-1,]+tmp[-length(A_new),])/2
  return(colSums(tmp_med * diff(A_new)))
}

expect_fun_ob2 <- function(mu,sigma_ini){
  den <- dtruncnorm(A_new,-1,1,mean=mu,sd=sigma_ini) * (A_new-mu)^2/(2*sigma_ini^4)
  tmp <- sp_new * den
  tmp_med <- (tmp[-1,]+tmp[-length(A_new),])/2
  return(colSums(tmp_med * diff(A_new)))
}

set.seed(2021)
n <- 2000
out <- matrix(NA,2000,1)
for(i in 1:nrow(out)){
  data <- generate_data(n,11.43,'OF')
  out[i,] <- mean(data$delta)
}
1-mean(out)

### start simulation (Type-I error) ###########################################
set.seed(as.numeric(args[1]) + 2020)
out <- matrix(NA,33,200)
n <- 2000

angle1 <- seq(0,pi,length.out = 50)
angle2 <- seq(0,2*pi,length.out = 100)
theta_save <- NA
for(a in 1:length(angle1)){
  for(b in 1:length(angle2)){
    theta_save <- rbind(theta_save,
                        c(cos(angle1[a]),sin(angle1[a])*cos(angle2[b]),sin(angle1[a])*sin(angle2[b])))
  }
}
theta_save <- theta_save[-1,]

num.theta <- nrow(theta_save)
num.pert <- 1000
i <- 1
while(i <= ncol(out)){
  #### generate data
  data <- generate_data(n,5.12,'OC')
  X <- as.matrix(data[,3:5])
  
  #### fit h model using Cox regression
  #model_mu <- coxph(Surv(obs_time, delta) ~ x1+x2,data=data)
  model_mu <- coxph(Surv(obs_time, delta) ~ x1+I(x2^2),data=data)
  beta_est <- coef(model_mu)
  X0 <- model.matrix(model_mu)
  ty <- residuals(model_mu,type = 'martingale')
 

  #### estimated expectation of B(A)
  X_gps <- X[,1:2]
  est <- newton(gamma_ini = c(0.25,-0.5),
                sigma_ini = 0.5)
  gamma_ini <- est[1:2]
  sigma_ini <- est[3]
  
  A_new <- data$A
  sp_new <- sweep(bs(c(0,A_new),degree = 2,df=5,Boundary.knots=c(-1,1))[-1,],
                  2,
                  bs(c(0,A_new),degree = 2,df=5,Boundary.knots=c(-1,1))[1,])
  
  mu_hat <- X_gps %*% gamma_ini
  diff_A <- data$A-mu_hat
  denom <- pnorm(1,mu_hat,sigma_ini)-pnorm(-1,mu_hat,sigma_ini)
  cum1 <- apply(matrix(mu_hat,ncol=1),1,cumfun1,sigma_ini = sigma_ini)
  cum2 <- apply(matrix(mu_hat,ncol=1),1,cumfun2,sigma_ini = sigma_ini)
  cum3 <- apply(matrix(mu_hat,ncol=1),1,cumfun3,sigma_ini = sigma_ini)
  cum4 <- apply(matrix(mu_hat,ncol=1),1,cumfun4,sigma_ini = sigma_ini)
  expect <- t(apply(matrix(mu_hat,ncol = 1),1,expect_fun_ob,sigma_ini=sigma_ini))

  ## E(B(A)|X_i), in randomized trial, f(A|X_i) is given
  sp <- sp_new
  ta <- sp - expect
  

  #### create table for test statistics
  teststat <- -Inf
  sup_ts <- 0
  teststat_p <- matrix(NA,nrow = num.pert,ncol=num.theta)
  prop <- 0
  
  teststat2 <- -Inf
  sup_ts2 <- 0
  teststat_p2 <- matrix(NA,nrow = num.pert,ncol=num.theta)
  
  #### generate perturbation r.v.
  G <- matrix(rnorm(n*num.pert),nrow = n,ncol = num.pert)
  
  #### check each theta
  for(k in 1:num.theta){
    subgroup <- drop(X %*% theta_save[k,]) >= 0.0
    if( sum(subgroup)<= (0.1*n) | sum(subgroup)>= (0.9*n) ) next
    prop <- prop + 1
    
    spsi1 <- ta * c(ty * subgroup)
    psi_t <- spsi1 
    V_t <- ginv(crossprod(psi_t,psi_t))
    
    temp <- drop(colSums(spsi1) %*% V_t %*% colSums(spsi1))
    temp2 <- drop(colSums(spsi1) %*% colSums(spsi1))
    
    if( temp > teststat ) {
      teststat <- temp
      sup_ts <- k
    }
    
    if( temp2 > teststat2 ) {
      teststat2 <- temp2
      sup_ts2 <- k
    }
    
    pert_tmp <- t(psi_t) %*% G
    teststat_p[,k] <- apply(pert_tmp,2,function(x){drop(x %*% V_t %*% x)})
    teststat_p2[,k] <- apply(pert_tmp,2,function(x){drop(x  %*% x)})
  }
  
  dissup_ts <- apply(X = teststat_p,
                     MARGIN = 1,
                     FUN = max,
                     na.rm = TRUE)
  p_value <- sum(dissup_ts >= teststat) / length(dissup_ts)
  
  dissup_ts2 <- apply(X = teststat_p2,
                      MARGIN = 1,
                      FUN = max,
                      na.rm = TRUE)
  p_value2 <- sum(dissup_ts2 >= teststat2) / length(dissup_ts2)
  
  ## check if reject null hypothesis
  out[1,i] <- p_value
  out[2,i] <- p_value2
  out[3,i] <- sum(drop(X %*% theta_save[sup_ts,]) >= 0.0)
  out[4,i] <- sum(drop(X %*% theta_save[sup_ts2,]) >= 0.0)
  out[5,i] <- teststat
  out[6,i] <- teststat2

  subgroup <- drop(X %*% theta_save[sup_ts,]) >= 0.0
  spsi1 <- ta * c(ty * subgroup)
  psi_t <- spsi1 
  V_t <- ginv(crossprod(psi_t,psi_t))
  #V_t <- ginv(crossprod(  sweep(psi_t,2,colMeans(psi_t)) , sweep(psi_t,2,colMeans(psi_t)) ))
  
  #out[7:15,i] <- c(V_t[1,],V_t[2,],V_t[3,])
  out[16:18,i] <- theta_save[sup_ts,]
  out[19:21,i] <- theta_save[sup_ts2,]
  theta0 <- c(-0.15,0.3,0.942)
  out[22,i] <- mean(abs( ((X %*% theta_save[sup_ts,])>=0) - ((X %*% theta0)>=0) ))
  out[23,i] <- mean(abs( ((X %*% theta_save[sup_ts2,])>=0) - ((X %*% theta0)>=0) ))
  
  #filename <- paste0('OC_power_lin0.25_cen20_10perboth_',n,'_',as.numeric(args[1]),'.csv')
  #filename <- paste0('OC_power_misBase_noexp_lin0.8_cen20_10perboth_',n,'_',as.numeric(args[1]),'.csv')
  filename <- paste0('OC_power_misGPS_noexp_logit1.2_cen20_10perboth_',n,'_',as.numeric(args[1]),'.csv')
  write.csv(out,filename)
  
  print(i)
  print(out[,i])
  i <- i+1
}


