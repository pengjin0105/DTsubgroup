library(MASS)
library(splines)
library(truncnorm)
library(randomForestSRC)
library(survival)

args=commandArgs(TRUE)

#### data generation ##########################################################
generate_data <- function(n,cen,model){
  if(model=='RF'){
    x1 <- rbinom(n,1,0.5)
    x2 <- runif(n,-1,1)
    ## parameter for h(X)
    beta0 <- c(0.1,0,-0.3)
    X <- cbind(1,x1,x2)
    Xtmp <- cbind(x1,x2,x2^2)
    theta0 <- c(-0.15,0.3,0.942)
    ## exposure
    A <- runif(n,-1,1)
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


expect_fun <- function(A_new,sp_new,den){
  tmp <- sp_new * den
  tmp_med <- (tmp[-1,]+tmp[-length(A_new),])/2
  return(colSums(tmp_med * diff(A_new)))
}

## estimate loglikelihood parameter
cumfun1 <- function(x,sigma_ini){
  A_new <- seq(-1,1,length.out = 500)
  int_tmp <- dnorm(A_new,x,sigma_ini)* (A_new-x) 
  return(sum((int_tmp[-length(A_new)]+int_tmp[-1])/2 * diff(A_new)))
}

cumfun2 <- function(x,sigma_ini){
  A_new <- seq(-1,1,length.out = 500)
  int_tmp <- dnorm(A_new,x,sigma_ini)* (A_new-x)^2 
  return(sum((int_tmp[-length(A_new)]+int_tmp[-1])/2 * diff(A_new)))
}

cumfun3 <- function(x,sigma_ini){
  A_new <- seq(-1,1,length.out = 500)
  int_tmp <- dnorm(A_new,x,sigma_ini)* (A_new-x)^3 
  return(sum((int_tmp[-length(A_new)]+int_tmp[-1])/2 * diff(A_new)))
}

cumfun4 <- function(x,sigma_ini){
  A_new <- seq(-1,1,length.out = 500)
  int_tmp <- dnorm(A_new,x,sigma_ini)* (A_new-x)^4
  return(sum((int_tmp[-length(A_new)]+int_tmp[-1])/2 * diff(A_new)))
}



newton <- function(gamma_ini,sigma_ini,iter_max=1000,tol=0.0001){
  for(j in 1:iter_max){
    mu_hat <- X_gps%*%gamma_ini
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

get_index <- function(ind.surv.time){
  sum(ind.surv.time > data$obs_time[data$delta==1])
}


### start simulation (Type-I error) ###########################################
set.seed(as.numeric(args[1]) + 2020)
out <- matrix(NA,28,100)
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
  data <- generate_data(n,5.48,'RF')
  X <- as.matrix(data[,3:5])
  
  #### fit u model: misBase
  model.misBase <- coxph(Surv(obs_time, delta) ~ x1+x2,data=data)
  beta_est.misBase <- coef(model.misBase)
  X0.misBase <- model.matrix(model.misBase)
  ty.misBase <- residuals(model.misBase,type = 'martingale')
  
  #### fit u model: correct
  model.true <- coxph(Surv(obs_time, delta) ~ x1+I(x2^2),data=data)
  beta_est.true <- coef(model.true)
  X0.true <- model.matrix(model.true)
  ty.true <- residuals(model.true,type = 'martingale')
  
  #### fit u model: nonparametric
  model.rf <- rfsrc(Surv(obs_time, delta) ~ x1+x2,data=data,ntree = 200,mtry = ncol(X))
  pred.event <- -log(model.rf$survival.oob)
  ind <- as.numeric(sapply(data$obs_time,get_index))
  comp <- as.numeric(mapply(function(i,j)pred.event[i,j],seq(n),ind))
  comp[is.na(comp)] <- 0
  ty.rf <- data$delta-comp
  ty.rf[ty.rf== -Inf] <- 0
  
  
  #### estimated expectation of B(A): true parametric
  X_gps <- X
  A_new <- data$A
  sp_new <- sweep(bs(c(0,A_new),degree = 2,df=5,Boundary.knots=c(-1,1))[-1,],
                  2,
                  bs(c(0,A_new),degree = 2,df=5,Boundary.knots=c(-1,1))[1,])
  expect <- expect_fun(A_new,sp_new,0.5)
  sp <- sp_new
  ta.para <- sweep(sp,2,expect)

  
  #### create table for test statistics
  teststat.misBase <- -Inf
  sup_ts.misBase <- 0
  teststat_p.misBase <- matrix(NA,nrow = num.pert,ncol=num.theta)
  prop <- 0
  
  teststat.true <- -Inf
  sup_ts.true <- 0
  teststat_p.true <- matrix(NA,nrow = num.pert,ncol=num.theta)
  
  teststat.rf <- -Inf
  sup_ts.rf <- 0
  teststat_p.rf <- matrix(NA,nrow = num.pert,ncol=num.theta)
  
  #### generate perturbation r.v.
  G <- matrix(rnorm(n*num.pert),nrow = n,ncol = num.pert)
  
  #### check each theta
  for(k in 1:num.theta){
    subgroup <- drop(X %*% theta_save[k,]) >= 0.0
    if( sum(subgroup)<= (0.1*n) | sum(subgroup)>= (0.9*n) ) next
    prop <- prop + 1
    
    
    spsi1.misBase <- ta.para * c(ty.misBase * subgroup)
    psi_t.misBase <- spsi1.misBase #- t(dpsi12.misBase %*% res1.misBase)
    V_t.misBase <- ginv(crossprod(psi_t.misBase,psi_t.misBase))
    temp.misBase <- drop(colSums(spsi1.misBase) %*% V_t.misBase %*% colSums(spsi1.misBase))
    
    spsi1.true <- ta.para * c(ty.true * subgroup)
    psi_t.true <- spsi1.true #- t(dpsi12.true %*% res1.true)
    V_t.true <- ginv(crossprod(psi_t.true,psi_t.true))
    temp.true <- drop(colSums(spsi1.true) %*% V_t.true %*% colSums(spsi1.true))
    
    spsi1.rf <- ta.para * c(ty.rf * subgroup)
    psi_t.rf <- spsi1.rf
    V_t.rf <- ginv(crossprod(psi_t.rf,psi_t.rf))
    temp.rf <- drop(colSums(spsi1.rf) %*% V_t.rf %*% colSums(spsi1.rf))
    
    if( temp.misBase > teststat.misBase ) {
      teststat.misBase <- temp.misBase
      sup_ts.misBase <- k
    }
    
    if( temp.true > teststat.true ) {
      teststat.true <- temp.true
      sup_ts.true <- k
    }
    
    if( temp.rf > teststat.rf ) {
      teststat.rf <- temp.rf
      sup_ts.rf <- k
    }
    
    pert_tmp.misBase <- t(psi_t.misBase) %*% G
    teststat_p.misBase[,k] <- apply(pert_tmp.misBase,2,function(x){drop(x %*% V_t.misBase %*% x)})
    
    pert_tmp.true <- t(psi_t.true) %*% G
    teststat_p.true[,k] <- apply(pert_tmp.true,2,function(x){drop(x %*% V_t.true %*% x)})
  
    pert_tmp.rf <- t(psi_t.rf) %*% G
    teststat_p.rf[,k] <- apply(pert_tmp.rf,2,function(x){drop(x %*% V_t.rf %*% x)})
    
  }
  
  dissup_ts.misBase <- apply(X = teststat_p.misBase,
                       MARGIN = 1,
                       FUN = max,
                       na.rm = TRUE)
  p_value.misBase <- sum(dissup_ts.misBase >= teststat.misBase) / length(dissup_ts.misBase)
  
  
  dissup_ts.true <- apply(X = teststat_p.true,
                          MARGIN = 1,
                          FUN = max,
                          na.rm = TRUE)
  p_value.true <- sum(dissup_ts.true >= teststat.true) / length(dissup_ts.true)
  
  dissup_ts.rf <- apply(X = teststat_p.rf,
                          MARGIN = 1,
                          FUN = max,
                          na.rm = TRUE)
  p_value.rf <- sum(dissup_ts.rf >= teststat.rf) / length(dissup_ts.rf)
  
  
  ## check if reject null hypothesis
  out[1,i] <- p_value.misBase
  out[2,i] <- p_value.true
  out[3,i] <- p_value.rf
  out[4,i] <- sum(drop(X %*% theta_save[sup_ts.misBase,]) >= 0.0)
  out[5,i] <- sum(drop(X %*% theta_save[sup_ts.true,]) >= 0.0)
  out[6,i] <- sum(drop(X %*% theta_save[sup_ts.rf,]) >= 0.0)
  out[7,i] <- teststat.misBase
  out[8,i] <- teststat.true
  out[9,i] <- teststat.rf

  out[10:12,i] <- theta_save[sup_ts.misBase,]
  out[13:15,i] <- theta_save[sup_ts.true,]
  out[16:18,i] <- theta_save[sup_ts.rf,]
  theta0 <- c(-0.15,0.3,0.942)
  out[19,i] <- mean(abs( ((X %*% theta_save[sup_ts.misBase,])>=0) - ((X %*% theta0)>=0) ))
  out[20,i] <- mean(abs( ((X %*% theta_save[sup_ts.true,])>=0) - ((X %*% theta0)>=0) ))
  out[21,i] <- mean(abs( ((X %*% theta_save[sup_ts.rf,])>=0) - ((X %*% theta0)>=0) ))
  
  #filename <- paste0('RF_power_expU_all_misAsX1only_logit1.2_regV_10perboth_',n,'_',as.numeric(args[1]),'.csv')
  #write.csv(out,filename)
  
  print(i)
  print(out[,i])
  i <- i+1
}

