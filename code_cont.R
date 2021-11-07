library(MASS)
library(splines)
library(truncnorm)

#### data generation ##########################################################
generate_data <- function(n,model){
  if(model=='RC'){
    x1 <- rbinom(n,1,0.5)
    x2 <- runif(n,-1,1)
    ## parameter for h(X)
    beta0 <- c(1,1,1)
    X <- cbind(1,x1,x2)
    Y <- X%*%beta0 + rnorm(n,0,0.5)
    ## exposure
    A <- runif(n,-1,1)
    data <- data.frame(Y,X,A)
    colnames(data) <- c('Y','x0','x1','x2','A')
    data <- data[order(data$A),]
  }
  if(model=='RF'){
    x1 <- rbinom(n,1,0.5)
    x2 <- runif(n,-1,1)
    ## parameter
    beta0 <- c(1,0.5,0,0.3)
    X <- cbind(1,x1,x2)
    Xtmp <- cbind(1,x1,x2,x2^2)
    Y <- Xtmp%*%beta0 + rnorm(n,0,0.5)
    ## exposure
    A <- runif(n,-1,1)
    data <- data.frame(Y,X,A)
    colnames(data) <- c('Y','x0','x1','x2','A')
    data <- data[order(data$A),]
  }
  return(data)
}

expect_fun <- function(A_new,sp_new,den){
  tmp <- sp_new * den
  tmp_med <- (tmp[-1,]+tmp[-length(A_new),])/2
  return(colSums(tmp_med * diff(A_new)))
}


### start simulation (Type-I error) ###########################################
set.seed(as.numeric(args[1]) + 2020)
out <- matrix(NA,26,500)
n <- 4000

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
  data <- generate_data(n,'RC')
  X <- as.matrix(data[,2:4])
  
  #### fit h model
  model_mu <- lm(Y~x1+x2,data=data)
  beta_est <- coef(model_mu)
  mu_fit <- fitted(model_mu)
  X0 <- model.matrix(model_mu)
  ty <- data$Y - mu_fit
  
  #### estimated expectation of B(A)
  A_new <- data$A
  sp_new <- sweep(bs(c(0,A_new),degree = 2,df=5,Boundary.knots=c(-1,1))[-1,],
                  2,
                  bs(c(0,A_new),degree = 2,df=5,Boundary.knots=c(-1,1))[1,])
  expect <- expect_fun(A_new,sp_new,0.5)
  sp <- sp_new
  ta <- sweep(sp,2,expect)
  
  #### psi2: score for beta
  spsi2 <- t(X0 * c(ty))
  
  #### C1 in the paper
  dpsi22 <- -crossprod(X0, X0) / n
  
  res1 <- solve(dpsi22) %*% spsi2
  
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
    
    dpsi12 <- -crossprod(subgroup * ta,X0) / n
    
    psi_t <- spsi1 - t(dpsi12 %*% res1)
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
  dpsi12 <- -crossprod(subgroup * ta,X0) / n
  
  psi_t <- spsi1 - t(drop(dpsi12 %*% res1)) 
  V_t <- ginv(crossprod(psi_t,psi_t))
  
  out[16:18,i] <- theta_save[sup_ts,]
  out[19:21,i] <- theta_save[sup_ts2,]
  theta0 <- c(-0.15,0.3,0.942)
  
  #filename <- paste0('RC_type1_df5_newV_10perboth_',n,'_',as.numeric(args[1]),'.csv')
  #write.csv(out,filename)
  
  print(i)
  print(out[,i])
  i <- i+1
}


