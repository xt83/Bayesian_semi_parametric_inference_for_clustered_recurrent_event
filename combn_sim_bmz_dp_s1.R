 # grace: setwd("/gpfs/loomis/home.grace/jc3625/project/competing_risk_simulation")678900
# farnam:

setwd("/gpfs/ysm/project/li_fan/xt83")
library(LaplacesDemon)

args <- commandArgs(trailingOnly = TRUE)
seed<- as.numeric(args[1])
set.seed(seed)

n<- 5400
num_hos<- 180
pr<- 3
pd<- 3

#decide which person belongs to which hospital
observed<- matrix(0,nrow=n,ncol=500)
observed[,1]<- seq(from=1,to=n,by=1)
observed[,2]<- rep(seq(from=1,to=num_hos,by=1),times=(n/num_hos))

#recurrent covariates zr
zr<- matrix(0,nrow=n,ncol=pr)
for (i in 1:pr){
  zr[,i]<- scale(rnorm(n,mean=0,sd=0.1))
}

observed[,(pd+6):(pd+pr+5)]<- zr

#survival covariates zd
zd<- matrix(0,nrow=n,ncol=pd)
zd[,1:2]<- zr[,2:3]
zd[,3]<- scale(rnorm(n,mean=0,sd=0.1))

observed[,6:(pd+5)]<- zd

#parameters
betar_true<- c(0.4,0.3,0.2)
kesi_true<- 1.5

betad_true<- c(0.2, 0.3, 0.4)
zeta_true<- 0.1
mu_intercept_true<- 0.15
tao_true<- -0.5

set<-rep(0,n)
for(i in 1:n){
  ind<-sample(1:4,1)
  set[i]<-ifelse(ind==1,1,0)*0.7+ifelse(ind==2,1,0)*2.2+ifelse(ind==3,1,0)*5.2+ifelse(ind==4,1,0)*8.2
}
sigma_true<-set


#latent indicator for cured:y 
num_eta<- 4
eta_true<- rep(1,times=num_eta)
logit_prob_true<- cbind(zr,zd[,3])%*%eta_true
prob_true<-  1/(1 + exp(-logit_prob_true))
y_true<- rbinom(n,1,prob_true)

#gamma_ij
epislon_true<- rep(0.5^2,times=num_hos)
gamma_true<- rep(0,times=n)
for (i in 1:num_hos){
  gamma_true[observed[,2]==i]<- rlnorm((n/num_hos),meanlog = 0,sdlog = sqrt(epislon_true[i]))
}

#mu_j: mixture normal distribution of 5 components
mu_hos_true<-rnormm(num_hos, p=rep(0.2,times=5), mu=seq(from=-0.4,to=0.4,by=0.2), sigma=rep(0.1,times=5)) 

#death time d
death_time<- rep(0,times=n)
scale_para<- rep(0,times=n)
for (i in 1:n){
  shape_para<- sigma_true[i]
  scale_para[i]<- exp(mu_intercept_true+zeta_true*log(gamma_true[i])+betad_true%*%zd[i,]+tao_true*mu_hos_true[(observed[i,2])])
  death_time[i]<- rweibull(1, shape=shape_para, scale = scale_para[i])
}


#for y_true=1: cured, no death observed (death censored, delta=1) and no recurrent event observed
for (i in 1:n){
  if (y_true[i]==1){
    observed[i,3]<- 0 #recurrent process indicator
    observed[i,4]<-  rbinom(1,1,0.5) #survival process indicator
    observed[i,(pd+pr+6)]<- 0 #recurrent event as 0
    
  }
  
  #for y_true=0: delta randomly generated
  else{
    
    observed[i,4]<- rbinom(1,1,0.5)
  }
  
}

#compute last observed time based on delta_true

for (i in 1:n){
  if (observed[i,4]==0){
    observed[i,5]<- death_time[i]
  }
  else{
    observed[i,5]<- runif(1,min=0,max=death_time[i])
  }
}

#compute recurrent events for people with y=0

for (i in 1:n){
  if (y_true[i]==0){
    
    t0<- 0
    t_max<- observed[i,5]
    lambda<- function(t){
      val<- gamma_true[i]*kesi_true*t^(kesi_true-1)*exp(betar_true%*%zr[i,]+mu_hos_true[observed[i,2]])
      return (val)
    }
    
    t1<- t0
    lambda_star <- max(sapply(seq(t0, t_max,length.out=1000), lambda))
    X <- numeric()
    while(t1 <= t_max){
      e<-rexp(1,lambda_star)
      u <- runif(1)
      t1<-t1+e
      if (t1>t_max) {
        break
      }
      if(u < lambda(t1)/lambda_star) {
        X <- c(X,t1)
      }
    }
    
    if (length(X)>0){
      observed[i,(pd+pr+6)]<- length(X)
      observed[i,(pd+pr+7):(pd+pr+length(X)+6)]<- X
    } 
    
    else {
      observed[i,(pd+pr+6)]<- 0
    }
  }
  
  if (observed[i,(pd+pr+6)]==0){
    observed[i,3]<- 0
  }
  
  else{
    observed[i,3]<- 1
  }
}

###check the simulation
sum(y_true==1)/n
sum(observed[,12]==0 & observed[,4]==1)/n
######################################################################################################################

### compute individual likelihood 
## basic blocks
#log_hazard function for recurrent process
log_h_ijr<- function(t, zr_ij,gamma_ij,kesi,beta_r,mu_j){
  val<- log(kesi)+(kesi-1)*log(t)+log(gamma_ij)+beta_r%*%zr_ij+mu_j
  return (val)
}

#log_survival function for recurrent process
log_s_ijr<- function(t, zr_ij,gamma_ij,kesi,beta_r,mu_j){
  lambda_ijr<- gamma_ij*exp(c(beta_r%*%zr_ij)+mu_j)
  val<- (-lambda_ijr*(t^kesi))
  return(val)
}

#survival function for recurrent process
s_ijr<- function(t, zr_ij,gamma_ij,kesi,beta_r,mu_j){
  lambda_ijr<- gamma_ij*exp(c(beta_r%*%zr_ij)+mu_j)
  val<- exp(-lambda_ijr*t^kesi)
  return(val)
}
#log_hazard function for survival process
log_h_ijd<- function(t,zd_ij,sigma, gamma_ij,zeta,mu_intercept,beta_d,tao,mu_j){
  val<- log(sigma)+ ((sigma)-1)*log(t)- (zeta*sigma)*log(gamma_ij)- ((mu_intercept+beta_d%*%zd_ij+tao*mu_j)*sigma)
  return(val)
}

#log_survival function for survival process
log_s_ijd<- function(t,zd_ij,sigma, gamma_ij,zeta,mu_intercept,beta_d,tao,mu_j){
  lambda_ijd<- (gamma_ij)^(-zeta*sigma)*exp(-(mu_intercept+c(beta_d%*%zd_ij)+tao*mu_j)*sigma)
  val<- (-lambda_ijd*(t^(sigma)))
  return(val)
}

#survival function for survival process
s_ijd<- function(t,zd_ij,sigma, gamma_ij,zeta,mu_intercept,beta_d,tao,mu_j){
  lambda_ijd<- (gamma_ij)^(-zeta*sigma)*exp(-(mu_intercept+c(beta_d%*%zd_ij)+tao*mu_j)*sigma)
  val<- exp(-lambda_ijd*t^(sigma))
  return(val)
}

##log_ likelihood function for one individual based on above blocks
#t_vector: observed time for recurrent event; default is 0 (no recurrent event observed)
log_l_ij<- function(u_ij,delta_ij,d_ij_min,t_vector,y_ij,zd_ij,sigma, gamma_ij,zeta,mu_intercept,beta_d,tao,mu_j,zr_ij,kesi,beta_r){
  if (u_ij==1 & delta_ij==1){
    sum_result<- 0
    for (i in 1:length(t_vector)){
      sum_result<- sum_result+log_h_ijr(t_vector[i],zr_ij,gamma_ij,kesi,beta_r,mu_j)
    }
    
    val<- sum_result+log_s_ijr(d_ij_min, zr_ij,gamma_ij,kesi,beta_r,mu_j)+log_s_ijd(d_ij_min,zd_ij,sigma, gamma_ij,zeta,mu_intercept,beta_d,tao,mu_j)
  }
  
  else if (u_ij==1 & delta_ij==0){
    sum_result<- 0
    for (i in 1:length(t_vector)){
      sum_result<- sum_result+log_h_ijr(t_vector[i],zr_ij,gamma_ij,kesi,beta_r,mu_j)
    }
    
    val<- sum_result+log_s_ijr(d_ij_min, zr_ij,gamma_ij,kesi,beta_r,mu_j)+log_s_ijd(d_ij_min,zd_ij,sigma, gamma_ij,zeta,mu_intercept,beta_d,tao,mu_j)+log_h_ijd(d_ij_min,zd_ij,sigma, gamma_ij,zeta,mu_intercept,beta_d,tao,mu_j)
  }
  
  else if (u_ij==0 & delta_ij==0){
    val<- (log_s_ijr(d_ij_min, zr_ij,gamma_ij,kesi,beta_r,mu_j)+log_s_ijd(d_ij_min,zd_ij,sigma, gamma_ij,zeta,mu_intercept,beta_d,tao,mu_j)+log_h_ijd(d_ij_min,zd_ij,sigma, gamma_ij,zeta,mu_intercept,beta_d,tao,mu_j))*(1-y_ij)
  }
  
  #no recurrent event + no death
  else if (u_ij==0 & delta_ij==1){
    val<- log(y_ij+ (1-y_ij)*s_ijd(d_ij_min,zd_ij,sigma,gamma_ij,zeta,mu_intercept,beta_d,tao,mu_j)*s_ijr(d_ij_min,zr_ij,gamma_ij,kesi,beta_r,mu_j))
  }
  
  return(val)
}

#global settings 
#format of the input data: there's an input matrix containing all people's information called "observed"
#from the first column to the last are as follows:
#1. id (i from 1 to n), 
#2. hospital (j from 1 to num_hos),
#3. u_ij(recurrent event observed?), (3)
#4. delta_ij(death observed?),       (4)
#5. d_ij_min(lase observed time?),   (5)
#6. zd_ij(covariates related to the survival process), (6: pd+5)
#7. zr_ij(covariates related to the recurrent process), (pd+6:pd+pr+5)
#8. # of recurrent events(0,1,2...), (pd+pr+6)
#9. t_vector(only need to fill in time of happened recurrent events; for example, if no recurrent event happens, then don't need to fill in anything; if 2 recurrent event happens, then fill in 2 time points)
#(pd+pr+7:pd+pr+num_rec_l+6)
#10. z_ij: temporarily assume the same as zd_ij
#need to type by hand
#n<- 1000# num of subjects
#num_hos<- 100 # num of hospitals
#observed<- matrix(0,nrow=n,ncol=50)# input matrix containing all the information 
#pr<- 5
#pd<- 5


#initialization of mcmc parameters
mcmc_samples<- 5000

kesi<- rep(0,times=mcmc_samples)
kesi[1]<- kesi_true
kesi_accept<- 0
metrop_var_kesi<- 0.01
akesi<- 1
bkesi<- 1
sigma2l<- rep(10^2,times=mcmc_samples)


betar<- matrix(0,nrow=mcmc_samples,ncol=pr)
betar[1,]<- betar_true
betar_accept<- rep(0,times=pr)
metrop_var_betar<- rep(0.05,times=pr)

sigma2r<- rep(10^2,times=mcmc_samples)


gamma_sub<- matrix(0,nrow=mcmc_samples,ncol=n)
gamma_sub[1,]<- gamma_true
gamma_sub_accept<- rep(0,times=n)
metrop_var_gamma_sub<- rep(50,times=n)

epislon<- matrix(0,nrow=mcmc_samples,ncol=num_hos)
epislon[1,]<- epislon_true
a0<- 0.1
b0<- 0.1

pc<- num_eta
eta<- matrix(0,nrow=mcmc_samples,ncol=pc)
eta[1,]<-  eta_true
eta_temp<- eta[1,]
eta_accept<- rep(0,times=pc)
metrop_var_eta<- rep(0.8^2,times=pc)
#eta<- matrix(0,nrow=mcmc_samples,ncol=(pc+1))# need to consider intercept
#eta_accept<- rep(0,times=(pc+1))
#metrop_var_eta<- rep(0.05,times=(pc+1))

y_index<- which(observed[,3]==0)# people need to decide between y=0 and y=1
y_not_index<- which (observed[,3]!=0) # people don't need to decide between y=0 and y=1; y always 0 (anyway its value won't affect the likelihood, so what value y takes on doesn't really matter)
p<- matrix(0,nrow=mcmc_samples,ncol=length(y_index)) # only these people need to consider the probability of p
#z_matrix<- cbind(rep(1,times=length(y_index)),observed[y_index,6:(pd+5)])
z_matrix<- cbind(zr[y_index,],zd[y_index,3])
logit_prob_initial<- z_matrix%*%eta[1,]
p[1,] <- 1/(1 + exp(-logit_prob_initial))

y<- matrix(0,nrow=mcmc_samples,ncol=n)
y[1,y_index]<- rbinom(length(y_index),1,p[1,])
y_temp_decision<- rep(0,times=n)

#survival process parameters
#sigma<- rep(0,times=mcmc_samples)
#sigma[1]<- sigma_true
#sigma_accept<- 0
#metrop_var_sigma<- 0.01
#asigma<- 0.1
#bsigma<- 0.1

zeta<- rep(0,times=mcmc_samples)
zeta[1]<- zeta_true
zeta_accept<- 0 
metrop_var_zeta<- 0.01
sigma2zeta<- 10^2

mu_intercept<- rep(0,times=mcmc_samples)
mu_intercept[1]<- mu_intercept_true
mu_intercept_accept<- 0 
metrop_var_mu_intercept<- 0.01
sigma2mu_intercept<- 10^2

betad<- matrix(0,nrow=mcmc_samples,ncol=pd)
betad[1,]<-  betad_true
betad_accept<- rep(0,times=pd)
metrop_var_betad<- rep(0.05,times=pd)

sigma2d<- rep(10^2,times=mcmc_samples)

tao<- rep(0,times=mcmc_samples)
tao[1]<- tao_true
tao_accept<- 0
metrop_var_tao<- 0.01
sigma2tao<- 10^2

mu_hos<- matrix(0,nrow=mcmc_samples,ncol=num_hos)
sigma<-matrix(0,nrow=mcmc_samples,ncol=n)
num_dp_cluster<- 10# need to think about a value? currently assign each hospital a cluster itself
num_dp_cluster_n<-5
mu_hos[1,]<- mu_hos_true
sigma[1,]<-sigma_true

s_label<- array(0,dim=c(mcmc_samples,num_hos,num_dp_cluster))
s_label_n<- array(0,dim=c(mcmc_samples,n,num_dp_cluster_n))
s_label[1,,]<- t(rmultinom(num_hos, size = 1, prob = rep((1/num_dp_cluster),times=num_dp_cluster)))
s_label_n[1,,]<- t(rmultinom(n, size = 1, prob = rep((1/num_dp_cluster_n),times=num_dp_cluster_n)))
pai<- rep(0,times=num_dp_cluster)
pai_n<- rep(0,times=num_dp_cluster_n)

b_dp_cluster<- matrix(0,nrow=mcmc_samples,ncol=num_dp_cluster)
b_dp_cluster_n<- matrix(3,nrow=mcmc_samples,ncol=num_dp_cluster_n)
b_dp_cluster_accept<- rep(0,times=num_dp_cluster)
b_dp_cluster_accept_n<- rep(0,times=num_dp_cluster_n)
metrop_var_b_dp_cluster<- rep(0.5,times=num_dp_cluster)
metrop_var_b_dp_cluster_n<- rep(1,times=num_dp_cluster_n)
#metrop_var_b_dp_cluster[c(15,17,18)]<- rep(2,times=3)
#metrop_var_b_dp_cluster[c(14,16,20)]<- rep(0.5,times=3)
#metrop_var_b_dp_cluster[c(1,4,5,11)]<- rep(3,times=4)
#metrop_var_b_dp_cluster[19]<- 0.5


w_dp_cluster<- matrix(0,nrow=mcmc_samples,ncol=num_dp_cluster)
w_dp_cluster_n<- matrix(0,nrow=mcmc_samples,ncol=num_dp_cluster_n)
w_dp_cluster[1,]<- rep((1/num_dp_cluster),times=num_dp_cluster)
w_dp_cluster_n[1,]<- rep((1/num_dp_cluster_n),times=num_dp_cluster_n)
w_inter_dp_cluster<- matrix(0,nrow=mcmc_samples,ncol=num_dp_cluster)
w_inter_dp_cluster_n<- matrix(0,nrow=mcmc_samples,ncol=num_dp_cluster_n)
w_inter_dp_cluster[1,]<- rep((1/num_dp_cluster),times=num_dp_cluster)
w_inter_dp_cluster_n[1,]<- rep((1/num_dp_cluster_n),times=num_dp_cluster_n)
alpha<- rep(1,times=mcmc_samples)

v2<- 1^2

#temp values used for updating, which means for parameters not directly affecting the updating process, temp values are not needed even with metroplis algorithm
betar_temp<- betar[1,]
gamma_temp<- gamma_sub[1,]
y_temp<- y[1,]
kesi_temp<-kesi[1]

sigma_temp<- sigma_true
zeta_temp<- zeta[1]
mu_intercept_temp<- mu_intercept[1]
betad_temp<- betad[1,]
tao_temp<- tao[1]
mu_hos_temp<- mu_hos[1,]


#log_h_functions
log_h_eta<- function(eta_val){
  eta_temp[k]<- eta_val
  logit_prob<- z_matrix%*%eta_temp
  prob<- 1/(1 + exp(-logit_prob))
  val<- sum(dbinom(x = y[(i-1),y_index],
                   size = rep(1, 
                              times = length(y_index)),
                   prob = prob,
                   log = TRUE)) - (eta_val^2)/(2*(10^2))
  return(val)
  
}

log_h_betar<- function(betar_val){
  sloglik<- 0
  betar_temp[k]<- betar_val
  
  for (l in 1:n){
    num_rec_l<- observed[l,(pd+pr+6)]
    
    if (num_rec_l>0){
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],observed[l,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],kesi_temp,betar_temp)
    }
    
    else{
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],0,y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],kesi_temp,betar_temp)
    }
    
  }
  
  val<- sloglik- betar_val^2/(2*sigma2r[(i-1)])
  
  return(val)
}

log_h_kesi<- function(kesi_val){
  sloglik<- 0
  kesi_temp<- kesi_val
  
  for (l in 1:n){
    num_rec_l<- observed[l,(pd+pr+6)]
    
    if (num_rec_l>0){
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],observed[l,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],kesi_temp,betar_temp)
    }
    
    else{
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],0,y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],kesi_temp,betar_temp)
    }
    
  }
  
  val<- sloglik+ (akesi-1)*log(kesi_val)-(kesi_val/bkesi)
  
  return(val)
}

log_h_gamma_sub<- function(gamma_sub_val){
  gamma_temp[k]<- gamma_sub_val
  num_rec_l<- observed[k,(pd+pr+6)]
  val<- log_l_ij(observed[k,3],observed[k,4],observed[k,5],observed[k,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[k],observed[k,6:(pd+5)],sigma_temp[k], gamma_temp[k],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[k,2]],observed[k,(pd+6):(pd+pr+5)],kesi_temp,betar_temp)
  val<- val- ((log(gamma_sub_val))^2)/(2*epislon[(i-1),observed[k,2]]) # observed[k,2] denotes jth hospital of this subject
  
  return(val)
  
}

#log_h_sigma<- function(sigma_val){
 # sloglik<- 0
  #sigma_temp<- sigma_val
  
  #for (l in 1:n){
    #num_rec_l<- observed[l,(pd+pr+6)]
    
    #if (num_rec_l>0){
      #sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],observed[l,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[l],observed[l,6:(pd+5)],sigma_temp, gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
    #}
    
    #else{
      #sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],0,y_temp[l],observed[l,6:(pd+5)],sigma_temp, gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
    #}
    
  #}
  
  #val<- sloglik- (asigma+1)*log(sigma_val)-(sigma_val/bsigma)
  
  #return(val)
#}

log_h_zeta<- function(zeta_val){
  sloglik<- 0
  zeta_temp<- zeta_val
  
  for (l in 1:n){
    num_rec_l<- observed[l,(pd+pr+6)]
    
    if (num_rec_l>0){
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],observed[l,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],kesi_temp,betar_temp)
    }
    
    else{
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],0,y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],kesi_temp,betar_temp)
    }
    
  }
  
  val<- sloglik- (zeta_val)^2/(2*sigma2zeta)
  
  return(val)
}

log_h_mu_intercept<- function(mu_intercept_val){
  sloglik<- 0
  mu_intercept_temp<- mu_intercept_val
  
  for (l in 1:n){
    num_rec_l<- observed[l,(pd+pr+6)]
    
    if (num_rec_l>0){
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],observed[l,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],kesi_temp,betar_temp)
    }
    
    else{
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],0,y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],kesi_temp,betar_temp)
    }
    
  }
  
  val<- sloglik- (mu_intercept_val)^2/(2*sigma2mu_intercept)
  
  return(val)
}

log_h_betad<- function(betad_val){
  sloglik<- 0
  betad_temp[k]<- betad_val
  
  for (l in 1:n){
    num_rec_l<- observed[l,(pd+pr+6)]
    
    if (num_rec_l>0){
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],observed[l,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],kesi_temp,betar_temp)
    }
    
    else{
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],0,y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],kesi_temp,betar_temp)
    }
    
  }
  
  val<- sloglik- betad_val^2/(2*sigma2d[(i-1)])
  
  return(val)
}


log_h_tao<- function(tao_val){
  sloglik<- 0
  tao_temp<- tao_val
  
  for (l in 1:n){
    num_rec_l<- observed[l,(pd+pr+6)]
    
    if (num_rec_l>0){
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],observed[l,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],kesi_temp,betar_temp)
    }
    
    else{
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],0,y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],kesi_temp,betar_temp)
    }
    
  }
  
  val<- sloglik- (tao_val)^2/(2*sigma2tao)
  
  return(val)
}

log_h_b_dp_cluster<- function(b_dp_cluster_val){
  sloglik<- 0
  mu_hos_temp[s_label[i,,k]==1]<- b_dp_cluster_val
  
  for (l in 1:(n-1)){
    num_rec_l<- observed[l,(pd+pr+6)]
    
    if (num_rec_l>0){
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],observed[l,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],kesi_temp,betar_temp)
    } else{
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],0,y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],kesi_temp,betar_temp)
    }
    
    if (is.nan(sloglik)==T){
      print(l)
      break
    }
    
  }
  
  val<- sloglik- b_dp_cluster_val^2/(2*v2)
  
  return(val)
}

log_h_b_dp_cluster_n<- function(b_dp_cluster_val_n){
  sloglik<- 0
  sigma_temp[s_label_n[i,,k]==1]<- b_dp_cluster_val_n
  
  for (l in 1:(n-1)){
    num_rec_l<- observed[l,(pd+pr+6)]
    
    if (num_rec_l>0){
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],observed[l,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],kesi_temp,betar_temp)
    } else{
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],0,y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],kesi_temp,betar_temp)
    }
    
    #if (is.nan(sloglik)==T){
    #  print(l)
    #  break
    #}
    
  }
  val<- sloglik- b_dp_cluster_val_n^2/(2*v2)
  return(val)
}
###MCMC algorithm

for (i in 2:mcmc_samples){
  # betar
  betar[i,]<- betar[(i-1),]
  
  for (k in 1:pr){
    betar_proposed<- rnorm(n=1,mean=betar[(i-1),k],sd=sqrt(metrop_var_betar[k]))
    
    ratio_1<- exp(log_h_betar(betar_proposed) - log_h_betar(betar[(i-1),k]))
    
    if(ratio_1>=runif(n=1,min=0,max=1)){
      betar[i,k]<- betar_proposed
      betar_accept[k]<- betar_accept[k]+1
    }
    
    betar_temp[k]<- betar[i,k] #make sure betar_temp always uses the latest updated value
  }
  
  
  kesi[i]<- kesi[(i-1)]
  
  kesi_proposed<- rnorm(n=1,mean=kesi[(i-1)],sd=sqrt(metrop_var_kesi))
  
  if (kesi_proposed>0) {
    ratio<- exp(log_h_kesi(kesi_proposed) - log_h_kesi(kesi[(i-1)]))
    
    if(ratio>=runif(n=1,min=0,max=1)){
      kesi[i]<- kesi_proposed
      kesi_accept<- kesi_accept+1
    }
  }
  
  kesi_temp<- kesi[i] 
  
  # gamma_sub
  gamma_sub[i,]<- gamma_sub[(i-1),]
  
  for (k in 1:n){
    gamma_sub_proposed<- rnorm(n=1,mean=gamma_sub[(i-1),k],sd=sqrt(metrop_var_gamma_sub[k]))
    
    if (gamma_sub_proposed>0){
      
      ratio<- exp(log_h_gamma_sub(gamma_sub_proposed) - log_h_gamma_sub(gamma_sub[(i-1),k]))
      
      if(ratio>=runif(n=1,min=0,max=1)){
        gamma_sub[i,k]<- gamma_sub_proposed
        gamma_sub_accept[k]<- gamma_sub_accept[k]+1
      }
      
    }
    
    gamma_temp[k]<- gamma_sub[i,k] 
    
    if ((gamma_sub_accept[k]/i)>0.5){
      metrop_var_gamma_sub[k]<- metrop_var_gamma_sub[k]+5
    } 
    
    if ((gamma_sub_accept[k]/i)<0.2) {
      metrop_var_gamma_sub[k] <- metrop_var_gamma_sub[k]-10
      if (metrop_var_gamma_sub[k] <= 0){
        metrop_var_gamma_sub[k]<- 0.00000001
      }
    }
  }
  
  #epislon
  for (k in 1:num_hos){
    hos_ind_chosen<- which(observed[,2]==k)
    epislon[i,k]<- 1/rgamma(1, (a0+ (sum(observed[,2]==k)/2)), (b0+ ((sum((log(gamma_sub[i,hos_ind_chosen]))^2))/2)))
  }
  
  # p 
  eta[i,]<- eta[(i-1),]
  
  for (k in 1:pc){
    eta_proposed<- rnorm(n=1,mean=eta[(i-1),k],sd=sqrt(metrop_var_eta[k]))
    
    ratio<- exp(log_h_eta(eta_proposed) - log_h_eta(eta[(i-1),k]))
    
    if(ratio>=runif(n=1,min=0,max=1)){
      eta[i,k]<- eta_proposed
      eta_accept[k]<- eta_accept[k]+1
    }
    
    eta_temp[k]<- eta[i,k] #make sure always uses the latest updated value
  }
  
  p[i,]<- 1/(1 + exp(-z_matrix%*%eta[i,]))
  
  # y: only people with delta=1 and num_rec=0 have this 
  correct_indicator_num<- 0
  
  for (k in y_index){
    #loglik for 1
    y_temp[k]<- 1
    val_1<- log_l_ij(observed[k,3],observed[k,4],observed[k,5],0,y_temp[k],observed[k,6:(pd+5)],sigma_temp[k], gamma_temp[k],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[k,2]],observed[k,(pd+6):(pd+pr+5)],kesi_temp,betar_temp)
    p_index<- which(y_index==k)
    val_1<- val_1+log (p[i,p_index]) # the kth person's corresponding probability is its index in y_index
    
    #loglik for 0
    y_temp[k]<- 0
    val_0<- log_l_ij(observed[k,3],observed[k,4],observed[k,5],0,y_temp[k],observed[k,6:(pd+5)],sigma_temp[k], gamma_temp[k],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[k,2]],observed[k,(pd+6):(pd+pr+5)],kesi_temp,betar_temp)
    val_0<- val_0+log (1-p[i,p_index])
    
    max_num<- max(val_1,val_0)
    val_1<- val_1- max_num
    val_0<- val_0- max_num
    
    val_1<- exp(val_1)/sum(exp(val_1)+exp(val_0))
    y[i,k]<- rbinom(1,1,val_1)
    y_temp[k]<-  y[i,k]
    
    #compute result for gamma[k]
    if (mean(y[1:i,k])>0.5){
      y_temp_decision[k]<- 1
    } else {
      y_temp_decision[k]<- 0
    }
    
    if (y_temp_decision[k]==y_true[k]){
      correct_indicator_num<- correct_indicator_num+1
    }
    
    
  }
  
  # sigma
  #sigma[i]<- sigma[(i-1)]
  
  #sigma_proposed<- rnorm(n=1,mean=sigma[(i-1)],sd=sqrt(metrop_var_sigma))
  
  #if (sigma_proposed>0) {
    #ratio<- exp(log_h_sigma(sigma_proposed) - log_h_sigma(sigma[(i-1)]))
    
    #if(ratio>=runif(n=1,min=0,max=1)){
      #sigma[i]<- sigma_proposed
      #sigma_accept<- sigma_accept+1
    #}
 # }
  
 # sigma_temp<- sigma[i] 
  
  for (l in 1:n){
    for (s in 1:num_dp_cluster_n){
      num_rec_l<- observed[l,(pd+pr+6)]
      sigma_temp[l]<- b_dp_cluster_n[(i-1),s]
      pai_n[s]<- log(w_dp_cluster_n[(i-1),s])+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],observed[l,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],kesi_temp,betar_temp)
    }
    pai_n<- pai_n- max(pai_n)
    pai_n<- exp(pai_n)/sum(exp(pai_n))
    s_label_n[i,l,]<- t(rmultinom(1, size = 1, prob = pai_n))
    
    # need to change corresponding b_dp_cluster
    sigma_temp[l]<- b_dp_cluster_n[(i-1),(s_label_n[i,l,]==1)]
  }
  
  #b_dp_cluster: updating of this is equivalent update of mu_hos
  b_dp_cluster_n[i,]<- b_dp_cluster_n[(i-1),]
  
  for (k in 1:num_dp_cluster_n){
    b_dp_propose_n<- abs(rnorm(n=1,mean= b_dp_cluster_n[(i-1),k],sd=0.01))
    ratio_2<- exp(log_h_b_dp_cluster_n(b_dp_propose_n) - log_h_b_dp_cluster_n(b_dp_cluster_n[(i-1),k]))
      if(ratio_2>=runif(n=1,min=0,max=1)){
        b_dp_cluster_n[i,k]<- b_dp_propose_n
        b_dp_cluster_accept_n[k]<- b_dp_cluster_accept_n[k]+1
      }
    sigma[i,(s_label_n[i,,k]==1)]<- b_dp_cluster_n[i,k] 
    sigma_temp[(s_label_n[i,,k]==1)]<- b_dp_cluster_n[i,k]
    
  }
  
  #weights
  w_inter_dp_cluster_n[i,num_dp_cluster_n]<- 1
  for (k in 1:(num_dp_cluster_n-1)){
    s_larger_cluster<- 0
    for (s in (k+1):num_dp_cluster_n){
      s_larger_cluster<- s_larger_cluster+ sum(s_label_n[i,,s]==1)
    }
    w_inter_dp_cluster_n[i,k]<- rbeta(1, (1+sum(s_label_n[i,,k]==1)),alpha[(i-1)]+ s_larger_cluster)
  }
  
  w_dp_cluster_n[i,1]<- w_inter_dp_cluster_n[i,1]
  
  for (k in 2:num_dp_cluster_n){
    w_dp_cluster_n[i,k]<- w_inter_dp_cluster_n[i,k]
    for (s in 1:(k-1)){
      w_dp_cluster_n[i,k]<- w_dp_cluster_n[i,k]*(1-w_inter_dp_cluster_n[i,s])
    }
  }
  
  
  
  
  #zeta
  
  zeta[i]<- zeta[(i-1)]
  
  zeta_proposed<- rnorm(n=1,mean=zeta[(i-1)],sd=sqrt(metrop_var_zeta))
  
  ratio<- exp(log_h_zeta(zeta_proposed) - log_h_zeta(zeta[(i-1)]))
  
  if(ratio>=runif(n=1,min=0,max=1)){
    zeta[i]<- zeta_proposed
    zeta_accept<- zeta_accept+1
  }
  
  zeta_temp<- zeta[i] 
  
  #mu_intercept
  
  mu_intercept[i]<- mu_intercept[(i-1)]
  
  mu_intercept_proposed<- rnorm(n=1,mean= mu_intercept[(i-1)],sd=sqrt(metrop_var_mu_intercept))
  
  ratio<- exp(log_h_mu_intercept(mu_intercept_proposed) - log_h_mu_intercept(mu_intercept[(i-1)]))
  
  if(ratio>=runif(n=1,min=0,max=1)){
    mu_intercept[i]<- mu_intercept_proposed
    mu_intercept_accept<-  mu_intercept_accept+1
  }
  
  mu_intercept_temp<-  mu_intercept[i]
  
  # betad
  betad[i,]<- betad[(i-1),]
  
  for (k in 1:pd){
    betad_proposed<- rnorm(n=1,mean=betad[(i-1),k],sd=sqrt(metrop_var_betad[k]))
    
    ratio<- exp(log_h_betad(betad_proposed) - log_h_betad(betad[(i-1),k]))
    
    if(ratio>=runif(n=1,min=0,max=1)){
      betad[i,k]<- betad_proposed
      betad_accept[k]<- betad_accept[k]+1
    }
    
    betad_temp[k]<- betad[i,k] 
  }
  
  
  #sigma2d
  #sigma2d[i]<- 1/rgamma(1,(asigmad+(pd/2)),(bsigmad+ (sum(betad[i,]^2)/2)))
  
  #tao
  
  tao[i]<- tao[(i-1)]
  
  tao_proposed<- rnorm(n=1,mean= tao[(i-1)],sd=sqrt(metrop_var_tao))
  
  ratio<- exp(log_h_tao(tao_proposed) - log_h_tao(tao[(i-1)]))
  
  if(ratio>=runif(n=1,min=0,max=1)){
    tao[i]<- tao_proposed
    tao_accept<- tao_accept+1
  }
  
  tao_temp<- tao[i] 
  
  #mu_hos
  #s_label
  
  for (k in 1:num_hos){
    for (s in 1:num_dp_cluster){
      mu_hos_temp[k]<- b_dp_cluster[(i-1),s]
      ploglik<- 0
      for (l in 1:n){ #all individuals in this hospital
        if (observed[l,2]==k){
          num_rec_l<- observed[l,(pd+pr+6)]
          ploglik<- ploglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],observed[l,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],kesi_temp,betar_temp)
        }
        
      }
      
      pai[s]<- s_label[(i-1),k,s]*log(w_dp_cluster[(i-1),s])+ ploglik
    }
    
    pai<- pai- max(pai)
    pai<- exp(pai)/sum(exp(pai))
    s_label[i,k,]<- t(rmultinom(1, size = 1, prob = pai))
    
    # need to change corresponding b_dp_cluster
    mu_hos_temp[k]<- b_dp_cluster[(i-1),(s_label[i,k,]==1)]
  }
  
  #b_dp_cluster: updating of this is equivalent update of mu_hos
  b_dp_cluster[i,]<- b_dp_cluster[(i-1),]
  
  for (k in 1:num_dp_cluster){
    b_dp_proposed<- rnorm(n=1,mean= b_dp_cluster[(i-1),k],sd=sqrt(metrop_var_b_dp_cluster[k]))
    
    ratio<- exp(log_h_b_dp_cluster(b_dp_proposed) - log_h_b_dp_cluster(b_dp_cluster[(i-1),k]))
    
    if(ratio>=runif(n=1,min=0,max=1)){
      b_dp_cluster[i,k]<- b_dp_proposed
      b_dp_cluster_accept[k]<- b_dp_cluster_accept[k]+1
    }
    
    mu_hos[i,(s_label[i,,k]==1)]<- b_dp_cluster[i,k] 
    mu_hos_temp[(s_label[i,,k]==1)]<- b_dp_cluster[i,k]
    
  }
  
  #weights
  w_inter_dp_cluster[i,num_dp_cluster]<- 1
  for (k in 1:(num_dp_cluster-1)){
    s_larger_cluster<- 0
    for (s in (k+1):num_dp_cluster){
      s_larger_cluster<- s_larger_cluster+ sum(s_label[i,,s]==1)
    }
    w_inter_dp_cluster[i,k]<- rbeta(1, (1+sum(s_label[i,,k]==1)),alpha[(i-1)]+ s_larger_cluster)
  }
  
  w_dp_cluster[i,1]<- w_inter_dp_cluster[i,1]
  
  for (k in 2:num_dp_cluster){
    w_dp_cluster[i,k]<- w_inter_dp_cluster[i,k]
    for (s in 1:(k-1)){
      w_dp_cluster[i,k]<- w_dp_cluster[i,k]*(1-w_inter_dp_cluster[i,s])
    }
  }
  
  print(c("Completion %:", 
          round(100*i/mcmc_samples, 2)))
  print(c("prop_of_correct_y:",correct_indicator_num/length(y_index)))
  print(c("betad Acceptance %:", 
          round(100*min(betad_accept/i), 2), 
          round(100*max(betad_accept/i), 2)))
  print(c("betad_true_1",mean(betad[1:i,1])))
  print(c("betad_true_2",mean(betad[1:i,2])))
  print(c("betad_true_3",mean(betad[1:i,3])))
  print(c("betar Acceptance %:", 
          round(100*min(betar_accept/i), 2), 
          round(100*max(betar_accept/i), 2))) 
  print(c("betar_true_1",mean(betar[1:i,1])))
  print(c("betar_true_2",mean(betar[1:i,2])))
  print(c("betar_true_3",mean(betar[1:i,3])))
  print(c("eta Acceptance %:", 
          round(100*min(eta_accept/i), 2), 
          round(100*max(eta_accept/i), 2))) 
  print(c("eta_true_1",mean(eta[1:i,1])))
  print(c("eta_true_2",mean(eta[1:i,2])))
  print(c("eta_true_3",mean(eta[1:i,3])))
  print(c("eta_true_4",mean(eta[1:i,4])))
  #print(c("sigma Acceptance %:",round(100*sigma_accept/i,2)))
  #print(c("sigma_est",mean(sigma[1:i])))
  print(c("b_dp_cluster Acceptance %:", 
          round(100*min(b_dp_cluster_accept/i), 2), 
          round(100*max(b_dp_cluster_accept/i), 2))) 
  print(c("gamma_sub Acceptance %:", 
          round(100*min(gamma_sub_accept/i), 2), 
          round(100*max(gamma_sub_accept/i), 2))) 
  print(c("zeta Acceptance %:", 100*(zeta_accept/i)))
  print(c("zeta_est",mean(zeta[1:i])))
  print(c("mu_intercept Acceptance %:", 100*(mu_intercept_accept/i)))
  print(c("mu_intercept_est",mean(mu_intercept[1:i])))
  print(c("tao Acceptance %:", 100*(tao_accept/i)))
  print(c("tao_est",mean(tao[1:i])))
  
}

save(mcmc_samples,betad, betad_true,betar, betar_true,sigma,sigma_true,y,y_true,y_index,eta,eta_true,tao,zeta,kesi,kesi_true,file = paste0("n_combn_sim_06_all_",n,"_", seed,"_partly_overlap_strong_all_logistic_larger_variance_t11", ".rdata")) 
