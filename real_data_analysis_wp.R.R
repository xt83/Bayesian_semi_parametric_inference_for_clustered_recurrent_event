# grace: setwd("/gpfs/loomis/home.grace/jc3625/project/competing_risk_simulation")
# farnam:

setwd("/gpfs/ysm/project/li_fan/xt83")
library(LaplacesDemon)

stride_data<- read.table("data.txt",sep="\t",header=T)

###distribution of fall_num: this is the fall_numth fall of this participant
#hist(stride_data$fallnum)
#sum(is.na(stride_data$fallnum))/dim(stride_data)[1]# 72% percent of participants don't have fall events

###keep only complete cases
#index_complete<-  complete.cases(stride_data)
#stride_data<-stride_data[index_complete,]

###summarize the distribution of the data
#for each event type: 1 in 2 in 3: in their original paper, they main interest lies in event type 1, adjudicated serious-fall
table(stride_data$endpoint)

###focus on people with adjudicated serious-fall
stride_data_type_1<- stride_data[stride_data$endpoint==1,]
stride_data_type_2<- stride_data[stride_data$endpoint==2,]
stride_data_type_3<- stride_data[stride_data$endpoint==3,]

#remove observations without race status
index_with_race<- which((stride_data_type_1$demo_race_allmiss==1)==FALSE)
stride_data_type_1<- stride_data_type_1[index_with_race,]

#transform the original dataset from long format to wide format in order to summarize the distribution information
stride_data_type_1$fallnum_new<- stride_data_type_1$fallnum
stride_data_type_1[is.na(stride_data_type_1$fallnum_new)==T,]$fallnum_new<- 0

#Add occasion index for each subject
stride_data_type_1$time_new <- 1
for (i in 2:nrow(stride_data_type_1)) {
  if (stride_data_type_1$screenid[i] == stride_data_type_1$screenid[i-1]) {stride_data_type_1$time_new[i] <- stride_data_type_1$time_new[i-1] + 1}
  else {stride_data_type_1$time_new[i]<- 1}
}

check_time_new<- stride_data_type_1[,c(3,115)]

stride_data_type_1_wide<- reshape(stride_data_type_1, timevar="time_new", idvar="screenid", v.names=c("fallnum_new"), direction="wide")
stride_data_type_1_wide$fall_summary<- 0

###summarize distribution using $ before transforming the data.frame into matrix
#distribution of the practice center
#barplot(table(stride_data_type_1_wide$practiceid))# this generates the distribution of participants in different practice centers: generally uniform
#class(table(stride_data_type_1_wide$practiceid))

#distribution of num_of_chronic_condition
#mean(stride_data_type_1_wide$numchron)
#sqrt(var(stride_data_type_1_wide$numchron))

#distribution of age
#mean(stride_data_type_1_wide$age)
#sqrt(var(stride_data_type_1_wide$age))

#distribution of sex
#sum(stride_data_type_1_wide$demo_sex_male==1)
#sum(stride_data_type_1_wide$demo_sex_female==1)

#distribution of race
#sum(stride_data_type_1_wide$demo_race_amind==1)
#sum(stride_data_type_1_wide$demo_race_asian==1)
#sum(stride_data_type_1_wide$demo_race_black==1)
#sum(stride_data_type_1_wide$demo_race_pacis==1)
#sum(stride_data_type_1_wide$demo_race_white==1)
#sum(stride_data_type_1_wide$demo_race_other==1)
#sum(stride_data_type_1_wide$demo_race_multi==1)

#distribution of surrogate consent
#table(stride_data_type_1_wide$surrogate_consent)

#distribution of fall
#table(stride_data_type_1_wide[,ncol(stride_data_type_1_wide)])
#hist(stride_data_type_1_wide[,ncol(stride_data_type_1_wide)])

#distribution of death
#table(stride_data_type_1_wide$dead)

#summarize # of recurrent events for each individual 
stride_data_type_1_wide<- as.matrix(stride_data_type_1_wide)
#colnames(stride_data_type_1_wide)[114]#check the start column is the fallnum_new.1
#colnames(stride_data_type_1_wide)[(ncol(stride_data_type_1_wide)-1)]#check this is the last fallnum_new
#colnames(stride_data_type_1_wide)[ncol(stride_data_type_1_wide)]

for (i in 1:nrow(stride_data_type_1_wide)){
  stride_data_type_1_wide[i,ncol(stride_data_type_1_wide)]<- max(stride_data_type_1_wide[i,114:(ncol(stride_data_type_1_wide)-1)],na.rm = T)
}

#check fallnum distribution 
hist(stride_data_type_1_wide[,"fall_summary"])
table(stride_data_type_1_wide[,"fall_summary"])
4828/6082# 80% 0 falls


#only keep necessary columns
stride_data_type_1_wide<- stride_data_type_1_wide[,c(2,3,5,6,7,8,10,11,14,92:102,ncol(stride_data_type_1_wide))]
colnames(stride_data_type_1_wide)#check 

#summarize time of recurrent falls
stride_data_type_1_recurrent_time<- stride_data_type_1[,c(3,9,ncol(stride_data_type_1))]
stride_data_type_1_recurrent_time_wide<- reshape(stride_data_type_1_recurrent_time, timevar="time_new", idvar="screenid", v.names=c("time"), direction="wide")

#merge
stride_data_type_1_wide_final<- merge(stride_data_type_1_wide, stride_data_type_1_recurrent_time_wide, by = "screenid")
stride_data_type_1_wide_final<- as.matrix(stride_data_type_1_wide_final)

#assign NA to no recurrent event time
stride_data_type_1_wide_final[stride_data_type_1_wide_final[,21]==0,22]<- NA
for (i in 1:nrow(stride_data_type_1_wide_final)){
  if (stride_data_type_1_wide_final[i,21]!=0){
    stride_data_type_1_wide_final[i,22+stride_data_type_1_wide_final[i,21]]<- NA
  }
}

###final dataset waiting to be used
##sex
stride_data_type_1_wide_final<- cbind(stride_data_type_1_wide_final,rep(1,times=nrow(stride_data_type_1_wide_final)))
dim(stride_data_type_1_wide_final)
colnames(stride_data_type_1_wide_final)
stride_data_type_1_wide_final[stride_data_type_1_wide_final[,12]==1,ncol(stride_data_type_1_wide_final)]<- 0
colnames(stride_data_type_1_wide_final)[ncol(stride_data_type_1_wide_final)]<- "sex"
stride_data_type_1_wide_final<- stride_data_type_1_wide_final[,-c(11,12)]
stride_data_type_1_wide_final<- stride_data_type_1_wide_final[,-11]

##race
stride_data_type_1_wide_final<- cbind(stride_data_type_1_wide_final,rep(1,times=nrow(stride_data_type_1_wide_final)),rep(1,times=nrow(stride_data_type_1_wide_final)))
dim(stride_data_type_1_wide_final)[2]

# for black people: 0 0 
colnames(stride_data_type_1_wide_final)[13]
stride_data_type_1_wide_final[stride_data_type_1_wide_final[,13]==1,dim(stride_data_type_1_wide_final)[2]]<- 0
stride_data_type_1_wide_final[stride_data_type_1_wide_final[,13]==1,(dim(stride_data_type_1_wide_final)[2]-1)]<- 0

#for white people: 1 0 
colnames(stride_data_type_1_wide_final)[15]
stride_data_type_1_wide_final[stride_data_type_1_wide_final[,15]==1,dim(stride_data_type_1_wide_final)[2]]<- 0

# for other people: 0 1  
stride_data_type_1_wide_final[(stride_data_type_1_wide_final[,13]==0)&(stride_data_type_1_wide_final[,15]==0) ,(dim(stride_data_type_1_wide_final)[2]-1)]<- 0
colnames(stride_data_type_1_wide_final)[(dim(stride_data_type_1_wide_final)[2]-1)]<- "race_white"
colnames(stride_data_type_1_wide_final)[(dim(stride_data_type_1_wide_final)[2])]<- "race_others"

#check race
stride_data_type_1_wide_final<- stride_data_type_1_wide_final[,-c(11,12,14,16,17)]
sum(stride_data_type_1_wide_final[,ncol(stride_data_type_1_wide_final)]==1) #161 others
sum(stride_data_type_1_wide_final[,(ncol(stride_data_type_1_wide_final)-1)]==1) # 4965 white 
colnames(stride_data_type_1_wide_final)
stride_data_type_1_wide_final<- stride_data_type_1_wide_final[,-c(11,12)]
#colnames(stride_data_type_1_wide_final)

##tertile
stride_data_type_1_wide_final<- cbind(stride_data_type_1_wide_final,rep(1,times=nrow(stride_data_type_1_wide_final)),rep(1,times=nrow(stride_data_type_1_wide_final)))
dim(stride_data_type_1_wide_final)[2]
colnames(stride_data_type_1_wide_final)

#size 1: 0 0 
stride_data_type_1_wide_final[stride_data_type_1_wide_final[,6]==1,ncol(stride_data_type_1_wide_final)]<- 0
stride_data_type_1_wide_final[stride_data_type_1_wide_final[,6]==1,(ncol(stride_data_type_1_wide_final)-1)]<- 0

#size 2: 1 0
stride_data_type_1_wide_final[stride_data_type_1_wide_final[,6]==2,ncol(stride_data_type_1_wide_final)]<- 0

#size 3: 0 1
stride_data_type_1_wide_final[stride_data_type_1_wide_final[,6]==3,(ncol(stride_data_type_1_wide_final)-1)]<- 0

colnames(stride_data_type_1_wide_final)[ncol(stride_data_type_1_wide_final)-1]<- "size_2"
colnames(stride_data_type_1_wide_final)[ncol(stride_data_type_1_wide_final)]<- "size_3"

#id
stride_data_type_1_wide_final[,1]<- 1:nrow(stride_data_type_1_wide_final)
colnames(stride_data_type_1_wide_final)[1]<- "id"

#u_ij
stride_data_type_1_wide_final<- cbind(stride_data_type_1_wide_final,rep(1,times=nrow(stride_data_type_1_wide_final)))
stride_data_type_1_wide_final[stride_data_type_1_wide_final[,"fall_summary"]==0,ncol(stride_data_type_1_wide_final)]<- 0
colnames(stride_data_type_1_wide_final)[ncol(stride_data_type_1_wide_final)]<- "recurrent_event_indicator"

#delta_ij
stride_data_type_1_wide_final[,"dead"]<- rep(1,times=nrow(stride_data_type_1_wide_final))-stride_data_type_1_wide_final[,"dead"]
sum(stride_data_type_1_wide_final[,"dead"]==0)
colnames(stride_data_type_1_wide_final)[8]<- "death_censoring_indicator"

#reorder as requested in the algorithm 
stride_data_type_1_wide_final<- stride_data_type_1_wide_final[,-6]
stride_data_type_1_wide_final_survival_cov<- stride_data_type_1_wide_final[,c("id","practiceid","recurrent_event_indicator","death_censoring_indicator","endtime","intv","urban","majwhite","size_2","size_3","numchron","age","sex","race_white","race_others")]
stride_data_type_1_wide_final_recurrent_cov<- stride_data_type_1_wide_final[,c("intv","urban","majwhite","size_2","size_3","numchron","age","sex","race_white","race_others","fall_summary")]
colnames(stride_data_type_1_wide_final)
stride_data_type_1_wide_final_recurrent_time<- stride_data_type_1_wide_final[,11:15]
stride_data_type_1_wide_final_combine<- cbind(stride_data_type_1_wide_final_survival_cov,stride_data_type_1_wide_final_recurrent_cov,stride_data_type_1_wide_final_recurrent_time)

#check the combined dataset
dim(stride_data_type_1_wide_final_combine)
colnames(stride_data_type_1_wide_final_combine)[15]#ensure delete 'others'
colnames(stride_data_type_1_wide_final_combine)[25]#ensure delete 'others'
stride_data_type_1_wide_final_combine<- stride_data_type_1_wide_final_combine[,-c(15,25)]

#check the final dataset
dim(stride_data_type_1_wide_final_combine)
colnames(stride_data_type_1_wide_final_combine)
sum(stride_data_type_1_wide_final_combine[,"race_white"]==1)
sum(stride_data_type_1_wide_final_combine[,"race_white"]==0)

#remove two variables: urban and majorwhite
stride_data_type_1_wide_final_combine<- stride_data_type_1_wide_final_combine[,-c(7,8,16,17)]
######################################################################################################################

### compute individual likelihood 
## basic blocks
#log_hazard function for recurrent process
log_h_ijr<- function(t,zr_ij,gamma_ij,beta_r,mu_j,lambda_n,time_q){
  tr<-0
  for(q in 1:5){
    tr=tr+lambda_n[q]*ifelse(t>time_q[q]&t<=time_q[q+1],1,0)
  }
  val<- log(tr)+log(gamma_ij)+beta_r%*%zr_ij+mu_j
  return (val)
}

#log_survival function for recurrent process
log_s_ijr<- function(t, zr_ij,gamma_ij,beta_r,mu_j,lambda_n,time_q){
  lambda_ijr<- gamma_ij*exp(c(beta_r%*%zr_ij)+mu_j)
  sr<-0
  for(q in 1:5){
    sr=sr+lambda_n[q]*max(0,min(time_q[q+1]-time_q[q],t-time_q[q]))
  }
  val<-(-sr*lambda_ijr)
  return(val)
}

#survival function for recurrent process
s_ijr<- function(t, zr_ij,gamma_ij,beta_r,mu_j,lambda_n,time_q){
  lambda_ijr<- gamma_ij*exp(c(beta_r%*%zr_ij)+mu_j)
  sr<-0
  for(q in 1:5){
    sr=sr+lambda_n[q]*max(0,min(time_q[q+1]-time_q[q],t-time_q[q]))
  }
  val<- exp(-lambda_ijr*sr)
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
log_l_ij<- function(u_ij,delta_ij,d_ij_min,t_vector,y_ij,zd_ij,sigma, gamma_ij,zeta,mu_intercept,beta_d,tao,mu_j,zr_ij,beta_r,lambda_n,time_q){
  if (u_ij==1 & delta_ij==1){
    sum_result<- 0
    for (i in 1:length(t_vector)){
      sum_result<- sum_result+log_h_ijr(t_vector[i],zr_ij,gamma_ij,beta_r,mu_j,lambda_n,time_q)
    }
    
    val<- sum_result+log_s_ijr(d_ij_min, zr_ij,gamma_ij,beta_r,mu_j,lambda_n,time_q)+log_s_ijd(d_ij_min,zd_ij,sigma, gamma_ij,zeta,mu_intercept,beta_d,tao,mu_j)
  }
  
  else if (u_ij==1 & delta_ij==0){
    sum_result<- 0
    for (i in 1:length(t_vector)){
      sum_result<- sum_result+log_h_ijr(t_vector[i],zr_ij,gamma_ij,beta_r,mu_j,lambda_n,time_q)
    }
    
    val<- sum_result+log_s_ijr(d_ij_min, zr_ij,gamma_ij,beta_r,mu_j,lambda_n,time_q)+log_s_ijd(d_ij_min,zd_ij,sigma, gamma_ij,zeta,mu_intercept,beta_d,tao,mu_j)+log_h_ijd(d_ij_min,zd_ij,sigma, gamma_ij,zeta,mu_intercept,beta_d,tao,mu_j)
  }
  
  else if (u_ij==0 & delta_ij==0){
    val<- (log_s_ijr(d_ij_min, zr_ij,gamma_ij,beta_r,mu_j,lambda_n,time_q)+log_s_ijd(d_ij_min,zd_ij,sigma, gamma_ij,zeta,mu_intercept,beta_d,tao,mu_j)+log_h_ijd(d_ij_min,zd_ij,sigma, gamma_ij,zeta,mu_intercept,beta_d,tao,mu_j))*(1-y_ij)
  }
  
  #no recurrent event + no death
  else if (u_ij==0 & delta_ij==1){
    val<- log(y_ij+ (1-y_ij)*s_ijd(d_ij_min,zd_ij,sigma,gamma_ij,zeta,mu_intercept,beta_d,tao,mu_j)*s_ijr(d_ij_min,zr_ij,gamma_ij,beta_r,mu_j,lambda_n,time_q))
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


n<- nrow(stride_data_type_1_wide_final_combine)
num_hos<- 86
observed<- stride_data_type_1_wide_final_combine
pr=7
pd=7
mcmc_samples<- 10000

t_0<-max(observed[,5])
time_1<-c(0,t_0/5,2*t_0/5,3*t_0/5,4*t_0/5,t_0+1)
time_q<-time_1
lambda_n<- matrix(0,nrow=mcmc_samples,ncol=5)
lambda_n[1,]<- c(1,1,1,1,1)
metrop_var_lambda_n<- rep(0.05,times=5)
sigma2l<- rep(10^2,times=mcmc_samples)
lambda_n_temp<-lambda_n[1,]

kesi<- rep(0,times=mcmc_samples)
kesi[1]<- 1
kesi_accept<- 0
metrop_var_kesi<- 0.01
akesi<- 1
bkesi<- 1
sigma2l<- rep(10^2,times=mcmc_samples)


betar<- matrix(0,nrow=mcmc_samples,ncol=pr)
betar_accept<- rep(0,times=pr)
metrop_var_betar<- rep(0.05,times=pr)
metrop_var_betar[4]<- 0.005
metrop_var_betar[5]<- 0.000005
metrop_var_betar[6]<- 0.1
metrop_var_betar[7]<- 0.1

sigma2r<- rep(10^2,times=mcmc_samples)


gamma_sub<- matrix(0,nrow=mcmc_samples,ncol=n)
gamma_sub[1,]<- rep(1^2,times=n)
gamma_sub_accept<- rep(0,times=n)
metrop_var_gamma_sub<- rep(30,times=n)

epislon<- matrix(0,nrow=mcmc_samples,ncol=num_hos)
epislon[1,]<-rep(1^2,times=num_hos)
a0<- 0.1
b0<- 0.1


#eta<- matrix(0,nrow=mcmc_samples,ncol=(pc+1))# need to consider intercept
#eta_accept<- rep(0,times=(pc+1))
#metrop_var_eta<- rep(0.05,times=(pc+1))

y_index<- which(observed[,3]==0)# people need to decide between y=0 and y=1
y_not_index<- which (observed[,3]!=0) # people don't need to decide between y=0 and y=1; y always 0 (anyway its value won't affect the likelihood, so what value y takes on doesn't really matter)
p<- matrix(0.5,nrow=mcmc_samples,ncol=length(y_index))

y<- matrix(0,nrow=mcmc_samples,ncol=n)
y[1,y_index]<- rbinom(length(y_index),1,p[1,])
y_temp_decision<- rep(0,times=n)
z_matrix<- cbind(rep(1,times=length(y_index)),stride_data_type_1_wide_final_combine[y_index,11:14])
num_eta<- dim(z_matrix)[2]-1
pc<- num_eta
eta<- matrix(0,nrow=mcmc_samples,ncol=(pc+1))
eta_temp<- eta[1,]
eta_accept<- rep(0,times=(pc+1))
metrop_var_eta<- rep(0.8^2,times=(pc+1))
#survival process parameters
#sigma<- rep(0,times=mcmc_samples)
#sigma[1]<- sigma_true
#sigma_accept<- 0
#metrop_var_sigma<- 0.01
#asigma<- 0.1
#bsigma<- 0.1

zeta<- rep(0,times=mcmc_samples)
zeta_accept<- 0 
metrop_var_zeta<- 0.01
sigma2zeta<- 10^2

mu_intercept<- rep(0,times=mcmc_samples)
metrop_var_mu_intercept<- 0.01
sigma2mu_intercept<- 10^2

betad<- matrix(0,nrow=mcmc_samples,ncol=pd)
betad_accept<- rep(0,times=pd)
metrop_var_betad<- rep(0.05,times=pd)
metrop_var_betad[4]<- 0.005
metrop_var_betad[5]<- 0.000005
metrop_var_betad[6]<- 0.1
metrop_var_betad[7]<- 0.1


sigma2d<- rep(10^2,times=mcmc_samples)

tao<- rep(0,times=mcmc_samples)
tao_accept<- 0
metrop_var_tao<- 0.01
sigma2tao<- 10^2

mu_hos<- matrix(0,nrow=mcmc_samples,ncol=num_hos)
sigma<-matrix(1,nrow=mcmc_samples,ncol=n)
num_dp_cluster<- 10# need to think about a value? currently assign each hospital a cluster itself
num_dp_cluster_n<-5

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

sigma_temp<- sigma[1,]
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
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],observed[l,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
    }
    
    else{
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],0,y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
    }
    
  }
  
  val<- sloglik- betar_val^2/(2*sigma2r[(i-1)])
  
  return(val)
}

log_h_lambda_n<- function(lambda_n_val){
  sloglik<- 0
  lambda_n_temp[s]<- lambda_n_val
  
  for (l in 1:n){
    num_rec_l<- observed[l,(pd+pr+6)]
    
    if (num_rec_l>0){
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],observed[l,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
    }
    
    else{
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],0,y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
    }
    
  }
  
  val<- sloglik-lambda_n_val^2/(2*sigma2l[(i-1)])
  
  return(val)
}

log_h_gamma_sub<- function(gamma_sub_val){
  gamma_temp[k]<- gamma_sub_val
  num_rec_l<- observed[k,(pd+pr+6)]
  val<- log_l_ij(observed[k,3],observed[k,4],observed[k,5],observed[k,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[k],observed[k,6:(pd+5)],sigma_temp[k], gamma_temp[k],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[k,2]],observed[k,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
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
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],observed[l,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
    }
    
    else{
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],0,y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
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
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],observed[l,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
    }
    
    else{
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],0,y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
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
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],observed[l,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
    }
    
    else{
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],0,y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
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
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],observed[l,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
    }
    
    else{
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],0,y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
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
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],observed[l,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
    } else{
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],0,y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
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
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],observed[l,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
    } else{
      sloglik<- sloglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],0,y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
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
  
  
  lambda_n[i,]<-lambda_n[(i-1),]
  for(s in 1:5){
    lambda_n_proposed<- abs(rnorm(n=1,mean=lambda_n[(i-1),s],sd=sqrt(metrop_var_lambda_n[s])))
    if (lambda_n_proposed>0) {
    ratio<- exp(log_h_lambda_n(lambda_n_proposed) - log_h_lambda_n(lambda_n[(i-1),s]))
    
    if(ratio>=runif(n=1,min=0,max=1)){
      lambda_n[i,s]<- lambda_n_proposed
    }
  }
  lambda_n_temp[s]<-lambda_n[i,s]
  }
  
  # gamma_sub
   gamma_sub[i,]<- gamma_sub[(i-1),]
  
  for (k in 1:n){
    gamma_sub_proposed<- rnorm(n=1,mean=gamma_sub[(i-1),k],sd=sqrt(metrop_var_gamma_sub[k]))
    
    if (gamma_sub_proposed>0){
      
      ratio<- exp(log_h_gamma_sub(gamma_sub_proposed) - log_h_gamma_sub(gamma_sub[(i-1),k]))
      if(is.nan(ratio)==F){
        if(ratio>=runif(n=1,min=0,max=1)){
          gamma_sub[i,k]<- gamma_sub_proposed
          gamma_sub_accept[k]<- gamma_sub_accept[k]+1
        }
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
  eta[i,]<- eta[(i-1),]
  
  for (k in 1:pc){
    eta_proposed<- rnorm(n=1,mean=eta[(i-1),k],sd=0.0001)
    
    ratio<- exp(log_h_eta(eta_proposed) - log_h_eta(eta[(i-1),k]))
    
    if(ratio>=runif(n=1,min=0,max=1)){
      eta[i,k]<- eta_proposed
    }
    
    eta_temp[k]<- eta[i,k] #make sure always uses the latest updated value
  }
  
  p[i,]<- 1/(1 + exp(-z_matrix%*%eta[i,]))
  
  # y: only people with delta=1 and num_rec=0 have this 
  correct_indicator_num<- 0
  
  for (k in y_index){
    #loglik for 1
    y_temp[k]<- 1
    val_1<- log_l_ij(observed[k,3],observed[k,4],observed[k,5],0,y_temp[k],observed[k,6:(pd+5)],sigma_temp[k], gamma_temp[k],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[k,2]],observed[k,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
    p_index<- which(y_index==k)
    val_1<- val_1+log (p[i,p_index]) # the kth person's corresponding probability is its index in y_index
    
    #loglik for 0
    y_temp[k]<- 0
    val_0<- log_l_ij(observed[k,3],observed[k,4],observed[k,5],0,y_temp[k],observed[k,6:(pd+5)],sigma_temp[k], gamma_temp[k],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[k,2]],observed[k,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
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
      pai_n[s]<- log(w_dp_cluster_n[(i-1),s])+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],observed[l,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
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
          ploglik<- ploglik+ log_l_ij(observed[l,3],observed[l,4],observed[l,5],observed[l,(pd+pr+7):(pd+pr+num_rec_l+6)],y_temp[l],observed[l,6:(pd+5)],sigma_temp[l], gamma_temp[l],zeta_temp,mu_intercept_temp,betad_temp,tao_temp,mu_hos_temp[observed[l,2]],observed[l,(pd+6):(pd+pr+5)],betar_temp,lambda_n_temp,time_q)
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
  if (i%%200==0){
    #save(mcmc_samples,betad,betar,sigma,kesi,y,y_index,file = paste0("stride_data_type_1_",as.numeric(args[1]), ".rdata")) 
    save.image(file="stride_data_type_1_analysis_informative_g0.rdata")
  }
}


