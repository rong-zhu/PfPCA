### Ground truth of f(s) ######
No<-50
s<-seq(-90,90,22.5); ## 9 stimlus
p<-length(s)
sigma0<-20;
TuningCurve<-function(s,sigma=sigma0){
  fs0<-dnorm(s/sigma)
  c<-5/max(fs0);
  fs<-c*fs0+0.5;#0.2,0.5
  return(fs);
}
fs<-TuningCurve(s);plot(s,fs)


################## Simulation setting ############################
######### Example: Case S4: sharpening ###########################
###  the setting of other three cases are attached at the end ####
### from 9 stimulus in [0,180]
phi4<-log(TuningCurve(s,sigma=1.2*sigma0))-log(TuningCurve(s,sigma=0.8*sigma0));
phi4n<-phi4/(sum(phi4^2))^0.5
aa<-exp(phi4+log(fs))
c_S4<-sum((log(aa)-log(fs))^2)

##### from 181 degrees from 0 to 180 ###
##### for plotting the groundtruth #####
s180<--90:90;fs180<-TuningCurve(s180)
###
phi180_4<-log(TuningCurve(s180,sigma=1.01*sigma0))-log(TuningCurve(s180,sigma=1*sigma0));
phi180_4n<-phi180_4/(sum(phi180_4^2))^0.5;
phi180_4M<-log(TuningCurve(s180,sigma=1.2*sigma0))-log(TuningCurve(s180,sigma=0.8*sigma0));
sigma_S4<-sum(phi180_4n*phi180_4M)


##################################################
## Simulate the spike counts and run Pf-PCA ###### 
##################################################
phi<-phi4n; 
sig2<-c_S4; 
rho<-0.8;
No<-50;p<-length(s); 
ra<-(1/rho-1)*sig2/(p-1) # control the noise
Sigma_mu<-sig2*phi%*%t(phi)+diag(rep(ra,p));

xi_bI<-matrix(rnorm(No*p,0,1),No,p);a.eig <- eigen(Sigma_mu);
Sig_L <- a.eig$vectors %*% diag(sqrt(a.eig$values));
xi_b<-xi_bI%*%t(Sig_L)+rep(1,No)%*%t(log(fs))
mu<-exp(xi_b)

## generate spike matrix N
N<-matrix(0,No,p);
for(i in 1:No){
  for(j in 1:p){
    N[i,j]<-rpois(n=1,lambda=mu[i,j])#+rpois(n=1,lambda=0.7);
  }
}

#### run Pf-PCA ########
#### recall the functions 
source("/Users/RONG/Documents/Neuro Projects/Neuronal Variability/Code/MainFunction.R")

##
post_results<-postPPs(N,k=length(s))
postlogmu<-post_results$postlogmu
##
FPCA_results<- FPCA(postlogmu,s,is_spline=T,period=180,nx=181,K=2)
meanmat<-FPCA_results$meanmat
pcmat<-FPCA_results$pcmat
scores<-FPCA_results$scores
varprop<-FPCA_results$varprop


##### save these results for S4: sharpening ####
N_S4_n50<-N; # the matrix of spike counts (50x9 matrix) 
alpha_S4_n50<-xi_bI[,1]; # true score (50-vector)
scores_S4_n50<-scores; # covered score (50-vector)
postlogmu_S4_n50<-postlogmu; #covered unobservable firing rates from spike counts in logspace (50x9 matrix)
meanmat_S4_n50<-meanmat; #covered mu_0 in logspace (181-vector)
pcmat_S4_n50<-pcmat; #covered phi(s) in logspace (181xK matrix)
varprop_S4_n50<-varprop; # the proportion of each components (K-vector)


###########################################
######## the format of saved data #########

#### combine all 4 cases and save them ####
postlogmu_4Cases<-rbind(postlogmu_S1_n50,postlogmu_S2_n50,postlogmu_S3_n50,postlogmu_S4_n50);
N_4Cases<-rbind(N_S1_n50,N_S2_n50,N_S3_n50,N_S4_n50);
alpha_4Cases<-rbind(alpha_S1_n50,alpha_S2_n50,alpha_S3_n50,alpha_S4_n50);
meanmat_4Cases<-rbind(meanmat_S1_n50,meanmat_S2_n50,meanmat_S3_n50,meanmat_S4_n50)
pcmat_4Cases<-rbind(pcmat_S1_n50,pcmat_S2_n50,pcmat_S3_n50,pcmat_S4_n50)
scores_4Cases<-rbind(scores_S1_n50,scores_S2_n50,scores_S3_n50,scores_S4_n50)
varprop_4Cases<-rbind(varprop_S1_n50,varprop_S2_n50,varprop_S3_n50,varprop_S4_n50,deparse.level = 0)

#### save the groundtruth of all 4 cases ####
sigma2_C<-c(c_S1,c_S2,c_S3,c_S4);
sigma2_180<-c(sigma_S1^2,sigma_S2^2,sigma_S3^2,sigma_S4^2);
phi_All<-rbind(phi1n,phi2n,phi3n,phi4n);
phi180_All<-rbind(phi180_1n,phi180_2n,phi180_3n,phi180_4n);



######################################################################
############ read results and plot them (Figure 2) ###################
No<-50;nx<-181
postlogmu_4Cases<-read.csv("/Users/rong/Documents/Neuro Projects/Neuronal Variability/Code/Simulation_Results/Simulation_postlogmu_4Cases_50x9.csv")
N_4Cases<-read.csv("/Users/RONG/Documents/Neuro Projects/Neuronal Variability/Code/Simulation_Results/SimulationData_N_4Cases_50x9.csv")
meanmat_4Cases<-read.csv("/Users/RONG/Documents/Neuro Projects/Neuronal Variability/Code/Simulation_Results/Simulation_meanmat_4Cases_50x9.csv")
pcmat_4Cases<-read.csv("/Users/RONG/Documents/Neuro Projects/Neuronal Variability/Code/Simulation_Results/SimulationData_pcmat_4Cases_50x9.csv")
alpha_4Cases<-read.csv("/Users/RONG/Documents/Neuro Projects/Neuronal Variability/Code/Simulation_Results/SimulationData_alpha_4Cases_50.csv")
scores_4Cases<-read.csv("/Users/RONG/Documents/Neuro Projects/Neuronal Variability/Code/Simulation_Results/SimulationData_scores_4Cases_50.csv")
varprop_4Cases<-read.csv("/Users/RONG/Documents/Neuro Projects/Neuronal Variability/Code/Simulation_Results/SimulationData_varprop_4Cases_50.csv")

sigma2_C<-read.csv("/Users/RONG/Documents/Neuro Projects/Neuronal Variability/Code/Simulation_Results/Simulation_Sigma2_S9.csv")
sigma2_180<-read.csv("/Users/RONG/Documents/Neuro Projects/Neuronal Variability/Code/Simulation_Results/Simulation_Sigma2_180.csv")
phi_All<-read.csv("/Users/RONG/Documents/Neuro Projects/Neuronal Variability/Code/Simulation_Results/Simulation_Phi_S9.csv")
phi180_All<-read.csv("/Users/RONG/Documents/Neuro Projects/Neuronal Variability/Code/Simulation_Results/Simulation_Phi_180.csv")
phi_All<-phi_All[,-1];sigma2_C<-sigma2_C[,-1];
phi180_All<-phi180_All[,-1];sigma2_180<-sigma2_180[,-1]


#### Bell-shape ####
sigma0<-20;
TuningCurve<-function(s,sigma=sigma0){
  fs0<-dnorm(s/sigma)
  c<-5/max(fs0);
  fs<-c*fs0+0.2;
  return(fs);
}
fs<-TuningCurve(s);
fs180<-TuningCurve(s180);
##############

s<-seq(-90,90,22.5);
s180<--90:90;
Fsize<-4; # four cases
ymax<-7.5; 
pdf(file="/Users/RONG/Documents/Neuro Projects/Neuronal Variability/Code/Simulation_Results/Simulation_fPCA_N50.pdf",width=4*1.7,height=Fsize*1.7);
par(mar=c(6,6,0.5,0.5),mfrow=c(Fsize,4),mex=.5,pty="s")
for(i in 1:Fsize){
  #true
  logfs<-log(fs180);phi1<-c(as.numeric(phi180_All[i,]));sigma<-sigma2_180[i]^0.5;
  if(i==Fsize){
    plot(s180,exp(logfs),xlab="stimulus",ylab="response (spikes/s)",type="l",lty="solid",cex.lab=1.4,pch=21,lwd=1.5,ylim=c(0,ymax),frame.plot=F,xaxt="n",yaxt="n");
    axis(side=1,at=seq(-90,90,45),labels=T,lwd=1.5,tcl=-0.2);
    axis(side=2,labels=T,lwd=1.5,tcl=-0.3);
  } else {
    plot(s180,exp(logfs),xlab="",ylab="response (spikes/s)",type="l",lty="solid",cex.lab=1.4,pch=21,lwd=1.5,ylim=c(0,ymax),frame.plot=F,xaxt="n",yaxt="n");
    axis(side=1,at=seq(-90,90,45),labels=F,lwd=1.5,tcl=-0.2);
    axis(side=2,labels=T,lwd=1.5,tcl=-0.3);
  }
  lines(s180,exp(logfs+sigma*phi1),type="l",lty="solid",col=2,lwd=1.5);
  lines(s180,exp(logfs-sigma*phi1),type="l",lty="solid",col=3,lwd=1.5);
  points(s,exp(log(fs)),pch=3)
  #fPCA
  f<-meanmat_4Cases[((i-1)*nx+1):(i*nx),2];
  phi<-pcmat_4Cases[((i-1)*nx+1):(i*nx),2];
  scores<-scores_4Cases[((i-1)*No+1):(i*No),2];
  prop<-varprop_4Cases[i,2]
  sigma<-(var(scores))^0.5
  sx<--90:90;
  if(i==Fsize){
    plot(sx,exp(f),xlab="stimulus",ylab="",type="l",lty="solid",cex.lab=1.4,pch=21,lwd=1.5,
         ylim=c(0,ymax),main=NULL,frame.plot=F,xaxt="n",yaxt="n");#paste("1st Component (",round(100*var[1]),"%)"));
    axis(side=1,at=seq(-90,90,45),labels=T,lwd=1.5,tcl=-0.2);
    axis(side=2,labels=F,lwd=1.5,tcl=-0.3);
  } else {
    plot(sx,exp(f),xlab="",ylab="",type="l",lty="solid",cex.lab=1.4,pch=21,lwd=1.5,
         ylim=c(0,ymax),main=NULL,frame.plot=F,xaxt="n",yaxt="n");#paste("1st Component (",round(100*var[1]),"%)"));
    axis(side=1,at=seq(-90,90,45),labels=F,lwd=1.5,tcl=-0.2)
    axis(side=2,labels=F,lwd=1.5,tcl=-0.3);
  }
  lines(sx,exp(f+sigma*phi),type="l",lty="solid",col=2,lwd=1.5);
  lines(sx,exp(f-sigma*phi),type="l",lty="solid",col=3,lwd=1.5);
  legend("topleft", legend=paste(round(100*prop,1),"%"),bty="n",cex=1.5);
  #PCA
  N_n<-as.matrix(N_4Cases[((i-1)*No+1):(i*No),-1])
  N.mean<-c(rep(1,No)%*%N_n)/No;
  svdN<-svd(var(N_n));
  svdN.u1<-svdN$u[,1];
  porp<-svd(var(N_n))$d
  svdvalue1<-round(100*porp[1]/sum(porp),1)
  if(i==Fsize){
    plot(s,N.mean,xlab="stimulus",ylab="",type="l",lty="solid",cex.lab=1.4,pch=21,lwd=1.5,
         ylim=c(0,ymax),main=NULL,frame.plot=F,xaxt="n",yaxt="n");
    axis(side=1,at=seq(-90,90,45),labels=T,lwd=1.5,tcl=-0.2);
    axis(side=2,labels=F,lwd=1.5,tcl=-0.3);
  } else {
    plot(s,N.mean,xlab="",ylab="",type="l",lty="solid",cex.lab=1.4,pch=21,lwd=1.5,
         ylim=c(0,ymax),main=NULL,frame.plot=F,xaxt="n",yaxt="n");
    axis(side=1,at=seq(-90,90,45),labels=F,lwd=1.5,tcl=-0.2);
    axis(side=2,labels=F,lwd=1.5,tcl=-0.3);
  }
  lines(s,N.mean+svdN.u1,type="l",lty="solid",col=3,lwd=1.5);
  lines(s,N.mean-svdN.u1,type="l",lty="solid",col=2,lwd=1.5);
  legend("topleft", legend=paste(svdvalue1,"%"),bty="n",cex=1.5);
  #Poisson-PCA
  postlogmu_n<-exp(as.matrix(postlogmu_4Cases[((i-1)*No+1):(i*No),-1]))
  N.mean<-c(rep(1,No)%*%postlogmu_n)/No;
  svdN.u1<-svd(var(postlogmu_n))$u[,1];
  porp<-svd(var(postlogmu_n))$d
  svdvalue1<-round(100*porp[1]/sum(porp),1)
  if(i==Fsize){
    plot(s,(N.mean),xlab="stimulus",ylab="",type="l",lty="solid",cex.lab=1.4,pch=21,lwd=1.5,
         ylim=c(0,ymax),main=NULL,frame.plot=F,xaxt="n",yaxt="n");
    axis(side=1,at=seq(-90,90,45),labels=T,lwd=1.5,tcl=-0.2);
    axis(side=2,labels=F,lwd=1.5,tcl=-0.3);
  } else {
    plot(s,(N.mean),xlab="",ylab="",type="l",lty="solid",cex.lab=1.4,pch=21,lwd=1.5,
         ylim=c(0,ymax),main=NULL,frame.plot=F,xaxt="n",yaxt="n");
    axis(side=1,at=seq(-90,90,45),labels=F,lwd=1.5,tcl=-0.2);
    axis(side=2,labels=F,lwd=1.5,tcl=-0.3);
  }
  lines(s,(N.mean+svdN.u1),type="l",lty="solid",col=3,lwd=1.5);
  lines(s,(N.mean-svdN.u1),type="l",lty="solid",col=2,lwd=1.5);
  legend("topleft", legend=paste(svdvalue1,"%"),bty="n",cex=1.5);
}
dev.off()


#####################################################
##################### plot scores ###################
pdf(file="/Users/RONG/Documents/Neuro Projects/Neuronal Variability/Code/Simulation_Results/Simulation_ScoreAlpha_N50.pdf",width=1*1.5,height=3*1.5);
par(mar=c(6,5,0.2,0.2)+.1,mfrow=c(3,1),mex=.4,pty="s")
##fPCA
No<-50;r2<-c();rho0<-c();
alpha_trueALL<-c()
score_predALL<-c()
for(i in 1:Fsize){
  score_pred<-scores_4Cases[((i-1)*No+1):(i*No),2]/(20)^0.5;
  alpha_true<-c(as.numeric(alpha_4Cases[i,-1]));
  score_predALL<-c(score_predALL,score_pred)
  alpha_trueALL<-c(alpha_trueALL,alpha_true);
  rho0[i]<-cor(score_pred,alpha_true)
  fit<-lm(score_pred~alpha_true);
  r2[i]<-round(summary(fit)$r.squared,2);
  b1<-summary(fit)$coefficient[2,1];
  if(b1<0){score_pred<--score_pred;
  fit<-lm(score_pred~alpha_true);
  r2[i]<-round(summary(fit)$r.squared,2);
  }
}
cor(score_predALL,alpha_trueALL);
ymin<-min(scores_4Cases[,2]/(20)^0.5);ymax<-max(scores_4Cases[,2]/(20)^0.5);
yscore<-c(as.matrix(alpha_4Cases[,-1]));xmin<-min(yscore);xmax<-max(yscore);mm<-max(-xmin,xmax)
score_pred<-scores_4Cases[1:No,2]/(20)^0.5;
alpha_true<-c(as.numeric(alpha_4Cases[1,-1]));
fit<-lm(score_pred~alpha_true);
r2[i]<-round(summary(fit)$r.squared,2);
b1<-summary(fit)$coefficient[2,1];
if(b1<0){score_pred<--score_pred;}
plot(alpha_true,score_pred,xlim=c(-mm,mm),ylim=c(-mm,mm),xlab="true score",ylab="predicted score",pch=21,bg = colors()[19],col="white",lwd = 0.5,frame.plot=F,xaxt="n",yaxt="n");
axis(side=1,labels=T,at=c(-2,0,2),tcl=-0.2)
axis(side=2,labels=T,at=c(-2,0,2),tcl=-0.2)
xx<--2+(0:10)/10*2*2;
lines(xx,xx,lty="dashed");
for(i in 2:Fsize){
  score_pred<-scores_4Cases[((i-1)*No+1):(i*No),2]/(20)^0.5;
  alpha_true<-c(as.numeric(alpha_4Cases[i,-1]));
  fit<-lm(score_pred~alpha_true);
  r2[i]<-round(summary(fit)$r.squared,2);
  b1<-summary(fit)$coefficient[2,1];
  if(b1<0){score_pred<--score_pred;}
  points(alpha_true,score_pred,pch=21,bg = colors()[18+i],lwd = 0.5,col="white");
}
legend("topleft", legend=c("M","A","T","W"),col=colors()[19:22], pch=19, cex=0.8);
## PCA
r2<-c();y1<-c();y2<-c();rho0<-c();
alpha_trueALL<-c()
score_predALL<-c()
for(i in 1:Fsize){
  N_n<-as.matrix(N_4Cases[((i-1)*No+1):(i*No),-1])
  prcN<-princomp(N_n,scores=T);
  score_pred<-prcN$scores[,1];
  alpha_true<-c(as.numeric(alpha_4Cases[i,-1]));
  score_predALL<-c(score_predALL,score_pred)
  alpha_trueALL<-c(alpha_trueALL,alpha_true);
  rho0[i]<-cor(score_pred,alpha_true)
  y1[i]<-min(score_pred);y2[i]<-max(score_pred);
  alpha_true<-c(as.numeric(alpha_4Cases[i,-1]));
  fit<-lm(score_pred~alpha_true);r2[i]<-round(summary(fit)$r.squared,2);
  b1<-summary(fit)$coefficient[2,1];
  if(b1<0){score_pred<--score_pred;
  fit<-lm(score_pred~alpha_true);r2[i]<-round(summary(fit)$r.squared,2);}
}
cor(score_predALL,alpha_trueALL);
N_n<-as.matrix(N_4Cases[1:No,-1])
prcN<-princomp(N_n,scores=T);
score_pred<-prcN$scores[,1]/sd(prcN$scores[,1]);
alpha_true<-c(as.numeric(alpha_4Cases[1,-1]));
fit<-lm(score_pred~alpha_true);
b1<-summary(fit)$coefficient[2,1];
if(b1<0){score_pred<--score_pred;}
ymin<-min(y1);ymax<-max(y2);
plot(alpha_true,score_pred,xlim=c(-mm,mm),ylim=c(-mm,mm),xlab="true score",ylab="predicted score",pch=21,lwd = 0.5,bg = colors()[19],col="white",frame.plot=F,xaxt="n",yaxt="n");
axis(side=1,labels=T,at=c(-2,0,2),tcl=-0.2)
axis(side=2,labels=T,at=c(-2,0,2),tcl=-0.2)
xx<--2+(0:10)/10*2*2;
lines(xx,xx,lty="dashed");
for(i in 2:Fsize){
  N_n<-as.matrix(N_4Cases[((i-1)*No+1):(i*No),-1])
  prcN<-princomp(N_n,scores=T);
  score_pred<-prcN$scores[,1]/sd(prcN$scores[,1]);
  alpha_true<-c(as.numeric(alpha_4Cases[i,-1]));
  fit<-lm(score_pred~alpha_true);
  b1<-summary(fit)$coefficient[2,1];
  if(b1<0){score_pred<--score_pred;}
  alpha_true<-c(as.numeric(alpha_4Cases[i,-1]));
  points(alpha_true,score_pred,bg = colors()[18+i],lwd = 0.5,col="white",pch=21);
}
## Poisson-PCA
r2<-c();y1<-c();y2<-c();rho0<-c();
alpha_trueALL<-c()
score_predALL<-c()
for(i in 1:Fsize){
  postlogmu_n<-as.matrix(postlogmu_4Cases[((i-1)*No+1):(i*No),-1])
  N_n<-postlogmu_n;
  prcN<-princomp(N_n,scores=T);
  score_pred<-prcN$scores[,1];
  alpha_true<-c(as.numeric(alpha_4Cases[i,-1]));
  score_predALL<-c(score_predALL,score_pred)
  alpha_trueALL<-c(alpha_trueALL,alpha_true);
  rho0[i]<-cor(score_pred,alpha_true)
  y1[i]<-min(score_pred);y2[i]<-max(score_pred);
  alpha_true<-c(as.numeric(alpha_4Cases[i,-1]));
  fit<-lm(score_pred~alpha_true);r2[i]<-round(summary(fit)$r.squared,2);
  b1<-summary(fit)$coefficient[2,1];
  if(b1<0){score_pred<--score_pred;
  fit<-lm(score_pred~alpha_true);r2[i]<-round(summary(fit)$r.squared,2);}
}
cor(score_predALL,alpha_trueALL);
postlogmu_n<-as.matrix(postlogmu_4Cases[1:No,-1])
N_n<-exp(postlogmu_n)
prcN<-princomp(N_n,scores=T);
score_pred<-prcN$scores[,1]/sd(prcN$scores[,1]);
alpha_true<-c(as.numeric(alpha_4Cases[1,-1]));
fit<-lm(score_pred~alpha_true);
b1<-summary(fit)$coefficient[2,1];
if(b1<0){score_pred<--score_pred;}
ymin<-min(y1);ymax<-max(y2);
plot(alpha_true,score_pred,xlim=c(-mm,mm),ylim=c(-mm,mm),xlab="true score",ylab="predicted score",pch=21,lwd = 0.5,bg = colors()[19],col="white",frame.plot=F,xaxt="n",yaxt="n");
axis(side=1,labels=T,at=c(-2,0,2),tcl=-0.2)
axis(side=2,labels=T,at=c(-2,0,2),tcl=-0.2)
xx<--2+(0:10)/10*2*2;
lines(xx,xx,lty="dashed");
for(i in 2:Fsize){
  postlogmu_n<-as.matrix(postlogmu_4Cases[((i-1)*No+1):(i*No),-1])
  N_n<-exp(postlogmu_n)
  prcN<-princomp(N_n,scores=T);
  score_pred<-prcN$scores[,1]/sd(prcN$scores[,1]);
  alpha_true<-c(as.numeric(alpha_4Cases[i,-1]));
  fit<-lm(score_pred~alpha_true);
  b1<-summary(fit)$coefficient[2,1];
  if(b1<0){score_pred<--score_pred;}
  alpha_true<-c(as.numeric(alpha_4Cases[i,-1]));
  points(alpha_true,score_pred,bg = colors()[18+i],lwd = 0.5,col="white",pch=21);
}
dev.off()


#########################################
### the setting of other three cases ####
### from 9 stimulus
#S1 Multplicative gain
phi1<-log(1.3*fs)-log(0.9*fs);
phi1n<-phi1/(sum(phi1^2))^0.5 
aa<-exp(phi1+log(fs))
c_S1<-sum((log(aa)-log(fs))^2)
### S2: additive gain
phi2<-log(0.4+fs)-log(fs);
phi2n<-phi2/(sum(phi2^2))^0.5
aa<-exp(phi2+log(fs))
c_S2<-sum((log(aa)-log(fs))^2)
### S3: tuning shift
phi3<-log(TuningCurve(s+5))-log(TuningCurve(s-5));
phi3n<-phi3/(sum(phi3^2))^0.5
aa<-exp(phi3+log(fs))
c_S3<-sum((log(aa)-log(fs))^2)


##### from 181 degrees from 0 to 180 ###
##### for plotting the groundtruth #####
s180<--90:90;fs180<-TuningCurve(s180)
### S1
phi180_1<-log(1.01*fs180)-log(0.99*fs180);
phi180_1n<-phi180_1/(sum(phi180_1^2))^0.5;
phi180_1M<-log(1.3*fs180)-log(0.9*fs180);
sigma_S1<-sum(phi180_1n*phi180_1M)
### S2
phi180_2<-log(0.2+fs180)-log(fs180-0.2);
phi180_2n<-phi180_2/(sum(phi180_2^2))^0.5;
phi180_2M<-log(fs180+0.4)-log(fs180);
sigma_S2<-sum(phi180_2n*phi180_2M)
### S3
phi180_3<-log(TuningCurve(s180+0.01))-log(TuningCurve(s180-0.01));
phi180_3n<-phi180_3/(sum(phi180_3^2))^0.5;
phi180_3M<-log(TuningCurve(s180+5))-log(TuningCurve(s180-5));
sigma_S3<-sum(phi180_3n*phi180_3M)





