library(nnet); # for MLR

#### model (relation between pc and mean) ####
#### regression parameter and width ##########
RegResults<-function(meanM,pc1M,ssize=8,onepeak=T){
  mean<-meanM;pc1<-pc1M;
  if(ssize==8){
    s<-round(22.5*(0:7))+1;
  } else {s<-round(30*(0:6))+1;
  }
  NeuronNo<-dim(mean)[2];
  beta_scaled<-matrix(0,NeuronNo,2);
  beta_unscaled<-beta_scaled;
  pvalue<-matrix(0,NeuronNo,2);
  R2<-c();R2b0<-c();R2Int<-c();
  fstat<-matrix(0,NeuronNo,3);
  LogmeanMax<-c();
  for(i in 1:NeuronNo){
    meand<-mean[,i]-max(mean[,i]); ## remove the maximum
    pc1d<-pc1[,i]; 
    pc1ds<-pc1d[s];meands<-meand[s]; 
    mean0<-mean[,i];
    fit_PC1toMeanUnScaled<-lm(pc1d[s]~mean0[s]);
    beta_unscaled[i,]<-summary(fit_PC1toMeanUnScaled)$coefficients[,1];
    ###
    LogmeanMax[i]<-max(mean[,i]);
    fit_PC1toMeanScaled<-lm(pc1d[s]~meand[s]);
    fit_PC1toMeanScaled_Summary<-summary(fit_PC1toMeanScaled)
    beta_scaled[i,]<-summary(fit_PC1toMeanScaled)$coefficients[,1];
    pvalue[i,]<-fit_PC1toMeanScaled_Summary$coefficients[,4];
    ###
    fit_PC1toMean_Summary<-summary(lm(pc1d~meand));
    R2[i]<-fit_PC1toMean_Summary$r.squared;
    fstat[i,]<-fit_PC1toMean_Summary$fstatistic;
    R2Int[i]<-1-sum((fit_PC1toMean_Summary$residuals)^2)/sum(pc1d^2)
  }
  ## list
  Reg<-list("fmax"=LogmeanMax,"beta_unscaled"=beta_unscaled,"beta_scaled"=beta_scaled,"pvalue"=pvalue,"R2"=R2,"R2Int"=R2Int,"fstat"=fstat);
  return(Reg)
}



###########################################
##### multinomial logistic regression  ####
MLR<-function(Counts_A,NumCV=5){
  Counts_NST0<-as.array(Counts_A,dim=c(Counts_A[1],Counts_A[2],Counts_A[3]))
  dims<-dim(Counts_NST0);
  X<-t(matrix(Counts_NST0,dims[1],prod(dims[2:3])));
  Y0<-c(1:dims[3])%x%rep(1,dims[2]);
  Y <- relevel(factor(Y0), ref = 1)
  Z = as.data.frame(cbind(X,Y));
  ndx<-sample(dim(Z)[1],dim(Z)[1]);
  prob<-c()
  for(cv in 1:NumCV){
    cvdx<-((dim(Z)[1]/NumCV)*(cv-1)+1):(cv*(dim(Z)[1]/NumCV));
    test_set<-ndx[cvdx];
    train_set<-ndx[-cvdx];
    Z_train<-Z[train_set,];
    Z_test<-Z[test_set,];
    X_test<-X[test_set,];
    Y_test<-Y[test_set];
    MLR_fit<-multinom(Y~.,data=Z, subset=train_set, maxit = 50000, trace=F)
    MLR_pred<-c(predict(MLR_fit,newdata=Z[test_set,]));
    prob[cv]<-sum(MLR_pred==Y_test)/length(Y_test);
  }
  return(mean(prob));
}


### function to divide high and low counts 
HighLow<-function(Counts_NST,bjs,orit){
  dims<-dim(Counts_NST);
  TrialNum<-c();
  for(t in 1:dims[2]){
    TrialNum[t]<-sum(Counts_NST[-bjs,t,orit])
  }
  Num_Neuron<-length(bjs)
  Counts_High<-array(0,c(Num_Neuron,dims[2]/2,dims[3]));
  Counts_Low<-array(0,c(Num_Neuron,dims[2]/2,dims[3]));
  Tra<-TrialNum;
  MFsort<-sort.int(Tra,index.return=T,decreasing=T);
  MFsort_idx<-MFsort$ix;
  for(i in 1:Num_Neuron){
    for(s in 1:dims[3]){
      iCounts<-Counts_NST[bjs[i],MFsort_idx,s]
      Counts_High[i,,s]<-iCounts[1:(dims[2]/2)];
      Counts_Low[i,,s]<-iCounts[(dims[2]/2+1):(dims[2])];
    }
  }
  ## list
  results<-list("High"=Counts_High,"Low"=Counts_Low);
  return(results)
}


############# flatness definition ##############
DELTA<-function(b,kappa){
  s<-(1:180-90)*pi/180;
  k<--log(0.2);
  mu<-exp(k*cos(s)-k)
  at<-0.2
  Y<-(mu^(1+kappa*at)*exp(b*at)-mu)-min(mu)*(exp(b*at)-1);
  ratio<-min(Y)/max(Y);
  return(ratio)
}


#### Fisher Information
FisherI<-function(meanmat_Neuron,pcmat_Neuron,scores_Neuron,a,b,model=TRUE,Support=FALSE){
  lengthTheta<-dim(meanmat_Neuron)[1]
  neuronNo<-dim(meanmat_Neuron)[2]
  trialNo<-dim(scores_Neuron)[1]
  K<-dim(pcmat_Neuron)[2]/dim(meanmat_Neuron)[2]
  #Ii_C<-array(0,dim=c(neuronNo,trialNo,lengthTheta-1))
  Ii_C<-matrix(0,neuronNo,trialNo)
  alpha1_matrix<-matrix(0,neuronNo,trialNo)
  for(i in 1:neuronNo){
    fi<-meanmat_Neuron[,i]
    pref_theta<-which.max(fi);
    phi<-pcmat_Neuron[,K*i-K+1];
    alpha<-scores_Neuron[,K*i-K+1];
    phi2<-pcmat_Neuron[,K*i-K+2];
    alpha2<-scores_Neuron[,K*i-K+2];
    phi3<-pcmat_Neuron[,K*i-K+3];
    alpha3<-scores_Neuron[,K*i-K+3];
    alpha1_matrix[i,]<-alpha
    for(k in 1:trialNo){
      if(model==FALSE){
        #a1<-diff(fi)+alpha[k]*diff(phi)+alpha2[k]*diff(phi2)+alpha3[k]*diff(phi3)
        #a2<-exp(fi[-1]+alpha[k]*phi[-1]+alpha2[k]*phi2[-1]+alpha3[k]*phi3[-1])
        a1_v0<-rep(0,lengthTheta-1);
        log_a2<-rep(0,lengthTheta-1);
        for(j in 1:3){
          a1_v0<-a1_v0+scores_Neuron[k,K*i-K+j]*diff(pcmat_Neuron[,K*i-K+j]);
          log_a2<-log_a2+scores_Neuron[k,K*i-K+j]*pcmat_Neuron[-1,K*i-K+j];
        }
        a1<-diff(fi)+a1_v0;
        a2<-exp(fi[-1]+log_a2);#exp(log_a2);#
      } else{
        a1<-diff(fi)*(1+alpha[k]*a[i])
        a2<-(exp(fi[-1]))^(alpha[k]*a[i]+1)*exp(alpha[k]*b[i])
      }
      aa<-a2*a1^2
      #Ii_C[i,k,]<-aa;
      if(Support==T){
        Ii_C[i,k]<-sum(aa[-c(1:22)]);#sum(aa[45:135])
      } else {Ii_C[i,k]<-sum(aa);}
    }
  }
  results<-list("Ii_C"=Ii_C,"alpha1_matrix"=alpha1_matrix)
  return(results)
}

##### shift the curve
ToShift <- function(x, lag) {
  n <- length(x)
  xnew <- rep(NA, n)
  if (lag < 0) {
    xnew[1:(n-abs(lag))] <- x[(abs(lag)+1):n]
    xnew[(n-abs(lag)+1):n] <- x[1:abs(lag)]
  } else if (lag > 0) {
    xnew[(lag+1):n] <- x[1:(n-lag)]
    xnew[1:abs(lag)] <- x[(n-abs(lag)+1):n]
  } else {
    xnew <- x
  }
  return(xnew)
}


############
############

