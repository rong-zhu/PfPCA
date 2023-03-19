## Main functions of Pf-PCA
library(fda)
##############################################################
######################### main functions ######################
## the first step: obtain the posterior mean from EM algorithm 
postPP<-function(N,k=NULL){
  ## N: Spike count of an individual neuron at each trial and each stimulus, given some time windown. 
  ##    It's a matrix with trial x stimulus. 
  ## k: rank. It is the number of stimuli by default. 
  ##############
  ### iteration 
  epsilon<-1
  t<-1 ## the number of iteration
  SIG<-var(log(N+1)) # the initial covariance matrix Sigma 
  No<-length(N[,1]) # the number of trials 
  p<-length(N[1,]) # the dimension of stimulus
  B<-3000 # Monte Carlo
  fs_it<-c(log(1/No*rep(1,No)%*%N))
  if(is.null(k)){k<-p}
  #####
  while(epsilon>0.001){
    ## Monte Carlo
    postlogmu<-matrix(0,No,p)
    var.cond<-matrix(0,p,p)
    for(j in 1:No){
      xi_bI<-matrix(rnorm(B*k,0,1),B,k)
      a.eig <- eigen(SIG);
      if(k>1){Sig_L <- a.eig$vectors[,1:k] %*% diag(sqrt(a.eig$values[1:k]));
      } else {Sig_L <- a.eig$vectors[,1]*sqrt(a.eig$values[1]);}
      xi_b<-xi_bI%*%t(Sig_L)+rep(1,B)%*%t(fs_it)
      xi_b[xi_b==-Inf]<--10
      pjb<-c()
      for(b in 1:B){
        pjib<-c()
        for(i in 1:p){
          lambd<-exp(xi_b[b,i])
          pjib[i]<-exp(-lambd)*lambd^(min(N[j,i],200))/factorial(min(N[j,i],200))
        }
        pjib[is.na(pjib)]<-0
        pjib[pjib==Inf]<-0;
        pjb[b]<-exp(sum(log(pjib)))
      }
      sum_pjb<-Re(sum(pjb));
      if(sum_pjb<exp(-740)) {sum_pjb<-exp(-740);}
      postlogmu[j,]<-t(pjb)%*%xi_b/sum_pjb;
      var.cond<-var.cond+1/No*(t(xi_b)%*%diag(pjb)%*%xi_b/sum_pjb-postlogmu[j,]%*%t(postlogmu[j,]))
    }
    ###
    SIG.new<-t(postlogmu-rep(1,No)%*%t(rep(1,No))%*%postlogmu/No)%*%(postlogmu-rep(1,No)%*%t(rep(1,No))%*%postlogmu/No)/No+var.cond
    fs_it.new<-c(1/No*rep(1,No)%*%postlogmu)
    epsilon<-sum((SIG.new-SIG)^2)/sum(SIG^2);
    SIG<-SIG.new
    t<-t+1
    print(t);print(epsilon)
    # control the maximum iteration 10 for  the rare cases
    # in general just about 4 iterations are enough. For rare cases, the iterations are beyond 10.
    if(t>=10){epsilon<-0.00001} 
  }
  my_list<-list("postlogmu"=postlogmu,"Sigma"=SIG)
  return(my_list)
}
##########


###########################################################
## the second step: extract the Pf-PCA from the posteriors  
FPCA<- function(postlogmu,is_spline,s,period=180,nx=181,K=3){
  ## s:         stimulus 
  ## is_spline: nonparametric method, bspline (is_spline=T) or fourier (is_spline=F)
  ## period:    period for fourier method
  ## nx:        length of the components 
  ## K:         rank, i.e., the number of fPCs. 
  
  if(is_spline){
    ### using bspline basis
    norder=4;nbasis=length(s)+norder-2
    rangeval=c(min(s),max(s))
    heightbasis=create.bspline.basis(rangeval,nbasis=nbasis,norder=norder,s)
    heightfdPar=fdPar(heightbasis,2,0.01)
  } else {
    ### fourier 
    period <- period
    rangeval=c(min(s),max(s))
    heightbasis <- create.fourier.basis(rangeval, nbasis=9, period=period)
    harmaccelLfd <- vec2Lfd(c(0,(2*pi/diff(rangeval))^2,0), rangeval)
  }
  
  logfs_gamma<-postlogmu
  
  if(is_spline){loglam   = seq(-2,10,0.1);
  } else {loglam = seq(-2,10,0.1);}
  dfsave  <- rep(0,length(loglam))
  gcvsave <- rep(0,length(loglam))
  for (i in 1:length(loglam)){
    lambdai   = 10^loglam[i];
    if(is_spline){hgtfdPari = fdPar(heightbasis,2,lambdai)
    } else {hgtfdPari = fdPar(heightbasis,harmaccelLfd,lambdai) } # fourier
    smooth.fit =smooth.basis(s, y=t(logfs_gamma), hgtfdPari);
    gcvsave[i]=sum(smooth.fit[[3]])
    dfsave[i]=smooth.fit[[2]]
  }
  
  ##
  logmin<-loglam[which.min(gcvsave)]; # the chosen lambda by GCV in log space
  lambda<-10^logmin;
  if(is_spline){hgtfdPari = fdPar(heightbasis,2,lambda)
  } else {hgtfdPari = fdPar(heightbasis,harmaccelLfd,lambda)} 
  smooth.fit =smooth.basis(s, y=t(logfs_gamma), hgtfdPari);
  smooth.fd <- smooth.fit$fd
  
  #K=3
  smooth.pcalist = pca.fd(smooth.fd,K)
  
  ####### list the results of Functional PCA
  pcafd<-smooth.pcalist
  harmfd <- pcafd[[1]]
  var_prop<-pcafd$varprop # the variance proportions of fPCs
  scores<-pcafd$scores # the scores of fPCs
  basisfd <- harmfd$basis
  rangex <- basisfd$rangeval
  argvals <- seq(rangex[1], rangex[2], length = nx)
  meanmat <- eval.fd(argvals, pcafd$meanfd) # the mean in log space
  pcmat<-eval.fd(argvals, pcafd$harmonics) # the fPCs in log space
  
  pca_list<-list("meanmat"=meanmat,"pcmat"=pcmat,"varprop"=var_prop,"loglambd"=logmin,"scores"=scores)
  return(pca_list)
}

