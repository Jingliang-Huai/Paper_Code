
rm(list=ls())
library("MTS")
library("MASS")
SummaryStatistics = function(data){
  moments = matrix(NA, ncol=ncol(data), nrow=14)
  colnames(moments)=colnames(data)
  rownames(moments)=c("Mean","Variance","Skewness","","Kurtosis","","JB","","ERS","","Q2(20)","","ARCH(20)","")
  for (i in 1:ncol(data)){
    moments[1,i] = mean(data[,i])
    moments[2,i] = var(data[,i])
    skew = moments::agostino.test(data[,i])
    moments[3,i] = skew$statistic[1]
    moments[4,i] = skew$p.value
    kurt = moments::anscombe.test(data[,i])
    moments[5,i] = kurt$statistic[1]-3
    moments[6,i] = kurt$p.value
    jb = moments::jarque.test(data[,i])
    moments[7,i] = jb$statistic
    moments[8,i] = jb$p.value
    ers = urca::ur.ers(data[,i],type="DF-GLS",model="constant")
    moments[9,i] = ers@teststat
    moments[10,i]= ers@testreg$coefficients[1,4]
    bt = WeightedPortTest::Weighted.Box.test(data[,i], type="Ljung-Box", lag=20)
    bt2 = WeightedPortTest::Weighted.Box.test(data[,i], type="Ljung-Box", lag=20, sqrd.res=T)
    moments[11,i] = bt2$statistic
    moments[12,i] = bt2$p.value
    bt3 = WeightedPortTest::Weighted.LM.test(data[,i], h.t=c(rep(var(data[,i]),nrow(data))), type="partial", lag=20)
    moments[13,i] = bt3$statistic
    moments[14,i] = bt3$p.value
  }
  
  cc=c(4,6,8,10,12,14)
  moments = round(moments,3)
  moments1 = moments
  for (j in 1:k){
    for (i in 1:length(cc)){
      i = cc[i]
      if (moments[i,j]<=0.01) {
        moments1[(i-1),j] = paste(format(round(moments[(i-1),j],3),nsmall=3),"***",sep="")
        moments1[i,j] = paste("(",format(round(moments[i,j],3),nsmall=3),")",sep="")
      } else if (moments[i,j]<=0.05) {
        moments1[(i-1),j] = paste(format(round(moments[(i-1),j],3),nsmall=3),"**",sep="")
        moments1[i,j] = paste("(",format(round(moments[i,j],3),nsmall=3),")",sep="")
      } else if (moments[i,j]<=0.10) {
        moments1[(i-1),j] = paste(format(round(moments[(i-1),j],3),nsmall=3),"*",sep="")
        moments1[i,j] = paste("(",format(round(moments[i,j],3),nsmall=3),")",sep="")
      } else {
        moments1[(i-1),j] = format(round(moments[(i-1),j],3),nsmall=3)
        moments1[i,j] = paste("(",format(round(moments[i,j],3),nsmall=3),")",sep="")
      }
    }
  }
  
  for (j in 1:k){
    i = 9
    if (moments[i,j]<=-2.57) {
      moments1[i,j] = paste(format(round(moments[i,j],3),nsmall=3),"***",sep="")
    } else if (moments[i,j]<=-1.96) {
      moments1[i,j] = paste(format(round(moments[i,j],3),nsmall=3),"**",sep="")
    } else if (moments[i,j]<=-1.62) {
      moments1[i,j] = paste(format(round(moments[i,j],3),nsmall=3),"*",sep="")
    } else {
      moments1[i,j] = format(round(moments[i,j],3),nsmall=3)
    }
  }
  moments1
}
Table = function(p){
  para = as.matrix(p)
  if (sum(is.nan(p)==TRUE)==0) {
    k = ncol(para)
    kk = nrow(para)/2
    for (i in 1:kk){
      for (j in 1:k){
        if (para[i,j]==0){
          next
        } else if (abs((para[2*i,j]))<=0.01) {
          para[i,j] = paste(format(round(as.numeric(para[i,j]),3),nsmall=3),"***",sep="")
          para[(kk+i),j] = paste("(",format(round(as.numeric(para[(kk+i),j]),3),nsmall=3),")",sep="")
        } else if (abs((para[2*i,j]))<=0.05) {
          para[i,j] = paste(format(round(as.numeric(para[i,j]),3),nsmall=3),"**",sep="")
          para[(kk+i),j] = paste("(",format(round(as.numeric(para[(kk+i),j]),3),nsmall=3),")",sep="")
        } else if (abs((para[2*i,j]))<=0.10) {
          para[i,j] = paste(format(round(as.numeric(para[i,j]),3),nsmall=3),"*",sep="")
          para[(kk+i),j] = paste("(",format(round(as.numeric(para[(kk+i),j]),3),nsmall=3),")",sep="")
        } else if (abs((para[2*i,j]))>0.1) {
          para[i,j] = paste(format(round(as.numeric(para[i,j]),3),nsmall=3),sep="")
          para[(kk+i),j] = paste("(",format(round(as.numeric(para[(kk+i),j]),3),nsmall=3),")",sep="")
        }
      }
    }
  }
  para
}
SignBias_WARCH = function(ugarch.fit,lag=20) {
  sign.bias = signbias(ugarch.fit)[1,][1:2]
  warch = Weighted.LM.test(residuals(ugarch.fit), sigma(ugarch.fit)^2, lag=lag, type=c("correlation"), fitdf=2, weighted=TRUE)
  pval = c(sign.bias[2]$prob,warch$p.value)
  x = Table(matrix(as.matrix(sign.bias),ncol=1))
  y = Table(matrix(c(warch$statistic,warch$p.value),ncol=1))
  statistic = rbind(x,y)
  rownames(statistic) = c("Sign Bias","",paste0("WARCH(",lag,")"),"")
  return = list(statistic=statistic, pval=pval)
}
ValueAtRisk = function(ugarch.fit, ugarch.spec, prob=0.05, conf.level=0.90){
  setfixed(ugarch.spec) = as.list(coef(ugarch.fit))
  data = ugarch.fit@model$modeldata$data
  ugarch.filter = ugarchfilter(ugarch.spec, data, n.old=length(data))
  VaR = fitted(ugarch.filter) + sigma(ugarch.filter)*qdist(ugarch.fit@model$modeldesc$distribution, p=prob, mu=ugarch.fit@fit$matcoef[1,1], sigma=1,
                                         skew=ifelse(is.na(coef(ugarch.fit)["skew"]),0,coef(ugarch.fit)["skew"]), shape=coef(ugarch.fit)["shape"])
  var_test = VaRTest(prob, as.numeric(data), as.numeric(VaR), conf.level=conf.level)
  vardur_test = VaRDurTest(p, as.numeric(data), as.numeric(VaR),conf.level=conf.level)
  f = function(x) {
    qdist(ugarch.fit@model$modeldesc$distribution, p=x, mu=0,sigma=1,skew=ifelse(is.na(coef(ugarch.fit)["skew"]),0,coef(ugarch.fit)["skew"]), shape=coef(ugarch.fit)["shape"])
  }
  ES = fitted(ugarch.filter) + sigma(ugarch.filter)*integrate(f,0,prob)$value/prob
  ES = ESTest(prob, as.numeric(data), as.numeric(ES), VaR, boot=TRUE, n.boot=1000, conf.level=conf.level)
  decision = ifelse(ES$boot.p.value>0.10,"H0","H1")
  x = Table(c(var_test$uc.LRstat, var_test$uc.LRp))
  y = Table(c(vardur_test$rLL,vardur_test$LRp))
  z = Table(c(1, ES$boot.p.value))
  z[1,1] = decision
  statistic = rbind(x,y,z)
  pval = c(var_test$uc.LRp,round(ES$boot.p.value,3),vardur_test$LRp)
  rownames(statistic) = c("VaR","","CVaR","","VaR Dur.","")
  return = list(statistic=statistic,pval=pval)
}
InformationCriterion = function (ugarch.fit, ugarch.spec, prob=0.05, conf.level=0.90, lag=20) {
  qprob = qnorm(1-prob)
  loss = sum(abs(ugarch.fit@fit$robust.tval[-c(1:2)])<=qprob)
  if (is.na(ugarch.fit@fit$robust.matcoef[1,2])==FALSE) {
    if ("skew" %in% rownames(ugarch.fit@fit$robust.matcoef)) {
      upper = ugarch.fit@fit$robust.matcoef["skew",1] + qprob*ugarch.fit@fit$robust.matcoef["skew",2]
      lower = ugarch.fit@fit$robust.matcoef["skew",1] - qprob*ugarch.fit@fit$robust.matcoef["skew",2]
      if (upper>1 && lower<1) {
        loss = loss + 100
      }
    }
  }
  var = ValueAtRisk(ugarch.fit, ugarch.spec, prob, conf.level)$pval
  sbwarch = SignBias_WARCH(ugarch.fit, lag=lag)$pval
  
  t = length(ugarch.fit@fit$z)
  IC = -2*likelihood(ugarch.fit) + loss*log(t)
  IC = IC + sum(c(var,sbwarch)<0.10)*10^5
  IC = ifelse(is.na(IC), 10^8, IC)
  IC
}
BestGARCH = function(distr=c("norm","snorm","std","sstd","ged","sged"), models=c("sGARCH","eGARCH","gjrGARCH","iGARCH","TGARCH","AVGARCH","NGARCH","NAGARCH","APARCH","ALLGARCH"), data, ar=0, ma=0, prob=0.05, conf.level=0.90, lag=20){
  data = matrix(data, ncol=1)
  GARCH_IC = matrix(10^7, nrow=length(distr), ncol=length(models))
  colnames(GARCH_IC) = models
  rownames(GARCH_IC) = distr
  spec_list = list()
  for (i in 1:length(models)) {
    spec_list[[i]] = list()
  }
  for (j in 1:length(models)) {
    print(paste0("-",models[j]))
    for (i in 1:length(distr)) {
      print(paste0("--",distr[i]))
      if (models[j] %in% c("AVGARCH","TGARCH","APARCH","NAGARCH","NGARCH","ALLGARCH")) {
        #if (models[j]=="AVGARCH") {
        #  fixed.pars=list(gamma1=0,gamma2=0,delta=1)
        #} else if (models[j]=="TGARCH") {
        #  fixed.pars=list(delta=1)
        #} else {
        #  fixed.pars=list()
        #}
        ugarch.spec = ugarchspec(mean.model=list(armaOrder=c(ar,ma)),
                                 variance.model=list(model="fGARCH", submodel=models[j], garchOrder=c(1,1), variance.targeting=FALSE), 
                                 distribution.model=distr[i])
                                 #fixed.pars=fixed.pars)
      } else {
        ugarch.spec = ugarchspec(mean.model=list(armaOrder=c(ar,ma)),
                                 variance.model=list(model=models[j], garchOrder=c(1,1), variance.targeting=FALSE),
                                 distribution.model=distr[i])
      }
      ugarch.fit = ugarchfit(ugarch.spec, data, solver="hybrid", solver.list=list(outer.iter=10, inner.iter=1000))
      if (ugarch.fit@fit$convergence==0) {
        GARCH_IC[i,j] = InformationCriterion(ugarch.fit, ugarch.spec, prob=prob, conf.level=conf.level, lag=lag)
        spec_list[[i]][[j]] = ugarch.spec
      }
    }
  }
  return=list(GARCH_IC=GARCH_IC, spec_list=spec_list)
}
BayesPrior = function(Y, nlag){
  k = ncol(Y)
  vars = MTS::VAR(Y, p=nlag, include.mean=TRUE, output=FALSE)
  varcoef = t(vars$Phi)
  SIGMA_OLS = vars$secoef
  Q_0 = vars$Sigma
  b_prior = varcoef
  beta_0.var = diag(c(vars$secoef[-(k+1),]))^2
  return=list(aprior=b_prior,Vprior=beta_0.var,Q_0=Q_0)
}

MinnesotaPrior = function(gamma, r, nlag){
  m = nlag*(r^2)
  A_prior = cbind(0*diag(r), matrix(0, ncol=(nlag-1)*r, nrow=r))
  aprior = c(A_prior)
  V_i = matrix(0, nrow=(m/r), ncol=r)
  for (i in 1:r){
    for (j in 1:(m/r)) {
      V_i[j,i] = gamma/(ceiling(j/r)^2)
    }
  }
  # Now V (MINNESOTA VARIANCE) is a diagonal matrix with diagonal elements the V_i'  
  V_i_T = t(V_i)
  Vprior = diag(c(V_i_T))
  diag(Vprior)
  return = list(aprior=aprior, Vprior=Vprior)
}
TVPVAR = function(Y, l, nlag, prior){
  beta_0.mean = prior$aprior
  beta_0.var = prior$Vprior
  Q_0 = prior$Q_0
  if (is.null(Q_0)) {
    Q_0 = cov(Y)
  }
  
  create_RHS_NI = function(templag, r, nlag, t){
    K = nlag*(r^2)
    x_t = matrix(0, (t-nlag)*r, K)
    for (i in 1:(t-nlag)){
      ztemp=NULL
      for (j in 1:nlag){
        xtemp = templag[i,((j-1)*r+1):(j*r)]
        xtemp = t(kronecker(diag(r),xtemp))
        ztemp = cbind(ztemp, xtemp)
      }
      x_t[((i-1)*r+1):(i*r),] = ztemp
    }
    return=list(x_t=x_t, K=K)
  }
  Y = scale(Y,TRUE,FALSE)
  y_true = 0
  FPC = Y
  YX = cbind(Y,Y)
  nfac = 0
  p = n = ncol(Y)
  r = nfac + p
  m = nlag*(r^2)
  k = nlag*r
  t = nrow(FPC)
  q = n + p
  Q_0 = Q_0
  
  # Initialize matrices
  beta_0_prmean = beta_0.mean
  beta_0_prvar = beta_0.var
  
  beta_pred = matrix(0,m,t)
  beta_update = matrix(0,m,t)
  
  Rb_t = array(0,c(m,m,t))
  Sb_t = array(0,c(m,m,t))
  
  beta_t = array(0, c(k,k,t))
  Q_t = array(0, c(r,r,t))
  
  # Decay and forgetting factors
  l_2 = l[2]
  l_4 = l[1]
  
  # Define lags of the factors to be used in the state (VAR) equation         
  yy = FPC[(nlag+1):t,]      
  xx = embed(FPC,nlag+1)[,-c(1:ncol(FPC))]
  templag = embed(FPC,nlag+1)[,-c(1:ncol(FPC))]
  RHS1 = create_RHS_NI(templag,r,nlag,t);  
  Flagtemp = RHS1$x_t
  m = RHS1$K
  Flag = rbind(matrix(0, k,m), Flagtemp)
  
  ###-----| 1. KALMAN FILTER
  for (irep in 1:t){
    #-----| Update the state covariances
    # 1. Get the variance of the factor
    
    # Update Q[t]
    if (irep==1){
      Q_t[,,irep] = Q_0
    } else if (irep > 1) {
      if (irep <= (nlag+1)) { 
        Gf_t = 0.1*(t(matrix(FPC[irep,],nrow=1))%*%(FPC[irep,]))
      } else {
        Gf_t = t(yy[(irep-nlag),]-xx[(irep-nlag),]%*%t(B[1:r,1:k])) %*% (yy[(irep-nlag),]-xx[(irep-nlag),]%*%t(B[1:r,1:k]))
      }
      Q_t[,,irep] = l_2*Q_t[,,(irep-1)] + (1-l_2)*Gf_t[1:r,1:r]
    }
    # -for beta
    if (irep <= (nlag+1)) {
      beta_pred[,irep] = beta_0_prmean
      beta_update[,irep] = beta_pred[,irep]
      Rb_t[,,irep] = beta_0_prvar
    } else if (irep > (nlag+1)) {
      beta_pred[,irep] = beta_update[,(irep-1)]
      Rb_t[,,irep] = (1/l_4)*Sb_t[,,(irep-1)]
    }
    
    # -for beta
    if (irep >= (nlag+1)) {
      # 2/ Update VAR coefficients conditional on Principal Componets estimates
      Rx = Rb_t[,,irep]%*%t(Flag[((irep-1)*r+1):(irep*r),])
      KV_b = Q_t[,,irep] + Flag[((irep-1)*r+1):(irep*r),]%*%Rx
      KG = Rx%*%MASS::ginv(KV_b)
      beta_update[,irep] = matrix(beta_pred[,irep], ncol=1) + (KG%*%(t(matrix(FPC[irep,], nrow=1))-Flag[((irep-1)*r+1):(irep*r),]%*%matrix(beta_pred[,irep], ncol=1)) )
      Sb_t[,,irep] = Rb_t[,,irep] - KG%*%(Flag[((irep-1)*r+1):(irep*r),]%*%Rb_t[,,irep])
    }
    
    # Assign coefficients
    bb = matrix(beta_update[,irep], ncol=1)
    splace = 0
    biga = matrix(0, r,r*nlag)
    for (ii in 1:nlag) {                                          
      for (iii in 1:r) {           
        biga[iii,((ii-1)*r+1):(ii*r)] = t(bb[(splace+1):((splace+r)),1])
        splace = splace + r
      }
    }
    
    B = rbind(biga, cbind(diag(r*(nlag-1)), matrix(0, nrow=r*(nlag-1), ncol=r)))
    
    if ((max(abs(eigen(B)$values))<=1)||(irep==1)){
      beta_t[,,irep] = B
    } else {
      beta_t[,,irep] = beta_t[,,(irep-1)]
      beta_update[,irep] = 0.99*beta_update[,(irep-1)]
    }
  }
  
  return = list(beta_t=beta_t[1:ncol(Y),,], Q_t=Q_t)
}
GFEVD = function(Phi, Sigma, n.ahead=10,normalize=TRUE,standardize=TRUE) {
  tvp.Phi = function (x, nstep = 10, ...) {
    nstep = abs(as.integer(nstep))
    K=nrow(x)
    p=floor(ncol(x)/K)
    A = array(0, c(K,K,nstep))
    for (i in 1:p){
      A[,,i]=x[,((i-1)*K+1):(i*K)]
    }
    
    Phi = array(0, dim = c(K, K, nstep + 1))
    Phi[, , 1] = diag(K)
    Phi[, , 2] = Phi[, , 1] %*% A[, , 1]
    if (nstep > 1) {
      for (i in 3:(nstep + 1)) {
        tmp1 = Phi[, , 1] %*% A[, , i - 1]
        tmp2 = matrix(0, nrow = K, ncol = K)
        idx = (i - 2):1
        for (j in 1:(i - 2)) {
          tmp2 = tmp2 + Phi[, , j + 1] %*% A[, , idx[j]]
        }
        Phi[, , i] = tmp1 + tmp2
      }
    }
    return(Phi)
  }
  A = tvp.Phi(Phi, (n.ahead-1))
  Sigma = Sigma
  gi = array(0, dim(A))
  sigmas = sqrt(diag(Sigma))
  for (j in 1:dim(A)[3]) {
    gi[,,j] = t(A[,,j]%*%Sigma%*%MASS::ginv(diag(sqrt(diag(Sigma)))))
  }
  if (standardize==TRUE){
    girf=array(NA, c(dim(gi)[1],dim(gi)[2], (dim(gi)[3])))
    for (i in 1:dim(gi)[3]){
      girf[,,i]=((gi[,,i])%*%MASS::ginv(diag(diag(gi[,,1]))))
    }
    gi=girf
  }
  
  num = apply(gi^2,1:2,sum)
  den = c(apply(num,1,sum))
  fevd = t(num)/den
  nfevd = fevd
  if (normalize==TRUE) {
    fevd=(fevd/apply(fevd, 1, sum))
  } else {
    fevd=(fevd)
  }
  return = list(GFEVD=fevd, GIRF=gi)
}


#ConnectednessTable= function (FEVD, digit = 2) 
#{
#  if (length(dim(FEVD)) <= 1) {
#    stop("FEVD needs to be at least a 2-dimensional matrix")
#  }
#  NAMES = colnames(FEVD)
#  k = dim(FEVD)[1]
#  if (is.null(NAMES)) {
#    NAMES = 1:k
#  }
#  CT = apply(FEVD, 1:2, mean) * 100
#  OWN = diag(diag(CT))
#  TO = colSums(CT - OWN)
#  FROM = rowSums(CT - OWN)
#  NET = TO - FROM
#  TCI = mean(TO)
#  cTCI = TCI * k/(k - 1)
#  NPDC = CT - t(CT)
#  NPT = rowSums(NPDC < 0)
#  INFLUENCE = 100 * abs(NPDC/t(t(CT) + CT))
#  table = format(round(cbind(CT, FROM), digit), nsmall = digit)
#  to = c(format(round(c(TO, sum(TO)), digit), nsmall = digit))
#  inc = c(format(round(colSums(CT), digit), nsmall = digit), 
#          "cTCI/TCI")
#  tci = paste0(format(round(cTCI, digit), nsmall = digit), 
#               "/", format(round(TCI, digit), nsmall = digit))
#  net = c(format(round(NET, digit), nsmall = digit))
#  net = c(net, tci)
#  npt = c(format(round(NPT, digit), nsmall = digit), "")
#  TABLE = rbind(table, to, inc, net, npt)
#  colnames(TABLE) = c(NAMES, "FROM")
#  rownames(TABLE) = c(NAMES, "TO", "Inc.Own", "NET", "NPT")
#  PCI = matrix(NA, k, k)
#  for (i in 1:k) {
#    for (j in 1:k) {
#      PCI[i, j] = 200 * (CT[i, j] + CT[j, i])/(CT[i, i] + 
#                                                 CT[i, j] + CT[j, i] + CT[j, j])
#    }
#  }
#  return = list(FEVD = CT, TCI = TCI, cTCI = cTCI, PCI = PCI, 
#                TO = TO, FROM = FROM, NET = NET, NPDC = NPDC, TABLE = TABLE, 
#                NPT = NPT, INFLUENCE = INFLUENCE)
#}


DCA = function(CV, digit=2){
  k = dim(CV)[1]
  CT = apply(CV,1:2,mean)*100 # spillover from others to one specific
  OWN = diag(diag(CT))
  TO = colSums(CT-OWN)
  FROM = rowSums(CT-OWN)
  NET = TO-FROM
  TCI = mean(TO)
  NPSO = CT-t(CT)
  NPDC = rowSums(NPSO>0)
  INFLUENCE = abs(NPSO/t(t(CT)+CT))
  table = format(round(cbind(CT,FROM),digit),nsmall=digit)
  to = c(format(round(c(TO,sum(TO)),digit),nsmall=digit))
  net = c(format(round(c(NET, TCI),digit),nsmall=digit))
  npdc = c(format(round(NPDC,digit),nsmall=digit), "")
  inc = c(format(round(colSums(CT), digit),nsmall=digit), "TCI")
  TABLE = rbind(table,to,inc,net,npdc)
  colnames(TABLE) = c(rownames(CV),"FROM others")
  rownames(TABLE) = c(rownames(CV),"TO others","Inc. own","NET","NPDC")
  PCI = matrix(NA, k, k)
  for (i in 1:k) {
    for (j in 1:k) {
      PCI[i,j] = 2*(CT[i,j]+CT[j,i])/(CT[i,i]+CT[i,j]+CT[j,i]+CT[j,j])
    }
  }
  return = list(CT=CT,TCI=TCI,TCI_corrected=TCI*k/(k-1),PCI=PCI,TO=TO,FROM=FROM,NET=NET,NPSO=NPSO,NPDC=NPDC,INFLUENCE=INFLUENCE,TABLE=TABLE)
}


##################################################################################################
##################################################################################################
##################################################################################################

library("abind")
library("rmgarch")
library("openxlsx")
library("parallel")
library("RColorBrewer")
library("WeightedPortTest")
library("psych")

#source("functions.R")

options(warn=-1)
options("mc.cores"=detectCores())
colors = c("black","tan4","springgreen4","springgreen2","steelblue4","steelblue1","maroon4","maroon1","orangered4","orangered","tan1","yellow3")
palette(colors)

DATA = read.xlsx("./green_and_non.xlsx", "Sheet1", detectDates=TRUE)
DATE = DATA[,1]
RAW = DATA[,-1]  ## dropping first column  
k = ncol(RAW)

### BASE ALL ON OTHER CURRENCY
#colnames(RAW)[1] = c("EUR(DM)")
NAMES=colnames(RAW)

# CALCULATE FIRST DIFFERENCED 1Y IRS SERIES
data = RAW[-1,]  ## dropping first observation 
colnames(data)=NAMES
for (i in 1:k) {
  data[,i] = diff(log(RAW[,i]))*100
}



### PRINCIPAL COMPONENT ANALYSIS
fa = principal(data, nfactor=1, method="tenBerge")
Y = cbind(data,(fa$scores))
NAMES = c(NAMES,"ME")
colnames(Y) = NAMES
k = ncol(Y)
#t = length(DATE[-1])
date=DATE[-1]


split = 3
par(mfcol=c(ceiling((k-1)/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:(k-1)){
   plot(DATE[-1], Y[,k], type="l", main=NAMES[i], xaxs="i", las=1, col="orange",ylim=c(-5,5),xlab="",ylab="")
   grid(NA,NULL)
   lines(DATE[-1], Y[,i], col="steelblue4")
   abline(h=0)
}

rownames(Y) = date
Y = as.data.frame(ts(na.omit(Y)))

digit = 4
SummaryStatistics_list = rbind(SummaryStatistics(na.omit(Y)), format(round(cor(na.omit(Y)),digit),nsmall=digit))
#options(max.print = 10000) 
#increase print limit to max allowed by your machine
options(max.print = .Machine$integer.max)
SummaryStatistics_list


##### DCC-COPULA ANALYSIS
lag = 20  #20 original --Lag must exceed fitted degrees of freedom
#lag = nlag
prob = 0.05
conf.level = 0.90
R = H = NULL
copula_fit = GARCH_spec = evaluation_list = list()

#ugarch.fit.results = list()
   spec = c()
   evaluation_matrix = matrix(NA, ncol=k, nrow=10)
   colnames(evaluation_matrix) = NAMES
   rownames(evaluation_matrix) = c("VaR","","CVaR","","VaR Dur.","","Sign Bias","",paste0("WARCH(",lag,")"),"")
   data = Y
   #data=data_list[[3]][-row(data_list[[3]])[data_list[[3]] == 0],]  ## dropping Zeros from the data
   data=data[-row(data)[data == 0],]  ## dropping Zeros from the data
   data = na.omit(data)
   for (j in 1:k) {
      print(NAMES[j])
      bestgarch = BestGARCH(data=data[,j],prob=prob,conf.level=conf.level,lag=lag,
                            distr=c("norm","snorm","std","sstd","ged","sged"), 
                            models=c("sGARCH","iGARCH","eGARCH","gjrGARCH","AVGARCH","TGARCH"))
      GARCH_selection = which(bestgarch$GARCH_IC==min(bestgarch$GARCH_IC),arr.ind=TRUE)
      print(paste(colnames(bestgarch$GARCH_IC)[GARCH_selection[2]], rownames(bestgarch$GARCH_IC)[GARCH_selection[1]]))
      print(bestgarch$GARCH_IC)
      ugarch.spec = bestgarch[[2]][[GARCH_selection[1]]][[GARCH_selection[2]]]
      ugarch.fit = ugarchfit(ugarch.spec,data=data[,j])
      print(ugarch.fit)
      evaluation_matrix[,j] = rbind(ValueAtRisk(ugarch.fit,ugarch.spec,prob=prob, conf.level=conf.level)$statistic,
                                    SignBias_WARCH(ugarch.fit,lag=lag)$statistic)
      spec = c(spec,ugarch.spec)
   }
   print(evaluation_matrix)
   evaluation_list = evaluation_matrix
   GARCH_spec = spec
   
   
   data = as.matrix(Y)
   #data=data_list[[3]][-row(data_list[[3]])[data_list[[3]] == 0],]  ## dropping Zeros from the data
   #data=data[-row(data)[data == 0],]  ## dropping Zeros from the data
   data = na.omit(data)
   mgarch.spec = cgarchspec(uspec=multispec(spec), dccOrder=c(1,1), asymmetric=FALSE,  
                            distribution.model=list(copula="mvt", method="Kendall", time.varying=TRUE, transformation="parametric"))
   copula_fit = cgarchfit(mgarch.spec, data=data, solver=c("hybrid","solnp"), fit.control=list(eval.se=TRUE))
   copula_fit
   H = abind(H,rcov(copula_fit))
   R = abind(R,rcor(copula_fit))
   evaluation_list
   summary(R)
   

   dcc.garch.spec = dccspec(uspec=multispec(spec), dccOrder=c(1,1), distribution = "mvt", model="DCC")
  
       # Model 1: Symmetric GARCH-DCC:
   dcc.garch.fit = dccfit(dcc.garch.spec, data = data)
   resid.dcc=dcc.garch.fit@mfit$stdresid
   residns.dcc=dcc.garch.fit@model$residuals
   rcov.dcc=rcov(dcc.garch.fit)
   

   #Diagnostic tests:
   colnames(resid.dcc)=colnames(data)
   
   diagnostic = function(x,residns,rcov) {
     sss = NULL
     for (i in 1:ncol(x)) {
       box=Box.test(x[,i], type="Ljung-Box", lag = 30)
       boxs =ifelse(box$p.value<0.01,paste(format(round(box$statistic,1),nsmall=1),"*"),format(box$statistic,nsmall=1))
       boxp=round(box$p.value,2)
       lm=Weighted.LM.test(residns[,i], rcov[i,i,], lag = 30, type = c("partial"), fitdf = 1, weighted = F)
       lms = ifelse(lm$p.value<0.01,paste(format(round(lm$statistic,1),nsmall=1),"*"),format(lm$statistic,nsmall=1))
       lmp=round(lm$p.value,2)
       ss = c(boxs, boxp, lms, lmp)
       sss = rbind(sss,ss)
     }
     rownames(sss) = colnames(x)
     colnames(sss) = c("Q(30)","p-value","Li-Mak(30)","p-value")
     sss
   }
   
   diag.dcc=diagnostic(resid.dcc, residns.dcc, rcov.dcc)
   diag.dcc
   
   
### DYNAMIC CONDITIONAL CORRELATIONS
t = length(DATE[-1])
date=DATE[-c(1:(length(Y[,1])-length(R[2,1,])+1))]

par(mfrow = c(ceiling(ncol(R)/4),6), oma = c(0,1,0,0) + 0.05, mar = c(1,1,1,1) + .05, mgp = c(0, 0.1, 0))
kk = k*(k-1)/2
for (i in 1:k) {
   for (j in 1:k) {
      if (i<j) {
         plot(date,R[j,i,], xlab="",ylab="",type="l",xaxs="i",col="steelblue4", las=1, main=paste0(NAMES[i],"-",NAMES[j]),tck=-0.02,yaxs="i",ylim=c(min(R),1))
         grid(NA,NULL)
         lines(date,R[j,i,],col="steelblue4")
         #print(colnames(mean(R[j,i,])) = paste0(NAMES[i],"-",NAMES[j]))
         t.test(R[j,i,])
         abline(h=0,lty=2)
         box()
      }
   }
}


dcc_corr=list()
jk = 1
for (i in 1:k) {
  for (j in 1:k) {
    temp=R[j,i,]
    dcc_corr=cbind(dcc_corr,temp)
    colnames(dcc_corr)[jk] = paste0(colnames(Y)[j],"-",colnames(Y)[i])
    jk = jk + 1
  }
}
head(dcc_corr)




Y = matrix(NA,ncol=k, nrow=dim(H)[3])
for (i in 1:k) {
   Y[,i] = log(H[i,i,])
}
colnames(Y) = NAMES

#####################################################
# VAR with liberary vars 
#####################################################
#VAR(y, p = 1, type = c("const", "trend", "both", "none"),
#season = NULL, exogen = NULL, lag.max = NULL,
#ic = c("AIC", "HQ", "SC", "FPE"))

library(vars)
lags=VARselect(Y, lag.max = 30, type = "const")  ## "const", "trend", "both", "none"
nlag=lags$selection[[3]]  ## it selects based on AIC to chnage for example to SIC replace 1 by 3 because: AIC(n)  HQ(n)  SC(n) FPE(n) 
print("nlag")
nlag

### TVP-VAR Estimation
#devtools::install_github("krlmlr/ulimit")
#memory.limit(9999999999)
Sys.setenv('R_MAX_VSIZE'=32000000000)
nlag  = nlag
nfore = 10
#prior = MinnesotaPrior(0.1, k, nlag)
prior = BayesPrior(Y[1:200,], nlag)
tvp_var = TVPVAR(Y, l=c(0.99,0.99), nlag, prior)
B_t = tvp_var$beta_t
Q_t = tvp_var$Q_t


### DYNAMIC CONNECTEDNESS APPROACH
t = nrow(Y)
total = total_corr = matrix(NA, ncol=1, nrow=t)
gfevd = ct = npso = pci = influence = array(NA, c(k, k, t))
net = to = from = matrix(NA, ncol=k, nrow=t)
colnames(gfevd)=rownames(gfevd)=colnames(influence)=rownames(influence)=colnames(pci)=rownames(pci)=NAMES
for (i in 1:t){
  gfevd[,,i] = GFEVD(B_t[,,i], Q_t[,,i], n.ahead=nfore)$GFEVD
  dca = DCA(gfevd[,,i])
  ct[,,i] = dca$CT
  pci[,,i] = dca$PCI
  to[i,] = dca$TO
  from[i,] = dca$FROM
  net[i,] = dca$NET
  npso[,,i] = dca$NPSO
  influence[,,i] = dca$INFLUENCE
  pci[,,i] = dca$PCI
  total[i,] = dca$TCI
  total_corr[i,] = dca$TCI_corrected
  if (i%%100==0) print(paste0(round(100*i/t,2),"%"))
}

options(max.print = .Machine$integer.max)
### Results for retruns 
print(dca$TABLE)  ## for full sample 



#par(mfcol=c(1,1))
#npso.lre = net
#par(mfcol=c(1,1))
#library("qgraph")
#Graph_pcor <- qgraph(corMat, graph = "pcor", layout = "spring")
#Graph_pcor <- qgraph(npso.lre, layout = "spring",
#                     #                     threshold = 0.2,
#                     color = 2:ncol(npso.lre))
##                    edge.color=1:ncol(npso))
##                     color = c("green", "red", "blue", "orange", "pink"))
## apply a hard thresholding to make the major connections more clearly
##con = ifelse(abs(npso.lre) >= ((1/length(npso.lre[,1])) * sum(npso.lre[order(npso.lre, decreasing = T)[1:length(npso.lre[,1])]])), 
##             npso.lre, 0)
##Graph_pcor <- qgraph(con, layout = "spring",
##                     #                     threshold = 0.2,
##                     color = 2:ncol(con))
##

### DYNAMIC TOTAL CONNECTEDNESS
par(mfcol=c(1,1), oma=c(0.5,1.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
plot(date,total_corr, type="l",xaxs="i",col="steelblue4", las=1, main="",ylab="",ylim=c(0,100),yaxs="i",xlab="",tck=-0.02)
grid(NA,NULL)
polygon(c(date,rev(date)),c(c(rep(0,t)),rev(total_corr)),col="steelblue4", border="steelblue4")
lines(date,total,col="red")
box()
legend("topright", c("TCI corrected","TCI"), fill=c("steelblue4", "red"))

### TOTAL DIRECTIONAL CONNECTEDNESS TO OTHERS
split = 4
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,1.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
   plot(date,to[,i], xlab="",ylab="",type="l",xaxs="i",col="steelblue4", las=1, main=paste(NAMES[i],"TO all others"),ylim=c(0,ceiling(max(to))),tck=-0.02,yaxs="i")
   grid(NA,NULL)
   polygon(c(date,rev(date)),c(c(rep(0,t)),rev(to[,i])),col="steelblue4", border="steelblue4")
   box()
}

### TOTAL DIRECTIONAL CONNECTEDNESS FROM OTHERS
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,1.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
   plot(date,from[,i], xlab="",ylab="",type="l",xaxs="i",col="steelblue4", las=1, main=paste(NAMES[i],"FROM all others"),ylim=c(0,100),tck=-0.02,yaxs="i")
   grid(NA,NULL)
   polygon(c(date,rev(date)),c(c(rep(0,t)),rev(from[,i])),col="steelblue4", border="steelblue4")
   box()
}

### NET TOTAL DIRECTIONAL CONNECTEDNESS
split=4
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,1.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
   plot(date,net[,i], xlab="",ylab="",type="l",xaxs="i",col="steelblue4", las=1, main=paste("NET",NAMES[i]),ylim=c(floor(min(net)),ceiling(max(net))),tck=-0.02,yaxs="i")
   grid(NA,NULL)
   polygon(c(date,rev(date)),c(c(rep(0,t)),rev(net[,i])),col="steelblue4", border="steelblue4")
   box()
}



### PAIRWISE CONNECTEDNESS INDEX
#par(mfrow=c(ceiling((k-1)/split),split), oma=c(0.5,1,0.05,0.05), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
par(mfrow = c(ceiling(ncol(npso)/4),6), oma = c(0,1,0,0) + 0.05, mar = c(1,1,1,1) + .05, mgp = c(0, 0.1, 0))
for (i in 1:ncol(Y)){
  for (j in 1:k) {
    if (i!=j) {
      plot(date,pci[j,i,], xlab="",ylab="",type="l",xaxs="i",col="steelblue4", las=1, main=paste0(NAMES[i],"-",NAMES[j]),tck=-0.02,yaxs="i",ylim=c(0,1))
      grid(NA,NULL)
      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(pci[j,i,])),col="steelblue4", border="steelblue4")
      box()
    }
  }
}
### NET PAIRWISE DIRECTIONAL CONNECTEDNESS
#i = 1
#print(NAMES[i])
#par(mfcol=c(ceiling((k-1)/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
#for (j in 1:k) {
#   if (i!=j) {
#      plot(date,npso[j,i,], xlab="",ylab="",type="l",xaxs="i",col="steelblue4", las=1, main=paste0(NAMES[i],"-",NAMES[j]),tck=-0.02,yaxs="i",ylim=c(floor(min(npso)),ceiling(max(npso))))
#      grid(NA,NULL)
#      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(npso[j,i,])),col="steelblue4", border="steelblue4")
#      abline(v=date[IND],lty=2)
#      box()
#   }
#}

### NET PAIRWISE DIRECTIONAL CONNECTEDNESS
#lk = ceiling(sqrt(2))
#split=6
#par(mfcol=c(ceiling(kk/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
kk = k*(k-1)/4
par(mfrow = c(ceiling(ncol(npso)/4),6), oma = c(0,1,0,0) + 0.05, mar = c(1,1,1,1) + .05, mgp = c(0, 0.1, 0))
for (i in 1:k) {
  for (j in 1:k) {
    if (i<j) {
#      plot(date,npso[j,i,], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste0(NAMES[i],"-",NAMES[j]),tck=-0.02,yaxs="i",ylim=c(floor(min(npso)),ceiling(max(npso))))
      plot(date,npso[j,i,], xlab="",ylab="",type="l",xaxs="i",col="steelblue4", las=1, main=paste0(NAMES[i],"-",NAMES[j]),tck=-0.02,yaxs="i",ylim=c(floor(min(npso)),ceiling(max(npso))))
      grid(NA,NULL)
#      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(npso[j,i,])),col="grey20", border="grey20")
      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(npso[j,i,])),col="steelblue4", border="steelblue4")
#      abline(v=date[IND],lty=2)
      box()
    }
  }
}

