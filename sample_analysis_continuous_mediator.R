#############
#This code implements the prospective likelihood approach for mediation analysis
#outlined in the Satten et al. manuscript


source('mediation_functions_binary_mediator.R') #Likelihood and gradient functions assuming a binary mediator with or without mediator-exposure interaction

source('mediation_functions_continuous_mediator.R') #Likelihood and gradient functions assuming a continuous mediator with or without mediator-exposure interaction

source('effect_functions.R') #Functions for calculating estimates (and standard errors) of total effects, direct effects, and indirect effects


###############
#Read in sample data file
# Column 1: Disease indicator (1=case, 0=control)
# Column 2: Continuous exposure
# Column 3: Continuous mediator
# Column 4: Binary covariate
# Column 5: Continuous covariate

datafile<-read.table('sample_continuous_mediator.dat',header=F)
colnames(datafile)<-c('Phenotype','Exposure','Mediator','Cov_bin','Cov_con')

Ydat<-datafile$Phenotype 
Adat<-datafile$Exposure
Mdat<-datafile$Mediator
Cdat<-data.frame(datafile$Cov_bin,datafile$Cov_con)
Cdat<-as.matrix(Cdat)

ncov<-ncol(Cdat) # Number of covariates

################
#Fit the likelihood shown in equation (8) of manuscript
#Likelihood function shown in source files

init.val.con.noint<-c(runif(5+2*ncov,-2,2),0.3) # Initial parameter values
lower.val.con.noint<-c(rep(-Inf,5+2*ncov),0.01) # Lower bounds of parameter estimates. Last element in list corresponds to the variance of the continuous mediator 
upper.val.con.noint<-c(rep(Inf,5+2*ncov),Inf) # Upper bounds of parameter estimates

#Optimize the likelihood
fit.prosp.con.noint<-optim(init.val.con.noint,method="L-BFGS-B",fn=like.prosp.mediator.con.noint,gr=grad.prosp.mediator.con.noint,Y=Ydat,M=Mdat,A=Adat,C=Cdat,hessian=T,lower=lower.val.con.noint,upper=upper.val.con.noint,control = list(maxit = 10000))

#Extract the parameter estimates and standard errors from likelihood

parest_prosp_noint<-fit.prosp.con.noint$par # Parameter estimates 

fisher_info_prosp_noint<-solve(fit.prosp.con.noint$hessian)
se_prosp_noint<-sqrt(diag(fisher_info_prosp_noint)) # Standard errors of estimates


#Calculate total effect, direct effect, and indirect effect and standard errors 
#Effect size functions shown in source files

effects_prosp_noint<-effect_con_noint(fit.prosp.con.noint$par,fisher_info_prosp_noint,1,0,0,rep(0,ncov))



###########################
#Now, fit the likelihood in (8) assuming a mediator-exposure interaction effect
#Likelihood function shown in source files


init.val.con.int<-c(runif(6+2*ncov,-2,2),0.3) # Initial parameter values
lower.val.con.int<-c(rep(-Inf,6+2*ncov),0.01) # Lower bound of parameter values. Last element in list corresponds to the variance of the continuous mediator
upper.val.con.int<-c(rep(Inf,6+2*ncov),Inf)  # Upper bound of parameter values

#Optimize the likelihood
fit.prosp.con.int<-optim(init.val.con.int,method="L-BFGS-B",fn=like.prosp.mediator.con.int,gr=grad.prosp.mediator.con.int,Y=Ydat,M=Mdat,A=Adat,C=Cdat,hessian=T,lower=lower.val.con.int,upper=upper.val.con.int,control = list(maxit = 10000)) # Optimize the likelihood

#Extract the parameter estimates and standard errors from likelihood

parest_prosp_int<-fit.prosp.con.int$par   #Parameter estimates

fisher_info_prosp_int<-solve(fit.prosp.con.int$hessian)
se_prosp_int<-sqrt(diag(fisher_info_prosp_int)) # Standard errors

#Calculate total effect, direct effect, and indirect effect and standard errors 
#Effect size functions shown in source files

effects_prosp_int<-effect_con_int(fit.prosp.con.int$par,fisher_info_prosp_int,1,0,0,rep(0,ncov))






