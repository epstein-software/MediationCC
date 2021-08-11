
like.prosp.mediator.bin.noint<-function(omega,Y,M,A,C){

#Likelihood function for a binary mediator with no interaction (gAM=0)

ncov<-ncol(C) # Number of covariates

#See equation (1) of paper for definition of parameters below

g0<-omega[1] # Disease intercept
gA<-omega[2] # Effect of exposure on disease
gM<-omega[3] # Effect of mediator on disease
gAM<-0
gC<-omega[4:(3+ncov)] # Effect of covariates on disease
b0<-omega[(4+ncov)] # Mediator intercept
bA<-omega[(5+ncov)] # Effect of exposure on mediator
bC<-omega[(6+ncov):(5+ncov+ncov)] # Effect of covariates on mediator


Nsub<-length(Y) # Total number of subjects


#See equation (7) of paper for definition of parameters below

theta_AC_1<-exp(g0+gA*A+C%*%gC)/(1+exp(b0+bA*A+C%*%bC))

theta_AC_2<-exp(g0+gA*A+gM+gAM*A+C%*%gC)*(exp(b0+bA*A+C%*%bC)/(1+exp(b0+bA*A+C%*%bC)))

theta_AC<-theta_AC_1+theta_AC_2


#See equation (8) of paper for derivation of likelihood below

l1_fin<-sum(M*(b0+bA*A+C%*%bC))
l2_fin<-sum(-log(1+exp(b0+bA*A+C%*%bC)))

l3_fin<-sum(Y*(g0+gA*A+gM*M+gAM*A*M+C%*%gC))

l4_fin<-sum(-log(1+theta_AC))

logl<-l1_fin+l2_fin+l3_fin+l4_fin

return(-logl)

}

grad.prosp.mediator.bin.noint<-function(omega,Y,M,A,C){

ncov<-ncol(C)

g0<-omega[1] # Disease intercept
gA<-omega[2] # Effect of exposure on disease
gM<-omega[3] # Effect of mediator on disease
gAM<-0
gC<-omega[4:(3+ncov)] # Effect of covariates on disease
b0<-omega[(4+ncov)] # Mediator intercept
bA<-omega[(5+ncov)] # Effect of exposure on mediator
bC<-omega[(6+ncov):(5+ncov+ncov)] # Effect of covariates on mediator


theta_AC_1<-exp(g0+gA*A+C%*%gC)/(1+exp(b0+bA*A+C%*%bC))

theta_AC_2<-exp(g0+gA*A+gM+gAM*A+C%*%gC)*(exp(b0+bA*A+C%*%bC)/(1+exp(b0+bA*A+C%*%bC)))

theta_AC<-theta_AC_1+theta_AC_2

theta_fun<-theta_AC/(1+theta_AC)

Nsub<-length(Y)

# Derive the derivatives of the likelihood

score<-numeric(length(omega))

score_g0<-sum(Y-theta_fun)
score_gA<-sum(A*(Y-theta_fun))
score_gM<-sum(Y*M-(theta_AC_2/(1+theta_AC)))
score_gAM<-sum(Y*A*M-(A*theta_AC_2/(1+theta_AC)))

cov_fun1<-apply(C,2,function(x){x*(Y-theta_fun)})
score_gC<-apply(cov_fun1,2,sum)

pi<-exp(b0+bA*A+C%*%bC)/(1+exp(b0+bA*A+C%*%bC))
c1<-exp(g0+gA*A+C%*%gC)
c2<-exp(g0+gA*A+gM+gAM*A+C%*%gC)

tempfun_num<--c1*pi*(1-pi)+c2*pi*(1-pi)
tempfun_den<-1+theta_AC

tempfun<-tempfun_num/tempfun_den

score_b0<-sum(M-pi-tempfun)
score_bA<-sum(A*(M-pi-tempfun))

cov_fun2<-apply(C,2,function(x){x*(M-pi-tempfun)})
score_bC<-apply(cov_fun2,2,sum)

score[1]<-score_g0
score[2]<-score_gA
score[3]<-score_gM
score[4:(3+ncov)]<-score_gC
score[(4+ncov)]<-score_b0
score[(5+ncov)]<-score_bA
score[(6+ncov):(5+ncov+ncov)]<-score_bC

return(-score)

}



like.tvw.mediator.bin.noint<-function(omega,Y,M,A,C){

#Implements the VW approach for a binary mediator with no interaction (gAM=0)

ncov<-ncol(C)

g0<-omega[1] # Disease intercept
gA<-omega[2] # Effect of exposure on disease
gM<-omega[3] # Effect of mediator on disease
gAM<-0
gC<-omega[4:(3+ncov)] # Effect of covariates on disease
b0<-omega[(4+ncov)] # Mediator intercept
bA<-omega[(5+ncov)] # Effect of exposure on mediator
bC<-omega[(6+ncov):(5+ncov+ncov)] # Effect of covariates on mediator

Nsub<-length(Y)

#See equation (6) for derivation of below

theta_AMC<-exp(g0+gA*A+gM*M+gAM*A*M+C%*%gC)

#Fit the VW approach using only controls to fit the mediator model

l1_fin<-sum(Y*log(theta_AMC))
l2_fin<-sum(-(log(1+theta_AMC)))
l3_fin<-sum((1-Y)*(M)*(b0+bA*A+C%*%bC)) # Fit the mediator model using controls only
l4_fin<-sum((1-Y)*(-log(1+exp(b0+bA*A+C%*%bC)))) # Fit the mediator model using controls only

logl<-l1_fin+l2_fin+l3_fin+l4_fin

return(-logl)

}

grad.tvw.mediator.bin.noint<-function(omega,Y,M,A,C){

#Derives the score functions for the VW approach for a binary mediator with no interaction (gAM=0)


ncov<-ncol(C)

g0<-omega[1] # Disease intercept
gA<-omega[2] # Effect of exposure on disease
gM<-omega[3] # Effect of mediator on disease
gAM<-0
gC<-omega[4:(3+ncov)] # Effect of covariates on disease
b0<-omega[(4+ncov)] # Mediator intercept
bA<-omega[(5+ncov)] # Effect of exposure on mediator
bC<-omega[(6+ncov):(5+ncov+ncov)] # Effect of covariates on mediator


#See equation (6) for derivation of below


theta_AMC<-exp(g0+gA*A+gM*M+gAM*A*M+C%*%gC)

theta_fun<-theta_AMC/(1+theta_AMC)

Nsub<-length(Y)

score<-numeric(length(omega))

score_g0<-sum(Y-theta_fun)
score_gA<-sum(A*(Y-theta_fun))
score_gM<-sum(M*(Y-theta_fun))
score_gAM<-sum(A*M*(Y-theta_fun))

cov_fun1<-apply(C,2,function(x){x*(Y-theta_fun)})
score_gC<-apply(cov_fun1,2,sum)

theta2_AC<-exp(b0+bA*A+C%*%bC)

theta2_fun<-theta2_AC/(1+theta2_AC)

score_b0<-sum((1-Y)*(M-theta2_fun))
score_bA<-sum((1-Y)*A*(M-theta2_fun))

cov_fun2<-apply(C,2,function(x){x*(1-Y)*(M-theta2_fun)})
score_bC<-apply(cov_fun2,2,sum)

score[1]<-score_g0
score[2]<-score_gA
score[3]<-score_gM
score[4:(3+ncov)]<-score_gC
score[(4+ncov)]<-score_b0
score[(5+ncov)]<-score_bA
score[(6+ncov):(5+ncov+ncov)]<-score_bC

return(-score)

}



like.prosp.mediator.bin.int<-function(omega,Y,M,A,C){

#Likelihood function for a binary mediator with mediator-exposure interaction (gAM)

ncov<-ncol(C)

g0<-omega[1] # Disease intercept
gA<-omega[2] # Effect of exposure on disease
gM<-omega[3] # Effect of mediator on disease
gAM<-omega[4] # Effect of mediator-exposure interaction on disease
gC<-omega[5:(4+ncov)] # Effect of covariates on disease
b0<-omega[(5+ncov)] # Mediator intercept
bA<-omega[(6+ncov)] # Effect of exposure on mediator
bC<-omega[(7+ncov):(6+ncov+ncov)] # Effect of covariates on mediator


Nsub<-length(Y)

#See equation (7) of paper for definition of parameters below

theta_AC_1<-exp(g0+gA*A+C%*%gC)/(1+exp(b0+bA*A+C%*%bC))

theta_AC_2<-exp(g0+gA*A+gM+gAM*A+C%*%gC)*(exp(b0+bA*A+C%*%bC)/(1+exp(b0+bA*A+C%*%bC)))

theta_AC<-theta_AC_1+theta_AC_2


#See equation (8) of paper for derivation of likelihood below

l1_fin<-sum(M*(b0+bA*A+C%*%bC))
l2_fin<-sum(-log(1+exp(b0+bA*A+C%*%bC)))

l3_fin<-sum(Y*(g0+gA*A+gM*M+gAM*A*M+C%*%gC))

l4_fin<-sum(-log(1+theta_AC))

logl<-l1_fin+l2_fin+l3_fin+l4_fin

return(-logl)

}

grad.prosp.mediator.bin.int<-function(omega,Y,M,A,C){

ncov<-ncol(C)

g0<-omega[1] # Disease intercept
gA<-omega[2] # Effect of exposure on disease
gM<-omega[3] # Effect of mediator on disease
gAM<-omega[4] # Effect of mediator-exposure interaction on disease
gC<-omega[5:(4+ncov)] # Effect of covariates on disease
b0<-omega[(5+ncov)] # Mediator intercept
bA<-omega[(6+ncov)] # Effect of exposure on mediator
bC<-omega[(7+ncov):(6+ncov+ncov)] # Effect of covariates on mediator


# Derive the derivatives of the likelihood

theta_AC_1<-exp(g0+gA*A+C%*%gC)/(1+exp(b0+bA*A+C%*%bC))

theta_AC_2<-exp(g0+gA*A+gM+gAM*A+C%*%gC)*(exp(b0+bA*A+C%*%bC)/(1+exp(b0+bA*A+C%*%bC)))

theta_AC<-theta_AC_1+theta_AC_2

theta_fun<-theta_AC/(1+theta_AC)

Nsub<-length(Y)

score<-numeric(length(omega))

score_g0<-sum(Y-theta_fun)
score_gA<-sum(A*(Y-theta_fun))
score_gM<-sum(Y*M-(theta_AC_2/(1+theta_AC)))
score_gAM<-sum(Y*A*M-(A*theta_AC_2/(1+theta_AC)))

cov_fun1<-apply(C,2,function(x){x*(Y-theta_fun)})
score_gC<-apply(cov_fun1,2,sum)

#score_gC<-sum(C*(Y-theta_fun))

pi<-exp(b0+bA*A+C%*%bC)/(1+exp(b0+bA*A+C%*%bC))
c1<-exp(g0+gA*A+C%*%gC)
c2<-exp(g0+gA*A+gM+gAM*A+C%*%gC)

tempfun_num<--c1*pi*(1-pi)+c2*pi*(1-pi)
tempfun_den<-1+theta_AC

tempfun<-tempfun_num/tempfun_den

score_b0<-sum(M-pi-tempfun)
score_bA<-sum(A*(M-pi-tempfun))

cov_fun2<-apply(C,2,function(x){x*(M-pi-tempfun)})
score_bC<-apply(cov_fun2,2,sum)


score[1]<-score_g0
score[2]<-score_gA
score[3]<-score_gM
score[4]<-score_gAM
score[5:(4+ncov)]<-score_gC
score[(5+ncov)]<-score_b0
score[(6+ncov)]<-score_bA
score[(7+ncov):(6+ncov+ncov)]<-score_bC

return(-score)

}



like.tvw.mediator.bin.int<-function(omega,Y,M,A,C){

#Implements the VW approach for a binary mediator with mediator-exposure interaction (gAM)

ncov<-ncol(C)

g0<-omega[1] # Disease intercept
gA<-omega[2] # Effect of exposure on disease
gM<-omega[3] # Effect of mediator on disease
gAM<-omega[4] # Effect of mediator-exposure interaction on disease
gC<-omega[5:(4+ncov)] # Effect of covariates on disease
b0<-omega[(5+ncov)] # Mediator intercept
bA<-omega[(6+ncov)] # Effect of exposure on mediator
bC<-omega[(7+ncov):(6+ncov+ncov)] # Effect of covariates on mediator

Nsub<-length(Y)

theta_AMC<-exp(g0+gA*A+gM*M+gAM*A*M+C%*%gC)


#Fit the VW approach using only controls to fit the mediator model

l1_fin<-sum(Y*log(theta_AMC)) 
l2_fin<-sum(-(log(1+theta_AMC)))
l3_fin<-sum((1-Y)*(M)*(b0+bA*A+C%*%bC)) # Fit the mediator model using only controls
l4_fin<-sum((1-Y)*(-log(1+exp(b0+bA*A+C%*%bC)))) # Fit the mediator model using only controls

logl<-l1_fin+l2_fin+l3_fin+l4_fin

return(-logl)

}

grad.tvw.mediator.bin.int<-function(omega,Y,M,A,C){

ncov<-ncol(C)

g0<-omega[1] # Disease intercept
gA<-omega[2] # Effect of exposure on disease
gM<-omega[3] # Effect of mediator on disease
gAM<-omega[4] # Effect of mediator-exposure interaction on disease
gC<-omega[5:(4+ncov)] # Effect of covariates on disease
b0<-omega[(5+ncov)] # Mediator intercept
bA<-omega[(6+ncov)] # Effect of exposure on mediator
bC<-omega[(7+ncov):(6+ncov+ncov)] # Effect of covariates on mediator

theta_AMC<-exp(g0+gA*A+gM*M+gAM*A*M+C%*%gC)

theta_fun<-theta_AMC/(1+theta_AMC)

Nsub<-length(Y)

score<-numeric(length(omega))

score_g0<-sum(Y-theta_fun)
score_gA<-sum(A*(Y-theta_fun))
score_gM<-sum(M*(Y-theta_fun))
score_gAM<-sum(A*M*(Y-theta_fun))

cov_fun1<-apply(C,2,function(x){x*(Y-theta_fun)})
score_gC<-apply(cov_fun1,2,sum)

theta2_AC<-exp(b0+bA*A+C%*%bC)
theta2_fun<-theta2_AC/(1+theta2_AC)

score_b0<-sum((1-Y)*(M-theta2_fun))
score_bA<-sum((1-Y)*A*(M-theta2_fun))

cov_fun2<-apply(C,2,function(x){x*(1-Y)*(M-theta2_fun)})
score_bC<-apply(cov_fun2,2,sum)


score[1]<-score_g0
score[2]<-score_gA
score[3]<-score_gM
score[4]<-score_gAM
score[5:(4+ncov)]<-score_gC
score[(5+ncov)]<-score_b0
score[(6+ncov)]<-score_bA
score[(7+ncov):(6+ncov+ncov)]<-score_bC

return(-score)

}





