
like.prosp.mediator.con.noint<-function(omega,Y,M,A,C){

#Likelihood function for a continuous (normal) mediator assuming no mediator-exposure interaction (gAM=0)

ncov<-ncol(C)

g0<-omega[1] # Disease intercept
gA<-omega[2] # Effect of exposure on disease
gM<-omega[3] # Effect of mediatoron disease
gAM<-0   # Effect of mediator-exposure interaction on disease (assumed 0)
gC<-omega[4:(3+ncov)] # Effect of covariates on disease
b0<-omega[(4+ncov)] # Mediator intercept
bA<-omega[(5+ncov)] # Effect of exposure on mediator
bC<-omega[(6+ncov):(5+ncov+ncov)] # Effect of covariates on mediator
sig2<-omega[6+ncov+ncov] # Residual variance of mediator

Nsub<-length(Y)

#See equation (7) of paper for definition of parameters below

theta_AC_1<-exp(g0+gA*A+(gM+gAM*A)*(b0+bA*A+C%*%bC)+C%*%gC)

theta_AC_2<-exp(0.5*sig2*(gM+gAM*A)**2)
theta_AC<-theta_AC_1*theta_AC_2

#See equation (8) of paper for derivation of likelihood below

l1_fin<-(-Nsub/2)*log(2*pi*sig2)

l2_fin<-sum((-1/(2*sig2))*((M-b0-bA*A-C%*%bC)**2))

l3_fin<-sum(Y*(g0+gA*A+gM*M+gAM*A*M+C%*%gC))

l4_fin<-sum(-log(1+theta_AC))

logl<-l1_fin+l2_fin+l3_fin+l4_fin

return(-logl)

}

grad.prosp.mediator.con.noint<-function(omega,Y,M,A,C){

ncov<-ncol(C)

g0<-omega[1] # Disease intercept
gA<-omega[2] # Effect of exposure on disease
gM<-omega[3] # Effect of mediatoron disease
gAM<-0   # Effect of mediator-exposure interaction on disease (assumed 0)
gC<-omega[4:(3+ncov)] # Effect of covariates on disease
b0<-omega[(4+ncov)] # Mediator intercept
bA<-omega[(5+ncov)] # Effect of exposure on mediator
bC<-omega[(6+ncov):(5+ncov+ncov)] # Effect of covariates on mediator
sig2<-omega[6+ncov+ncov] # Residual variance of mediator


theta_AC_1<-exp(g0+gA*A+(gM+gAM*A)*(b0+bA*A+C%*%bC)+C%*%gC)
theta_AC_2<-exp(0.5*sig2*(gM+gAM*A)**2)
theta_AC<-theta_AC_1*theta_AC_2


theta_fun<-theta_AC/(1+theta_AC)

Nsub<-length(Y)

score<-numeric(length(omega))

score_g0<-sum(Y-theta_fun)
score_gA<-sum(A*(Y-theta_fun))
score_gM<-sum(Y*M-theta_fun*(b0+bA*A+C%*%bC+(gM+gAM*A)*sig2))
score_gAM<-sum(Y*A*M-theta_fun*A*(b0+bA*A+C%*%bC+(gM+gAM*A)*sig2))

cov_fun1<-apply(C,2,function(x){x*(Y-theta_fun)})
score_gC<-apply(cov_fun1,2,sum)

temp_fun<-(M-b0-bA*A-C%*%bC)/sig2

score_b0<-sum(temp_fun-theta_fun*(gM+gAM*A))
score_bA<-sum(A*(temp_fun-theta_fun*(gM+gAM*A)))

cov_fun2<-apply(C,2,function(x){x*(temp_fun-theta_fun*(gM+gAM*A))})

score_bC<-apply(cov_fun2,2,sum)

score_sig2<-sum((-1/(2*sig2))-theta_fun*0.5*((gM+gAM*A)**2)+((temp_fun**2)/(2)))

score[1]<-score_g0
score[2]<-score_gA
score[3]<-score_gM
score[4:(3+ncov)]<-score_gC
score[(4+ncov)]<-score_b0
score[(5+ncov)]<-score_bA
score[(6+ncov):(5+ncov+ncov)]<-score_bC
score[6+ncov+ncov]<-score_sig2

return(-score)

}


like.tvw.mediator.con.noint<-function(omega,Y,M,A,C){

#Implements the VW approach for a continuous mediator with mediator-exposure interaction (gAM)

ncov<-ncol(C)

g0<-omega[1] # Disease intercept
gA<-omega[2] # Effect of exposure on disease
gM<-omega[3] # Effect of mediatoron disease
gAM<-0   # Effect of mediator-exposure interaction on disease (assumed 0)
gC<-omega[4:(3+ncov)] # Effect of covariates on disease
b0<-omega[(4+ncov)] # Mediator intercept
bA<-omega[(5+ncov)] # Effect of exposure on mediator
bC<-omega[(6+ncov):(5+ncov+ncov)] # Effect of covariates on mediator
sig2<-omega[6+ncov+ncov] # Residual variance of mediator

Nsub<-length(Y)

theta_AMC<-exp(g0+gA*A+gM*M+gAM*A*M+C%*%gC)

#Fit the VW approach using only controls to fit the mediator model

l1_fin<-sum(Y*log(theta_AMC))
l2_fin<-sum(-(log(1+theta_AMC)))
l3_fin<-sum((1-Y)*(-1/2)*log(2*pi*sig2))#Fit the mediator model using only controls
l4_fin<-sum((1-Y)*(-1/(2*sig2))*((M-b0-bA*A-C%*%bC)**2))# Fit the mediator model using only controls

logl<-l1_fin+l2_fin+l3_fin+l4_fin

return(-logl)

}

grad.tvw.mediator.con.noint<-function(omega,Y,M,A,C){

ncov<-ncol(C)

g0<-omega[1]
gA<-omega[2]
gM<-omega[3]
gAM<-0
gC<-omega[4:(3+ncov)]
b0<-omega[(4+ncov)]
bA<-omega[(5+ncov)]
bC<-omega[(6+ncov):(5+ncov+ncov)]
sig2<-omega[6+ncov+ncov]

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

temp_fun<-(M-b0-bA*A-C%*%bC)/sig2

score_b0<-sum((1-Y)*temp_fun)
score_bA<-sum((1-Y)*A*temp_fun)

cov_fun2<-apply(C,2,function(x){x*(1-Y)*temp_fun})

score_bC<-apply(cov_fun2,2,sum)

score_sig2<-sum((1-Y)*((-1/(2*sig2))+(1/2)*(temp_fun**2)))

score[1]<-score_g0
score[2]<-score_gA
score[3]<-score_gM
score[4:(3+ncov)]<-score_gC
score[(4+ncov)]<-score_b0
score[(5+ncov)]<-score_bA
score[(6+ncov):(5+ncov+ncov)]<-score_bC
score[6+ncov+ncov]<-score_sig2

return(-score)

}


like.prosp.mediator.con.int<-function(omega,Y,M,A,C){

#Likelihood function for a continuous (normal)  mediator with mediator-exposure interaction (gAM)

ncov<-ncol(C)

g0<-omega[1] # Disease intercept
gA<-omega[2] # Effect of exposure on disease
gM<-omega[3] # Effect of mediator on disease
gAM<-omega[4] # Effect of mediator-exposure interaction on disease
gC<-omega[5:(4+ncov)] # Effect of covariates on disease
b0<-omega[(5+ncov)] # Mediator intercept
bA<-omega[(6+ncov)] # Effect of exposure on mediator
bC<-omega[(7+ncov):(6+ncov+ncov)] # Effect of covariates on mediator
sig2<-omega[7+ncov+ncov] # Residual variance of mediator

Nsub<-length(Y)


#See equation (7) of paper for definition of parameters below

theta_AC_1<-exp(g0+gA*A+(gM+gAM*A)*(b0+bA*A+C%*%bC)+C%*%gC)

theta_AC_2<-exp(0.5*sig2*(gM+gAM*A)**2)
theta_AC<-theta_AC_1*theta_AC_2


#See equation (8) of paper for derivation of likelihood below

l1_fin<-(-Nsub/2)*log(2*pi*sig2)

l2_fin<-sum((-1/(2*sig2))*((M-b0-bA*A-C%*%bC)**2))

l3_fin<-sum(Y*(g0+gA*A+gM*M+gAM*A*M+C%*%gC))

l4_fin<-sum(-log(1+theta_AC))

logl<-l1_fin+l2_fin+l3_fin+l4_fin

return(-logl)

}

grad.prosp.mediator.con.int<-function(omega,Y,M,A,C){

ncov<-ncol(C)

g0<-omega[1] # Disease intercept
gA<-omega[2] # Effect of exposure on disease
gM<-omega[3] # Effect of mediator on disease
gAM<-omega[4] # Effect of mediator-exposure interaction on disease
gC<-omega[5:(4+ncov)] # Effect of covariates on disease
b0<-omega[(5+ncov)] # Mediator intercept
bA<-omega[(6+ncov)] # Effect of exposure on mediator
bC<-omega[(7+ncov):(6+ncov+ncov)] # Effect of covariates on mediator
sig2<-omega[7+ncov+ncov] # Residual variance of mediator

theta_AC_1<-exp(g0+gA*A+(gM+gAM*A)*(b0+bA*A+C%*%bC)+C%*%gC)
theta_AC_2<-exp(0.5*sig2*(gM+gAM*A)**2)
theta_AC<-theta_AC_1*theta_AC_2


theta_fun<-theta_AC/(1+theta_AC)

Nsub<-length(Y)

score<-numeric(length(omega))

score_g0<-sum(Y-theta_fun)
score_gA<-sum(A*(Y-theta_fun))
score_gM<-sum(Y*M-theta_fun*(b0+bA*A+C%*%bC+(gM+gAM*A)*sig2))
score_gAM<-sum(Y*A*M-theta_fun*A*(b0+bA*A+C%*%bC+(gM+gAM*A)*sig2))

cov_fun1<-apply(C,2,function(x){x*(Y-theta_fun)})
score_gC<-apply(cov_fun1,2,sum)

temp_fun<-(M-b0-bA*A-C%*%bC)/sig2

score_b0<-sum(temp_fun-theta_fun*(gM+gAM*A))
score_bA<-sum(A*(temp_fun-theta_fun*(gM+gAM*A)))

cov_fun2<-apply(C,2,function(x){x*(temp_fun-theta_fun*(gM+gAM*A))})

score_bC<-apply(cov_fun2,2,sum)

score_sig2<-sum((-1/(2*sig2))-theta_fun*0.5*((gM+gAM*A)**2)+((temp_fun**2)/(2)))

score[1]<-score_g0
score[2]<-score_gA
score[3]<-score_gM
score[4]<-score_gAM
score[5:(4+ncov)]<-score_gC
score[(5+ncov)]<-score_b0
score[(6+ncov)]<-score_bA
score[(7+ncov):(6+ncov+ncov)]<-score_bC
score[7+ncov+ncov]<-score_sig2

return(-score)

}


like.tvw.mediator.con.int<-function(omega,Y,M,A,C){

#Implements the VW approach for a continuous mediator with mediator-exposure interaction (gAM)

ncov<-ncol(C)

g0<-omega[1] # Disease intercept
gA<-omega[2] # Effect of exposure on disease
gM<-omega[3] # Effect of mediator on disease
gAM<-omega[4] # Effect of mediator-exposure interaction on disease
gC<-omega[5:(4+ncov)] # Effect of covariates on disease
b0<-omega[(5+ncov)] # Mediator intercept
bA<-omega[(6+ncov)] # Effect of exposure on mediator
bC<-omega[(7+ncov):(6+ncov+ncov)] # Effect of covariates on mediator
sig2<-omega[7+ncov+ncov] # Residual variance of mediator

Nsub<-length(Y)

theta_AMC<-exp(g0+gA*A+gM*M+gAM*A*M+C%*%gC)


#Fit the VW approach using only controls to fit the mediator model


l1_fin<-sum(Y*log(theta_AMC))
l2_fin<-sum(-(log(1+theta_AMC)))
l3_fin<-sum((1-Y)*(-1/2)*log(2*pi*sig2)) #Fit the mediator model using only controls
l4_fin<-sum((1-Y)*(-1/(2*sig2))*((M-b0-bA*A-C%*%bC)**2)) # Fit the mediator model using only controls

logl<-l1_fin+l2_fin+l3_fin+l4_fin

return(-logl)

}

grad.tvw.mediator.con.int<-function(omega,Y,M,A,C){

ncov<-ncol(C)

g0<-omega[1] # Disease intercept
gA<-omega[2] # Effect of exposure on disease
gM<-omega[3] # Effect of mediator on disease
gAM<-omega[4] # Effect of mediator-exposure interaction on disease
gC<-omega[5:(4+ncov)] # Effect of covariates on disease
b0<-omega[(5+ncov)] # Mediator intercept
bA<-omega[(6+ncov)] # Effect of exposure on mediator
bC<-omega[(7+ncov):(6+ncov+ncov)] # Effect of covariates on mediator
sig2<-omega[7+ncov+ncov] # Residual variance of mediator

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

temp_fun<-(M-b0-bA*A-C%*%bC)/sig2

score_b0<-sum((1-Y)*temp_fun)
score_bA<-sum((1-Y)*A*temp_fun)

cov_fun2<-apply(C,2,function(x){x*(1-Y)*temp_fun})
score_bC<-apply(cov_fun2,2,sum)

score_sig2<-sum((1-Y)*((-1/(2*sig2))+(1/2)*(temp_fun**2)))

score[1]<-score_g0
score[2]<-score_gA
score[3]<-score_gM
score[4]<-score_gAM
score[5:(4+ncov)]<-score_gC
score[(5+ncov)]<-score_b0
score[(6+ncov)]<-score_bA
score[(7+ncov):(6+ncov+ncov)]<-score_bC
score[7+ncov+ncov]<-score_sig2

return(-score)

}
