
effect_con_int<-function(omega,sigmat,A,Astar,M,C){

#Function assumes continuous mediator and that mediator-exposure interaction parameter (gAM) is included

ncov<-length(C)

g0<-omega[1]
gA<-omega[2]
gM<-omega[3]
gAM<-omega[4]
gC<-omega[5:(4+ncov)]
b0<-omega[(5+ncov)]
bA<-omega[(6+ncov)]
bC<-omega[(7+ncov):(6+ncov+ncov)]
sig2<-omega[7+ncov+ncov]

OR_CDE<-exp(gA+gAM*M)*((A-Astar)) #Conditional direct effect of A on Y

OR_NDE<-exp((gA+gAM*(b0+bA*Astar+C%*%bC+gM*sig2))*(A-Astar)+0.5*(gAM**2)*sig2*(A**2-Astar**2)) # Natural direct effect of A of Y

OR_NIE<-exp((gM*bA+gAM*bA*A)*(A-Astar))  # Natural indirect effect of A on Y through M

OR_TE<-OR_NDE*OR_NIE  # Total effect of A on Y

####################################
#Use delta method to construct standard errors of log(OR_CDE), log(OR_NDE), log(OR_NIE), and log(OR_NTE)

#Formula below implemented from appendix of Valeri and VanderWeele 

gamma_logCDE<-c(0,1,0,M,rep(0,ncov),0,0,rep(0,ncov),0)

gamma_logNDE<-c(0,1,gAM*sig2,b0+bA*Astar+C%*%bC+gM*sig2+gAM*sig2*(A+Astar),rep(0,ncov),gAM,gAM*Astar,gAM*C,gM*gAM+0.5*(gAM**2)*(A+Astar))

gamma_logNIE<-c(0,0,bA,bA*A,rep(0,ncov),0,gM+gAM*A,rep(0,ncov),0)

gamma_logTE<-c(0,1,gAM*sig2+bA,b0+bA*(A+Astar)+C%*%bC+gM*sig2+gAM*sig2*(A**2-Astar**2),rep(0,ncov),gAM,gAM*(A+Astar)+gM,gAM*C,0.5*(gAM**2)*(A**2-Astar**2))


var_logCDE<-t(gamma_logCDE)%*%sigmat%*%gamma_logCDE
var_logNDE<-t(gamma_logNDE)%*%sigmat%*%gamma_logNDE
var_logNIE<-t(gamma_logNIE)%*%sigmat%*%gamma_logNIE
var_logTE<-t(gamma_logTE)%*%sigmat%*%gamma_logTE

se_logCDE<-sqrt(var_logCDE)*abs(A-Astar)
se_logNDE<-sqrt(var_logNDE)*abs(A-Astar)
se_logNIE<-sqrt(var_logNIE)*abs(A-Astar)
se_logTE<-sqrt(var_logTE)*abs(A-Astar)


#################################


results<-numeric(8)

results[1]<-OR_TE
results[2]<-se_logTE
results[3]<-OR_NDE
results[4]<-se_logNDE
results[5]<-OR_NIE
results[6]<-se_logNIE
results[7]<-OR_CDE
results[8]<-se_logCDE

return(results)

}


effect_con_noint<-function(omega,sigmat,A,Astar,M,C){

#Function assumes continuous mediator and that mediator-exposure interaction parameter (gAM) is excluded (gAM=0)

ncov<-length(C)

g0<-omega[1]
gA<-omega[2]
gM<-omega[3]
gAM<-0 
gC<-omega[4:(3+ncov)]
b0<-omega[(4+ncov)]
bA<-omega[(5+ncov)]
bC<-omega[(6+ncov):(5+ncov+ncov)]
sig2<-omega[6+ncov+ncov]

OR_CDE<-exp(gA+gAM*M)*((A-Astar)) #Conditional direct effect of A on Y
OR_NDE<-exp((gA+gAM*(b0+bA*Astar+C%*%bC+gM*sig2))*(A-Astar)+0.5*(gAM**2)*sig2*(A**2-Astar**2))  # Natural direct effect of A of Y

OR_NIE<-exp((gM*bA+gAM*bA*A)*(A-Astar)) # Natural indirect effect of A on Y through M

OR_TE<-OR_NDE*OR_NIE # Total effect of A on Y

####################################
#Use delta method to construct standard errors of log(OR_CDE), log(OR_NDE), log(OR_NIE), and log(OR_NTE)

#Formula below implemented from appendix of Valeri and VanderWeele 


gamma_logCDE<-c(0,1,0,M,rep(0,ncov),0,0,rep(0,ncov),0)
gamma_logCDE<-gamma_logCDE[-c(4)]

gamma_logNDE<-c(0,1,gAM*sig2,b0+bA*Astar+C%*%bC+gM*sig2+gAM*sig2*(A+Astar),rep(0,ncov),gAM,gAM*Astar,gAM*C,gM*gAM+0.5*(gAM**2)*(A+Astar))
gamma_logNDE<-gamma_logNDE[-c(4)]

gamma_logNIE<-c(0,0,bA,bA*A,rep(0,ncov),0,gM+gAM*A,rep(0,ncov),0)
gamma_logNIE<-gamma_logNIE[-c(4)]


gamma_logTE<-c(0,1,gAM*sig2+bA,b0+bA*(A+Astar)+C%*%bC+gM*sig2+gAM*sig2*(A**2-Astar**2),rep(0,ncov),gAM,gAM*(A+Astar)+gM,gAM*C,0.5*(gAM**2)*(A**2-Astar**2))

gamma_logTE<-gamma_logTE[-c(4)]


var_logCDE<-t(gamma_logCDE)%*%sigmat%*%gamma_logCDE
var_logNDE<-t(gamma_logNDE)%*%sigmat%*%gamma_logNDE
var_logNIE<-t(gamma_logNIE)%*%sigmat%*%gamma_logNIE
var_logTE<-t(gamma_logTE)%*%sigmat%*%gamma_logTE

se_logCDE<-sqrt(var_logCDE)*abs(A-Astar)
se_logNDE<-sqrt(var_logNDE)*abs(A-Astar)
se_logNIE<-sqrt(var_logNIE)*abs(A-Astar)
se_logTE<-sqrt(var_logTE)*abs(A-Astar)

#################################3


results<-numeric(8)

results[1]<-OR_TE
results[2]<-se_logTE
results[3]<-OR_NDE
results[4]<-se_logNDE
results[5]<-OR_NIE
results[6]<-se_logNIE
results[7]<-OR_CDE
results[8]<-se_logCDE

return(results)

}





effect_bin_noint<-function(omega,sigmat,A,Astar,M,C){

#Function assumes binary mediator and that mediator-exposure interaction parameter (gAM) is excluded (gAM=0)

ncov<-length(C)

g0<-omega[1]
gA<-omega[2]
gM<-omega[3]
gAM<-0
gC<-omega[4:(3+ncov)]
b0<-omega[(4+ncov)]
bA<-omega[(5+ncov)]
bC<-omega[(6+ncov):(5+ncov+ncov)]

OR_CDE<-exp(gA+gAM*M)*((A-Astar)) #Controlled direct effect of A on Y

OR_NDE_num<-exp(gA*A)*(1+exp(gM+gAM*A+b0+bA*Astar+C%*%bC))
OR_NDE_den<-exp(gA*Astar)*(1+exp(gM+gAM*Astar+b0+bA*Astar+C%*%bC))

OR_NDE<-OR_NDE_num/OR_NDE_den # Natural direct effect of A on Y

OR_NIE_num<-(1+exp(b0+bA*Astar+C%*%bC))*(1+exp(gM+gAM*A+b0+bA*A+C%*%bC))
OR_NIE_den<-(1+exp(b0+bA*A+C%*%bC))*(1+exp(gM+gAM*A+b0+bA*Astar+C%*%bC))

OR_NIE<-OR_NIE_num/OR_NIE_den # Natural indirect effect of A on Y through M

OR_TE<-OR_NDE*OR_NIE # Total effect of A on Y



####################################
#Use delta method to construct standard errors of log(OR_CDE), log(OR_NDE), log(OR_NIE), and log(OR_NTE)

#Formula below implemented from appendix of Valeri and VanderWeele 


#gamma_logCDE<-c(0,(A-Astar),rep(0,ncov),0,0,0,rep(0,ncov))
gamma_logCDE<-c(0,(A-Astar),0,rep(0,ncov),0,0,rep(0,ncov))



Afun_num<-exp(gM+gAM*A+b0+bA*Astar+C%*%bC)
Afun_den<-1+exp(gM+gAM*A+b0+bA*Astar+C%*%bC)

Afun<-Afun_num/Afun_den

Bfun_num<-exp(gM+gAM*Astar+b0+bA*Astar+C%*%bC)
Bfun_den<-1+exp(gM+gAM*Astar+b0+bA*Astar+C%*%bC)

Bfun<-Bfun_num/Bfun_den


d1_nde<-Afun-Bfun
d2_nde<-Astar*(Afun-Bfun)
d3_nde<-as.vector(C)*(Afun-Bfun)
d4_nde<-0
d5_nde<-(A-Astar)
d6_nde<-Afun-Bfun
d7_nde<-A*Afun-Astar*Bfun
d8_nde<-rep(0,ncov)

gamma_logNDE<-c(d4_nde,d5_nde,d6_nde,d8_nde,d1_nde,d2_nde,d3_nde)


Afun_num<-exp(gM+gAM*A+b0+bA*A+C%*%bC)
Afun_den<-1+exp(gM+gAM*A+b0+bA*A+C%*%bC)

Afun<-Afun_num/Afun_den

Bfun_num<-exp(gM+gAM*A+b0+bA*Astar+C%*%bC)
Bfun_den<-1+exp(gM+gAM*A+b0+bA*Astar+C%*%bC)

Bfun<-Bfun_num/Bfun_den

Kfun_num<-exp(b0+bA*A+C%*%bC)
Kfun_den<-1+exp(b0+bA*A+C%*%bC)

Kfun<-Kfun_num/Kfun_den

Dfun_num<-exp(b0+bA*Astar+C%*%bC)
Dfun_den<-1+exp(b0+bA*Astar+C%*%bC)

Dfun<-Dfun_num/Dfun_den


d1_nie<-(Dfun+Afun)-(Kfun+Bfun)
d2_nie<-Astar*(Dfun-Bfun)+A*(Afun-Kfun)
d3_nie<-as.vector(C)*((Dfun+Afun)-(Kfun+Bfun))
d4_nie<-0
d5_nie<-0
d6_nie<-Afun-Bfun
d7_nie<-A*(Afun-Bfun)
d8_nie<-rep(0,ncov)

gamma_logNIE<-c(d4_nie,d5_nie,d6_nie,d8_nie,d1_nie,d2_nie,d3_nie)

d1_te<-d1_nde+d1_nie
d2_te<-d2_nde+d2_nie
d3_te<-d3_nde+d3_nie
d4_te<-d4_nde+d4_nie
d5_te<-d5_nde+d5_nie
d6_te<-d6_nde+d6_nie
d7_te<-d7_nde+d7_nie
d8_te<-d8_nde+d8_nie

gamma_logTE<-c(d4_te,d5_te,d6_te,d8_te,d1_te,d2_te,d3_te)

var_logCDE<-t(gamma_logCDE)%*%sigmat%*%gamma_logCDE
var_logNDE<-t(gamma_logNDE)%*%sigmat%*%gamma_logNDE
var_logNIE<-t(gamma_logNIE)%*%sigmat%*%gamma_logNIE
var_logTE<-t(gamma_logTE)%*%sigmat%*%gamma_logTE

se_logCDE<-sqrt(var_logCDE)
se_logNDE<-sqrt(var_logNDE)
se_logNIE<-sqrt(var_logNIE)
se_logTE<-sqrt(var_logTE)

############################


results<-numeric(8)

results[1]<-OR_TE
results[2]<-se_logTE
results[3]<-OR_NDE
results[4]<-se_logNDE
results[5]<-OR_NIE
results[6]<-se_logNIE
results[7]<-OR_CDE
results[8]<-se_logCDE

return(results)

}


effect_bin_int<-function(omega,sigmat,A,Astar,M,C){

#Function assumes binary mediator and that mediator-exposure interaction parameter (gAM) is included

ncov<-length(C)

g0<-omega[1]
gA<-omega[2]
gM<-omega[3]
gAM<-omega[4]
gC<-omega[5:(4+ncov)]
b0<-omega[(5+ncov)]
bA<-omega[(6+ncov)]
bC<-omega[(7+ncov):(6+ncov+ncov)]

OR_CDE<-exp(gA+gAM*M)*((A-Astar)) #Controlled direct effect of A on Y

OR_NDE_num<-exp(gA*A)*(1+exp(gM+gAM*A+b0+bA*Astar+C%*%bC))
OR_NDE_den<-exp(gA*Astar)*(1+exp(gM+gAM*Astar+b0+bA*Astar+C%*%bC))

OR_NDE<-OR_NDE_num/OR_NDE_den # Natural direct effect of A on Y

OR_NIE_num<-(1+exp(b0+bA*Astar+C%*%bC))*(1+exp(gM+gAM*A+b0+bA*A+C%*%bC))
OR_NIE_den<-(1+exp(b0+bA*A+C%*%bC))*(1+exp(gM+gAM*A+b0+bA*Astar+C%*%bC))

OR_NIE<-OR_NIE_num/OR_NIE_den # Natural indirect effect of A on Y through M

OR_TE<-OR_NDE*OR_NIE # Total effect of A on Y



####################################
#Use delta method to construct standard errors of log(OR_CDE), log(OR_NDE), log(OR_NIE), and log(OR_NTE)

#Formula below implemented from appendix of Valeri and VanderWeele 


#gamma_logCDE<-c(0,(A-Astar),rep(0,ncov),0,0,0,rep(0,ncov))
gamma_logCDE<-c(0,(A-Astar),0,M*(A-Astar),rep(0,ncov),0,0,rep(0,ncov))



Afun_num<-exp(gM+gAM*A+b0+bA*Astar+C%*%bC)
Afun_den<-1+exp(gM+gAM*A+b0+bA*Astar+C%*%bC)

Afun<-Afun_num/Afun_den

Bfun_num<-exp(gM+gAM*Astar+b0+bA*Astar+C%*%bC)
Bfun_den<-1+exp(gM+gAM*Astar+b0+bA*Astar+C%*%bC)

Bfun<-Bfun_num/Bfun_den


d1_nde<-Afun-Bfun
d2_nde<-Astar*(Afun-Bfun)
d3_nde<-as.vector(C)*(Afun-Bfun)
d4_nde<-0
d5_nde<-(A-Astar)
d6_nde<-Afun-Bfun
d7_nde<-A*Afun-Astar*Bfun
d8_nde<-rep(0,ncov)


gamma_logNDE<-c(d4_nde,d5_nde,d6_nde,d7_nde,d8_nde,d1_nde,d2_nde,d3_nde)

#gamma_logNDE<-c(d4_nde,d5_nde,d6_nde,d8_nde,d1_nde,d2_nde,d3_nde)


Afun_num<-exp(gM+gAM*A+b0+bA*A+C%*%bC)
Afun_den<-1+exp(gM+gAM*A+b0+bA*A+C%*%bC)

Afun<-Afun_num/Afun_den

Bfun_num<-exp(gM+gAM*A+b0+bA*Astar+C%*%bC)
Bfun_den<-1+exp(gM+gAM*A+b0+bA*Astar+C%*%bC)

Bfun<-Bfun_num/Bfun_den

Kfun_num<-exp(b0+bA*A+C%*%bC)
Kfun_den<-1+exp(b0+bA*A+C%*%bC)

Kfun<-Kfun_num/Kfun_den

Dfun_num<-exp(b0+bA*Astar+C%*%bC)
Dfun_den<-1+exp(b0+bA*Astar+C%*%bC)

Dfun<-Dfun_num/Dfun_den


d1_nie<-(Dfun+Afun)-(Kfun+Bfun)
d2_nie<-Astar*(Dfun-Bfun)+A*(Afun-Kfun)
d3_nie<-as.vector(C)*((Dfun+Afun)-(Kfun+Bfun))
d4_nie<-0
d5_nie<-0
d6_nie<-Afun-Bfun
d7_nie<-A*(Afun-Bfun)
d8_nie<-rep(0,ncov)

gamma_logNIE<-c(d4_nie,d5_nie,d6_nie,d7_nie,d8_nie,d1_nie,d2_nie,d3_nie)

#gamma_logNIE<-c(d4_nie,d5_nie,d6_nie,d8_nie,d1_nie,d2_nie,d3_nie)


d1_te<-d1_nde+d1_nie
d2_te<-d2_nde+d2_nie
d3_te<-d3_nde+d3_nie
d4_te<-d4_nde+d4_nie
d5_te<-d5_nde+d5_nie
d6_te<-d6_nde+d6_nie
d7_te<-d7_nde+d7_nie
d8_te<-d8_nde+d8_nie

gamma_logTE<-c(d4_te,d5_te,d6_te,d7_te,d8_te,d1_te,d2_te,d3_te)

#gamma_logTE<-c(d4_te,d5_te,d6_te,d8_te,d1_te,d2_te,d3_te)



var_logCDE<-t(gamma_logCDE)%*%sigmat%*%gamma_logCDE
var_logNDE<-t(gamma_logNDE)%*%sigmat%*%gamma_logNDE
var_logNIE<-t(gamma_logNIE)%*%sigmat%*%gamma_logNIE
var_logTE<-t(gamma_logTE)%*%sigmat%*%gamma_logTE

se_logCDE<-sqrt(var_logCDE)
se_logNDE<-sqrt(var_logNDE)
se_logNIE<-sqrt(var_logNIE)
se_logTE<-sqrt(var_logTE)

############################


results<-numeric(8)

results[1]<-OR_TE
results[2]<-se_logTE
results[3]<-OR_NDE
results[4]<-se_logNDE
results[5]<-OR_NIE
results[6]<-se_logNIE
results[7]<-OR_CDE
results[8]<-se_logCDE

return(results)

}

