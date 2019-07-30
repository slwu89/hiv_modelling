#MSM model
rm(list=ls())
library(FME)
#library(ggplot2)
pars=list(S_Rec0=9000,S_DD0=10000,S_INS0=8000,pI_Rec0=0.01,pI_DD0=0.01,pI_INS0=0.01,theta_v=10,beta=0.01,delta_Rec=0.1,delta_DD=0.1, delta_INS=0.1,gamma1=2,
          gamma2=7.35,gamma3=5,gamma4=5, delta1=5, omega=0.1, mu=0.1,c_REC=300,c_DD=200,c_INS=100,p_REC=0.1,p_DD=0.3,p_INS=0.8,e=0.85,f=0.6,
          muREC=15,muDD=15,muINS=15,b_con=0.15, con_grad=0.07,con_max=0.78,treat_cov=0.1,art_factor=0.1,b_r=0.2)
FSW.ode<- function(pars,times=seq(1980,2017,by=1))
{derivs <- function(t,state,pars)
  
{ # returns rate of change
  with(as.list(c(state, pars)), 
{
  
  
  if(t<1998){
    con_use =b_con
  } else if(t>=1998&t<2006){
    con_use=b_con+con_grad*(t-1998)
  }else {
    con_use=con_max
  }
  cduse=(1-con_use*e)
  
  if(t<2006){
    omega=0
    }else {
    omega=treat_cov
  }
  
  #Retirement rates
 
  mu_DD=1/muDD
  mu_Rec=1/muREC
  mu_INS=1/muINS
 
   #progression rate
  gamma_1=12/gamma1
  gamma_2=1/gamma2
  gamma_3=12/gamma3
  gamma_4=12/gamma4
  delta=1/delta1
  
  #Population sizes
  N_REC=S_Rec+I_Rec1+I_Rec2+I_Rec3+I_Rec4+I_Rect
  N_DD=S_DD+I_DD1+I_DD2+I_DD3+I_DD4+I_DDt
  #Balancing contraints
  N_INS=(N_REC*c_REC*(1-2*p_REC)+N_DD*c_DD*(1-2*p_DD))/(c_INS*(2*p_INS-1))  
  #N_INS=S_INS+I_INS1+I_INS2+I_INS3+I_INS1+I_INSt
  
  #Total infection
  TI_REC=I_Rec1+I_Rec2+I_Rec3+I_Rec4+I_Rect
  TI_DD=I_DD1+I_DD2+I_DD3+I_DD4+I_DDt
  TI_INS=I_INS1+I_INS2+I_INS3+I_INS4+I_INSt
  TI=TI_REC+TI_DD+TI_INS
  PI_REC=TI_REC/TI*100
  PI_DD=TI_DD/TI*100
  PI_INS=TI_INS/TI*100
  
  #Mixing and transmission
  NT_REC=p_REC*N_REC+p_DD*N_DD+p_INS*N_INS
  pri_REC=p_REC*(theta_v*I_Rec1+I_Rec2+I_Rec3+I_Rec4+art_factor*I_Rect)/NT_REC
  pri_DD=p_DD*(theta_v*I_DD1+I_DD2+I_DD3+I_DD4+art_factor*I_DDt)/NT_REC
  pri_INS=p_INS*(theta_v*I_INS1+I_INS2+I_INS3+I_INS1+art_factor*I_INSt)/NT_REC
  # Insertive acts
  NT_REC1=(1-p_REC)*N_REC+(1-p_DD)*N_DD+(1-p_INS)*N_INS
  pri_REC1=(1-p_REC)*(theta_v*I_Rec1+I_Rec2+I_Rec3+I_Rec4+art_factor*I_Rect)/NT_REC1
  pri_DD1=(1-p_DD)*(theta_v*I_DD1+I_DD2+I_DD3+I_DD4+art_factor*I_DDt)/NT_REC1
  pri_INS1=(1-p_INS)*(theta_v*I_INS1+I_INS2+I_INS3+I_INS1+art_factor*I_INSt)/NT_REC1
  
  #Insertive contacts
  INS_REC=c_REC*p_REC*N_REC
  INS_DD=c_DD*p_DD*N_DD
  INS_INS=c_INS*p_INS*N_INS
  
  #Receptive contacts
  REC_REC=c_REC*(1-p_REC)*N_REC
  REC_DD=c_DD*(1-p_DD)*N_DD
  REC_INS=c_INS*(1-p_INS)*N_INS
  
  #Forces of infection
  lambda_Rec=cduse*beta*c_REC*(1-p_REC)*(pri_DD+pri_INS)+cduse*beta*b_r*c_REC*p_REC*(pri_DD1+pri_INS1)
  lambda_DD=cduse*beta*c_DD*(1-p_DD)*(pri_REC+pri_INS) +cduse*beta*b_r*c_DD*p_DD*(pri_REC1+pri_INS1) 
  lambda_INS=cduse*beta*c_INS*(1-p_INS)*(pri_REC+pri_DD) +cduse*beta*b_r*c_INS*p_INS*(pri_REC1+pri_DD1)
  
  #Equations for Receptive MSM (REC)
  
  dS_Rec=(delta+mu_Rec)*S_Rec0-lambda_Rec*S_Rec-mu_Rec*S_Rec
  dI_Rec1=lambda_Rec*S_Rec-(gamma_1+mu_Rec)*I_Rec1
  dI_Rec2=gamma_1*I_Rec1-(gamma_2+mu_Rec)*I_Rec2
  dI_Rec3=gamma_2*I_Rec2-(gamma_3+omega+mu_Rec)*I_Rec3
  dI_Rec4=gamma_3*I_Rec3-(gamma_4+omega+mu_Rec)*I_Rec4
  dI_Rect=omega*(I_Rec3+I_Rec4)-(delta+mu_Rec)*I_Rect
  
  #Equations for Double Decker MSM (DD)
  
  dS_DD=(delta+mu_DD)*S_DD0- lambda_DD*S_DD-mu_DD*S_DD
  dI_DD1=lambda_DD*S_DD-(gamma_1+mu_DD)*I_DD1
  dI_DD2=gamma_1*I_DD1-(gamma_2+mu_DD)*I_DD2
  dI_DD3=gamma_2*I_DD2-(gamma_3+omega+mu_DD)*I_DD3
  dI_DD4=gamma_3*I_DD3-(gamma_4+omega+mu_DD)*I_DD4
  dI_DDt=omega*(I_DD3+I_DD4)-(delta+mu_DD)*I_DDt
  
  #Equations for Insertive MSM (INS)
  
  dS_INS=(delta+mu_INS)*S_INS0-lambda_INS*S_INS-mu_INS*S_INS
  dI_INS1= lambda_INS*S_INS-(gamma_1+mu_INS)*I_INS1
  dI_INS2=gamma_1*I_INS1-(gamma_2+mu_INS)*I_INS2
  dI_INS3=gamma_2*I_INS2-(gamma_3+omega+mu_INS)*I_INS3
  dI_INS4=gamma_3*I_INS3-(gamma_4+omega+mu_INS)*I_INS4
  dI_INSt=omega*(I_INS3+I_INS4)-(delta+mu_INS)*I_INSt
  
  
  
  return(list(c(dS_Rec,dI_Rec1,dI_Rec2,dI_Rec3,dI_Rec4,dI_Rect,dS_DD,dI_DD1,dI_DD2,dI_DD3,dI_DD4,dI_DDt,dS_INS,dI_INS1,dI_INS2,dI_INS3,dI_INS4,dI_INSt),
              pr_REC=(I_Rec1+I_Rec2+I_Rec3+I_Rec4+I_Rect)*100/N_REC,
              pr_DD=(I_DD1+I_DD2+I_DD3+I_DD4+I_DDt)*100/N_DD,
              pr_INS=(I_INS1+I_INS2+I_INS3+I_INS4+I_INSt)*100/N_INS,
              TI_REC=I_Rec1+I_Rec2+I_Rec3+I_Rec4+I_Rect,
              TI_DD=I_DD1+I_DD2+I_DD3+I_DD4+I_DDt,
              TI_INS=I_INS1+I_INS2+I_INS3+I_INS4+I_INSt,
              TI=TI_REC+TI_DD+TI_INS,
              PI_REC=TI_REC/TI*100,
              PI_DD=TI_DD/TI*100,
              PI_INS=TI_INS/TI*100,
              INCI_REC=lambda_Rec*S_Rec,
              INCI_DD=lambda_DD*S_DD,
              INCI_INS=lambda_INS*S_INS,
              N_REC=S_Rec+I_Rec1+I_Rec2+I_Rec3+I_Rec4+I_Rect,
              N_DD=S_DD+I_DD1+I_DD2+I_DD3+I_DD4+I_DDt,
              N_INS=(N_REC*c_REC*(1-2*p_REC)+N_DD*c_DD*(1-2*p_DD))/(c_INS*(2*p_INS-1))))})}
#initial conditions
state<-c(S_Rec=pars$S_Rec0,I_Rec1=pars$pI_Rec0*pars$S_Rec0,I_Rec2=0,I_Rec3=0,I_Rec4=0,I_Rect=0,S_DD=pars$S_DD0,
         I_DD1=pars$pI_DD0*pars$S_DD0,I_DD2=0,I_DD3=0,I_DD4=0,I_DDt=0,
         S_INS=((pars$S_Rec0*pars$c_REC*(1-2*pars$p_REC)+pars$S_DD0*pars$c_DD*(1-2*pars$p_DD))/(pars$c_INS*(2*pars$p_INS-1))-pars$pI_INS0*((pars$S_Rec0*pars$c_REC*(1-2*pars$p_REC)+pars$S_DD0*pars$c_DD*(1-2*pars$p_DD))/(pars$c_INS*(2*pars$p_INS-1)))),
         I_INS1=(pars$pI_INS0*((pars$S_Rec0*pars$c_REC*(1-2*pars$p_REC)+pars$S_DD0*pars$c_DD*(1-2*pars$p_DD))/(pars$c_INS*(2*pars$p_INS-1)))),
         I_INS2=0,I_INS3=0,I_INS4=0,I_INSt=0) 
# ode solves the model by integration
return(as.data.frame(ode(y = state,times = times, func = derivs, parms = pars)))
}

out=FSW.ode(pars)
#Time series plots
R1<-cbind(out$S_Rec,out$I_Rec1,out$I_Rec2,out$I_Rec3,out$I_Rec4,out$I_Rect)
R2<-cbind(out$S_DD,out$I_DD1,out$I_DD2,out$I_DD3,out$I_DD4,out$I_DDt)
R3<-cbind(out$S_INS,out$I_INS1,out$I_INS2,out$I_INS3,out$I_INS4,out$I_INSt)

t<-out[,1]

par(mfrow=c(1,1))
matplot(t,R1,type="l",lwd = c(2,2,2,2,2,2), col=c("blue","red","yellow","black","green","brown"), xlab = "Time (years)", ylab = "Population size")
legend("topright", c("S_Rec","I_Rec1","I_Rec2","I_Rec3","I_Rec4","I_Rect"),lty = 1:6,bty="n",col =c("blue","red","yellow","black","green","brown"), lwd = c(2,2,2,2,2,2),cex = 0.6)

par(mfrow=c(1,1))
matplot(t,R2,type="l",lwd = c(2,2,2,2,2,2), col=c("blue","red","yellow","black","green","brown"), xlab = "Time (years)", ylab = "Population size")
legend("topright", c("S_DD","I_DD1","I_DD2","I_DD3","I_DD4","I_DDt"),lty = 1:6,bty="n",col =c("blue","red","yellow","black","green","brown"), lwd = c(2,2,2,2,2,2),cex = 0.6)

par(mfrow=c(1,1))
matplot(t,R3,type="l",lwd = c(2,2,2,2,2,2), col=c("blue","red","yellow","black","green","brown"), xlab = "Time (years)", ylab = "Population size")
legend("topright", c("S_INS","I_INS1","I_INS2","I_INS3","I_INS4","I_INSt"),lty = 1:6,bty="n",col =c("blue","red","yellow","black","green","brown"), lwd = c(2,2,2,2,2,2),cex = 0.6)

#parameters
#        beta  "delta_REC","delta_DD", "delta_INS","gamma1","gamma2", "gamma3","gamma4", "delta", "omega", "mu","c_REC",  "c_DD","c_INS"     S_Rec0,     S_DD0,   pI_Rec0, pI_DD0, pI_INS0,  muREC,  muDD,   muINS  p_REC, p_DD, p_INS   e      b_con,     con_grad,      con_max,    treat_cov, art_factor, b_r)
min<-c( 0.002,   0.0006,     0.001,      0.001,      1.23,     3,        3 ,     4 ,        0.8,       5,     5,     156 ,   104,     2,       6230,    2490,    0,     0,      0,      9.3,    11.4,   11.8,  0.03,  0.29,  0.7,  0.61,  0 ,   0.0667,     0.766,       0.10,       0,          0.003)
max<-c( 0.014,   0.004,      0.006,      0.30,       2.9,     12.1,      13,     14 ,       8.4,       15,   10,     442 ,   312,    208,      12460,   15565,  0.04,   0.04,   0.04,   18.6,   22.8,    23.6,  0.17,  0.46,  0.88, 0.775, 0.25,  0.0857,   0.987,     0.50,      0.3,         0.5)
parRanges<-cbind(min,max)
rownames(parRanges)= c("beta","delta_Rec","delta_DD","delta_INS","gamma1","gamma2","gamma3","gamma4","delta1","omega","mu","c_REC","c_DD","c_INS","S_Rec0","S_DD0","pI_Rec0","pI_DD0","pI_INS0","muREC","muDD","muINS","p_REC","p_DD","p_INS","e","b_con","con_grad","con_max","treat_cov","art_factor","b_r")
parRanges

# Sensitivity analysis of all parameters (LHS)
SA0=sensRange(func =FSW.ode, parms = pars, dist = "latin",sensvar = c("INCI_REC","INCI_DD","INCI_INS","PI_REC","PI_DD","PI_INS","N_REC","N_DD","N_INS","pr_REC","pr_DD","pr_INS"),parRange = parRanges, num=50000)
write.table(SA0,"SA04.txt",sep="\t")
SA1 <- summary(SA0)
par(mfrow=c(1,1))
plot(SA1, xlab = "time (years)",ylab = "HIV prevalence", mfrow = NULL,quant = F, col = c("lightblue","darkblue"), legpos = "topright")
mtext(outer = TRUE, line = -1.5, side = 1, "",cex = 1.25)

#selecting output parameters
d1 = read.table("SA04.txt", sep="\t",header = T)
#d2=subset(d1, select=c(beta,delta_Rec,	delta_DD,delta_INS,gamma1,gamma2,gamma3,gamma4,delta1,	omega,mu,	c_REC,	c_DD,	c_INS,S_Rec0,S_DD0,pI_Rec0,pI_DD0,pI_INS0,muREC,muDD,muINS,p_REC,	p_DD,p_INS,e,b_con,con_grad,con_max,treat_cov,art_factor,b_r,N_REC1980,N_DD1980, N_INS1980,INCI_REC2009,INCI_DD2009,INCI_INS2009,PI_REC2009,PI_DD2009,PI_INS2009,pr_REC,pr_DD, pr_INS))
d3=subset(d1,(pr_REC2006>=16.2&pr_REC2006<=29.2)&(pr_REC2009>=16.2&pr_REC2009<=28.8)&(pr_DD2006>=5.6&pr_DD2006<=20.0)&(pr_DD2009>=6.8&pr_DD2009<=17.4)&(pr_INS2006>=4.7&pr_INS2006<=20.6)&(pr_INS2009>=4.4&pr_INS2009<=21.8))

write.table(d3,"FitsA1.txt",sep="\t")

SS<-(d3$pr_REC2006-22.7)^2+(d3$pr_REC2009-22.5)^2+(d3$pr_DD2006-12.8)^2+(d3$pr_DD2009-12.1)^2+(d3$pr_INS2006-12.7)^2+(d3$pr_INS2009-13.1)^2
attach(d3) 
d3$sums<- SS 
detach(d3) 
d4=d3[order(d3$sums), ]
d44=subset(d4,select=c(pr_REC1980,pr_REC1981,pr_REC1982,pr_REC1983,pr_REC1984,pr_REC1985,pr_REC1986,pr_REC1987,pr_REC1988,pr_REC1989,pr_REC1990, pr_REC1991,pr_REC1992,pr_REC1993,pr_REC1994,pr_REC1995,pr_REC1996,pr_REC1997,pr_REC1998,pr_REC1999,pr_REC2000,pr_REC2001,pr_REC2002,pr_REC2003,pr_REC2004,pr_REC2005,pr_REC2006,pr_REC2007,pr_REC2008,pr_REC2009,pr_REC2010,pr_REC2011))
d4p1=subset(d4,select=c(pr_DD1980,pr_DD1981,pr_DD1982,pr_DD1983,pr_DD1984,pr_DD1985,pr_DD1986,pr_DD1987,pr_DD1988,pr_DD1989,pr_DD1990,pr_DD1991,pr_DD1992,pr_DD1993,pr_DD1994,pr_DD1995,pr_DD1996,pr_DD1997,pr_DD1998,pr_DD1999,pr_DD2000,pr_DD2001,pr_DD2002,pr_DD2003,pr_DD2004,pr_DD2005,pr_DD2006,pr_DD2007,pr_DD2008,pr_DD2009,pr_DD2010,pr_DD2011))
d4p2=subset(d4,select=c(pr_INS1980,pr_INS1981,pr_INS1982,pr_INS1983,pr_INS1984,pr_INS1985,pr_INS1986,pr_INS1987,pr_INS1988,pr_INS1989,pr_INS1990,pr_INS1991,pr_INS1992,pr_INS1993,pr_INS1994,pr_INS1995,pr_INS1996,pr_INS1997,pr_INS1998,pr_INS1999,pr_INS2000,pr_INS2001,pr_INS2002,pr_INS2003,pr_INS2004,pr_INS2005,pr_INS2006,pr_INS2007,pr_INS2008,pr_INS2009,pr_INS2010,pr_INS2011))

write.table(d4,"FitsA11.txt",sep="\t")

d4 = read.table("FitsA11.txt", sep="\t",header = T)

years <- c(1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990, 1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011)
d5 <- rbind(years, d4[1:nrow(d4), ])
d6=d44
d7=d4p1
d8=d4p2

par(mfrow = c(3,3)) 
#########################################################
plot(years,d6[1,], type = "l", ylim = c(0,60),  xlab = "Time (years)", ylab = "Receptive HIV prevalence",main="Approach 1",col = 1) 
for (i in 1:137){
  lines(years,d6[i,], col = "blue")
} 

segments(2006,16.2,2006,29.2)
arrows(2006,16.2,2006,29.2,length=0.03,angle=90,code=3)
points(2006,22.7,pch=0,cex=0.6)

segments(2009,16.2,2009,28.8)
arrows(2009,16.2,2009,28.8,length=0.03,angle=90,code=3)
points(2009,22.5,pch=0,cex=0.6)

plot(years,d7[1,], type = "l", ylim = c(0, 50),  xlab = "Time (years)", ylab = "Double Decker HIV prevalence",main="Approach 1",col = 1) 
for (i in 1:137){
  lines(years,d7[i,], col = "blue") 
} 

segments(2006,5.6,2006,20)
arrows(2006,5.6,2006,20,length=0.03,angle=90,code=3)
points(2006,12.8,pch=0,cex=0.6)

segments(2009,6.8,2009,17.4)
arrows(2009,6.8,2009,17.4,length=0.03,angle=90,code=3)
points(2009,12.1,pch=0,cex=0.6)


plot(years,d8[1,], type = "l", ylim = c(0, 50),  xlab = "Time (years)", ylab = "Insertive HIV prevalence",main="Approach 1",col = 1) 
for (i in 1:137){
  lines(years,d8[i,], col = "blue" )
}
segments(2006,4.7,2006,20.6)
arrows(2006,4.7,2006,20.6,length=0.03,angle=90,code=3)
points(2006,12.7,pch=0,cex=0.6)

segments(2009,4.4,2009,21.8)
arrows(2009,4.4,2009,21.8,length=0.03,angle=90,code=3)
points(2009,13.1,pch=0,cex=0.6)




#Approach 2
#selecting output parameters
DA = d1
b3=subset(DA,(pr_REC2006>=16.2&pr_REC2006<=29.2)&(pr_REC2009>=16.2&pr_REC2009<=28.8)&(pr_DD2006>=5.6&pr_DD2006<=20.0)&(pr_DD2009>=6.8&pr_DD2009<=17.4))

write.table(b3,"FitsA2.txt",sep="\t")

SS<-(b3$pr_REC2006-22.7)^2+(b3$pr_REC2009-22.5)^2+(b3$pr_DD2006-12.8)^2+(b3$pr_DD2009-12.1)^2
attach(b3) 
b3$sums<- SS 
detach(b3) 
b4=b3[order(b3$sums), ]
b44=subset(b4,select=c(pr_REC1980,pr_REC1981,pr_REC1982,pr_REC1983,pr_REC1984,pr_REC1985,pr_REC1986,pr_REC1987,pr_REC1988,pr_REC1989,pr_REC1990, pr_REC1991,pr_REC1992,pr_REC1993,pr_REC1994,pr_REC1995,pr_REC1996,pr_REC1997,pr_REC1998,pr_REC1999,pr_REC2000,pr_REC2001,pr_REC2002,pr_REC2003,pr_REC2004,pr_REC2005,pr_REC2006,pr_REC2007,pr_REC2008,pr_REC2009,pr_REC2010,pr_REC2011))
b4p1=subset(b4,select=c(pr_DD1980,pr_DD1981,pr_DD1982,pr_DD1983,pr_DD1984,pr_DD1985,pr_DD1986,pr_DD1987,pr_DD1988,pr_DD1989,pr_DD1990,pr_DD1991,pr_DD1992,pr_DD1993,pr_DD1994,pr_DD1995,pr_DD1996,pr_DD1997,pr_DD1998,pr_DD1999,pr_DD2000,pr_DD2001,pr_DD2002,pr_DD2003,pr_DD2004,pr_DD2005,pr_DD2006,pr_DD2007,pr_DD2008,pr_DD2009,pr_DD2010,pr_DD2011))
b4p2=subset(b4,select=c(pr_INS1980,pr_INS1981,pr_INS1982,pr_INS1983,pr_INS1984,pr_INS1985,pr_INS1986,pr_INS1987,pr_INS1988,pr_INS1989,pr_INS1990,pr_INS1991,pr_INS1992,pr_INS1993,pr_INS1994,pr_INS1995,pr_INS1996,pr_INS1997,pr_INS1998,pr_INS1999,pr_INS2000,pr_INS2001,pr_INS2002,pr_INS2003,pr_INS2004,pr_INS2005,pr_INS2006,pr_INS2007,pr_INS2008,pr_INS2009,pr_INS2010,pr_INS2011))

write.table(b4,"FitsA22.txt",sep="\t")

b4 = read.table("FitsA22.txt", sep="\t",header = T)

years <- c(1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990, 1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011)
b5 <- rbind(years, b4[1:nrow(b4), ])
b6=b44
b7=b4p1
b8=b4p2

#par(mfrow = c(1,3)) 
#########################################################
plot(years,b6[1,], type = "l", ylim = c(0,100),  xlab = "Time (years)", ylab = "Receptive HIV prevalence",main="Approach 2",col = 1) 
for (i in 1:137){
  lines(years,b6[i,], col = "orange") 
} 

segments(2006,16.2,2006,29.2)
arrows(2006,16.2,2006,29.2,length=0.03,angle=90,code=3)
points(2006,22.7,pch=0,cex=0.6)

segments(2009,16.2,2009,28.8)
arrows(2009,16.2,2009,28.8,length=0.03,angle=90,code=3)
points(2009,22.5,pch=0,cex=0.6)

plot(years,b7[1,], type = "l", ylim = c(0, 100),  xlab = "Time (years)", ylab = "Double Decker HIV prevalence",main="Approach 2",col = 1) 
for (i in 1:137){
  lines(years,b7[i,], col = "orange") 
} 

segments(2006,5.6,2006,20)
arrows(2006,5.6,2006,20,length=0.03,angle=90,code=3)
points(2006,12.8,pch=0,cex=0.6)

segments(2009,6.8,2009,17.4)
arrows(2009,6.8,2009,17.4,length=0.03,angle=90,code=3)
points(2009,12.1,pch=0,cex=0.6)


plot(years,b8[1,], type = "l", ylim = c(0, 100),  xlab = "Time (years)", ylab = "Insertive HIV prevalence",main="Approach 2",col = 1) 
for (i in 1:137){
  lines(years,b8[i,], col = "orange") 
}


#Approach 3

#selecting output parameters
DA1 = d1
c3=subset(DA1,(c_REC>=c_INS)&(pr_REC2006>=16.2&pr_REC2006<=29.2)&(pr_REC2009>=16.2&pr_REC2009<=28.8)&(pr_DD2006>=5.6&pr_DD2006<=20.0)&(pr_DD2009>=6.8&pr_DD2009<=17.4))

write.table(c3,"FitsA3.txt",sep="\t")

SS<-(c3$pr_REC2006-22.7)^2+(c3$pr_REC2009-22.5)^2+(c3$pr_DD2006-12.8)^2+(c3$pr_DD2009-12.1)^2
attach(c3) 
c3$sums<- SS 
detach(c3) 
c4=c3[order(c3$sums), ]
c44=subset(c4,select=c(pr_REC1980,pr_REC1981,pr_REC1982,pr_REC1983,pr_REC1984,pr_REC1985,pr_REC1986,pr_REC1987,pr_REC1988,pr_REC1989,pr_REC1990, pr_REC1991,pr_REC1992,pr_REC1993,pr_REC1994,pr_REC1995,pr_REC1996,pr_REC1997,pr_REC1998,pr_REC1999,pr_REC2000,pr_REC2001,pr_REC2002,pr_REC2003,pr_REC2004,pr_REC2005,pr_REC2006,pr_REC2007,pr_REC2008,pr_REC2009,pr_REC2010,pr_REC2011))
c4p1=subset(c4,select=c(pr_DD1980,pr_DD1981,pr_DD1982,pr_DD1983,pr_DD1984,pr_DD1985,pr_DD1986,pr_DD1987,pr_DD1988,pr_DD1989,pr_DD1990,pr_DD1991,pr_DD1992,pr_DD1993,pr_DD1994,pr_DD1995,pr_DD1996,pr_DD1997,pr_DD1998,pr_DD1999,pr_DD2000,pr_DD2001,pr_DD2002,pr_DD2003,pr_DD2004,pr_DD2005,pr_DD2006,pr_DD2007,pr_DD2008,pr_DD2009,pr_DD2010,pr_DD2011))
c4p2=subset(c4,select=c(pr_INS1980,pr_INS1981,pr_INS1982,pr_INS1983,pr_INS1984,pr_INS1985,pr_INS1986,pr_INS1987,pr_INS1988,pr_INS1989,pr_INS1990,pr_INS1991,pr_INS1992,pr_INS1993,pr_INS1994,pr_INS1995,pr_INS1996,pr_INS1997,pr_INS1998,pr_INS1999,pr_INS2000,pr_INS2001,pr_INS2002,pr_INS2003,pr_INS2004,pr_INS2005,pr_INS2006,pr_INS2007,pr_INS2008,pr_INS2009,pr_INS2010,pr_INS2011))

write.table(c4,"FitsA33.txt",sep="\t")

c4 = read.table("FitsA33.txt", sep="\t",header = T)

years <- c(1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990, 1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011)
c5 <- rbind(years, c4[1:nrow(c4), ])
c6=c44
c7=c4p1
c8=c4p2

#par(mfrow = c(1,3)) 
#########################################################
plot(years,c6[1,], type = "l", ylim = c(0,100),  xlab = "Time (years)", ylab = "Receptive HIV prevalence",main="Approach 3",col = 1) 
for (i in 1:137){
  lines(years,c6[i,], col = "red") 
} 

segments(2006,16.2,2006,29.2)
arrows(2006,16.2,2006,29.2,length=0.03,angle=90,code=3)
points(2006,22.7,pch=0,cex=0.6)

segments(2009,16.2,2009,28.8)
arrows(2009,16.2,2009,28.8,length=0.03,angle=90,code=3)
points(2009,22.5,pch=0,cex=0.6)

plot(years,c7[1,], type = "l", ylim = c(0, 100),  xlab = "Time (years)", ylab = "Double Decker HIV prevalence",main="Approach 3",col = 1) 
for (i in 1:137){
  lines(years,c7[i,], col = "red") 
} 

segments(2006,5.6,2006,20)
arrows(2006,5.6,2006,20,length=0.03,angle=90,code=3)
points(2006,12.8,pch=0,cex=0.6)

segments(2009,6.8,2009,17.4)
arrows(2009,6.8,2009,17.4,length=0.03,angle=90,code=3)
points(2009,12.1,pch=0,cex=0.6)


plot(years,c8[1,], type = "l", ylim = c(0, 100),  xlab = "Time (years)", ylab = "Insertive HIV prevalence",main="Approach 3",col = 1) 
for (i in 1:137){
  lines(years,c8[i,], col = "red") 
}


par(mfrow=c(1,2),mar=c(5.0,6.5,4.0,2.0))
boxplot(d4$pr_REC2009,b4$pr_REC2009,c4$pr_REC2009,
        d4$pr_DD2009,b4$pr_DD2009,c4$pr_DD2009,
        d4$pr_INS2009,b4$pr_INS2009,c4$pr_INS2009,
        outline=FALSE, whisklty = 1,xaxt="n",yaxt="n",cex.axis=1,cex.lab=2.0,
        ylim = c(0,30),cex.lab=0.85, cex.axis=0.7, cex.sub=0.7,
        at =c(1,1.25,1.45,2.5,2.75,2.95,4,4.25,4.45),
        col=c("blue","red","orange",
              "blue","red","orange",
              "blue","red","orange"),
        pars = list(boxwex = 0.18, staplewex = 0.5, outwex = 0.3))
axis(1,at=c(1.35,2.85,4.55),c("Receptive","Versatile","Insertive"),tick=T,cex.axis=1,tck=-0.02)
axis(2,las=1, tck=-0.01,tick=T,cex.axis=1)
mtext(side = 2, "HIV prevalence(%)", line = 3.5,cex=1)
legend("topleft",fill=c("blue","red","orange"),bty="n",ncol=1,
legend=c("Approach 1","Approach 2","Approach 3"),cex=0.85)
       
#par(mar=c(5.0,6.5,4.0,2.0))
boxplot(d4$N_REC2009,b4$N_REC2009,c4$N_REC2009,
        d4$N_DD2009,b4$N_DD2009,c4$N_DD2009,
        d4$N_INS2009,b4$N_INS2009,c4$N_INS2009,
        outline=FALSE, whisklty = 1,xaxt="n",yaxt="n",cex.axis=1,cex.lab=2.0,
        ylim = c(0,500000),cex.lab=0.85, cex.axis=0.7, cex.sub=0.7,
        at =c(1,1.25,1.45,2.5,2.75,2.95,4,4.25,4.45),
        col=c("blue","red","orange",
              "blue","red","orange",
              "blue","red","orange"),
        pars = list(boxwex = 0.18, staplewex = 0.5, outwex = 0.3))
axis(1,at=c(1.35,2.85,4.55),c("Receptive","Versatile","Insertive"),tick=T,cex.axis=1,tck=-0.02)
axis(2,las=1, tck=-0.01,tick=T,cex.axis=1)
mtext(side = 2, "Population size", line = 4,cex=1)
legend("topleft",fill=c("blue","red","orange"),bty="n",ncol=1,
       legend=c("Approach 1","Approach 2","Approach 3"),cex=0.85)



Approach1=colMeans(cbind(d4$PI_INS2009,d4$PI_DD2009,d4$PI_REC2009))
Approach2=colMeans(cbind(b4$PI_INS2009,b4$PI_DD2009,b4$PI_REC2009))
Approach3=colMeans(cbind(c4$PI_INS2009,c4$PI_DD2009,c4$PI_REC2009))
counts <- cbind(Approach1, Approach2,Approach3)
par(mfrow=c(1,2),mar=c(5.0,6.5,5.0,4))
barplot(counts,xaxt="n",ylab = "Relative number of incident infections",
        ylim = c(0, 100), width=0.5, col=c("blue","red","orange"), beside=FALSE)
axis(1,at=c(0.3,0.9,1.5),c("Approach 1","Approach 2","Approach 3"),tick=F,cex.axis=1,tck=-0.02)
#par(xpd=TRUE)
legend("topright",inset=c(-0.14,0), xpd = TRUE,fill=c("blue","red","orange"),bty="n",ncol=1, 
       legend=c("Insertive","Versatile","Receptive"),cex=0.69)


Insertive=colMeans(cbind((d4$INCI_INS2009/d4$INCI_INS2009),(b4$INCI_INS2009/b4$INCI_INS2009),(c4$INCI_INS2009/c4$INCI_INS2009)))
Versatile=colMeans(cbind((d4$INCI_DD2009/d4$INCI_INS2009),(b4$INCI_DD2009/b4$INCI_INS2009),(c4$INCI_DD2009/c4$INCI_INS2009)))
Receptive=colMeans(cbind((d4$INCI_REC2009/d4$INCI_INS2009),(b4$INCI_REC2009/b4$INCI_INS2009),(c4$INCI_REC2009/c4$INCI_INS2009)))

counts <- cbind(Insertive, Versatile,Receptive)
barplot(counts, ylab = "Relative number of incident infections",
        ylim = c(0, 4), col=c("blue","red","orange"), beside=TRUE)
legend("topleft",fill=c("blue","red","orange"),bty="n",ncol=1, 
       legend=c("Approach 1","Approach 2","Approach 3"),cex=0.69)
