getwd()
library(FME)
library(ggplot2)

thetaNames <- c("P_transmission", "HighV_factor", "T_factor", "C_FGMG", "C_FSWC", "C_F2M2", "Prop_FSW", "Prop_MC", "Prop_IFG", "Prop_IMG", "Prop_IFSW", "Prop_IMC", "Prop_F2", "Prop_M2", "Prop_IF2", "Prop_IM2",  "Tau",     "mu",     "mu_I.2",  "mu_T.2",  "mu_A.2", "mu_I.1",  "mu_T.1",  "mu_A.1", "omega", "rho", "con_F2.2", "con_FSW.2", "con_FG.2", "con_MG.2", "con_M2.2", "con_MC.2","con_F2.1", "con_FSW.1", "con_FG.1", "con_MG.1", "con_M2.1", "con_MC.1","con_eff" ,"circum.2" , "circum.1", "circum_eff", "sup.2.rate", "perc.2nd.line", "unsup.1.rate", "unsup.2.rate", "fall.off.2", "birth", "sigma.1", "sigma.2", "sup.1.rate.pre2013", "sup.1.rate.pre2015", "sup.1.rate.post2015", "sup.1.rate.pre2015int", "sup.1.rate.post2015int", "fall.off.1.pre2013", "fall.off.1.pre2015", "fall.off.1.post2015", "fall.off.1.pre2015int", "fall.off.1.post2015int", "phiM", "phiF")

pars = list(NF_0 = 12864738, NM_0 = 12594866, P_transmission = 0.01, HighV_factor = 10,
            T_factor = 0.08, C_FGMG = 1, C_FSWC = 100, C_F2M2 = 3,
            #C_MGFG = 1, C_CFSW = 10, C_M2F2 = 3,
            Prop_FSW = 0.08, Prop_MC = 0.088, Prop_IFG = 0.05, Prop_IMG = 0.039, Prop_IFSW = 0.266, Prop_IMC = 0.06,
            Prop_F2 = 0.125, Prop_M2 = 0.125, Prop_IF2 = 0.06, Prop_IM2 = 0.06,
            Tau = 0.242, mu = 0.001,
            mu_I.2 = 0.06, mu_T.2 = 0.01, mu_A.2 = 0.0935, omega = 0.1, rho = 0.08,
            mu_I.1 = 0.1, mu_T.1 = 0.01, mu_A.1 = 0.3,
            con_F2.2 = 0.274, con_FSW.2 = 0.632, con_FG.2 = 0.576, con_MG.2 = 0.587, con_M2.2 = 0.213, con_MC.2 = 0.45,con_eff = 0.875,
            con_F2.1 = 0.1, con_FSW.1 = 0.05, con_FG.1 = 0.05, con_MG.1 = 0.05, con_M2.1 = 0.05, con_MC.1 = 0.05,
            circum.2 = 0.54, circum.1 = 0.8, circum_eff = 0.6, sup.2.rate = 0.08, perc.2nd.line = 0.1, unsup.1.rate = 0.05,
            unsup.2.rate = 0.05, fall.off.2 = 0.06, birth = 0.015, sigma.1 = 0.05, sigma.2 = 0.4,
            sup.1.rate.pre2013 = 0.554, sup.1.rate.pre2015 = 0.634, sup.1.rate.post2015 = 0.554,
            sup.1.rate.pre2015int = 0.85, sup.1.rate.post2015int = 0.749,
            fall.off.1.pre2013 = 0.173, fall.off.1.pre2015 = 0.109, fall.off.1.post2015 = 0.173,
            fall.off.1.pre2015int = 0.009, fall.off.1.post2015int = 0.067
            )

pars.2 = list(NF_0 = 1000000, NM_0 = 1000000, P_transmission = 0.01, HighV_factor = 10,
            T_factor = 0.08, C_FGMG = 1, C_FSWC = 100, C_F2M2 = 3,
            #C_MGFG = 1, C_CFSW = 10, C_M2F2 = 3,
            Prop_FSW = 0.08, Prop_MC = 0.088, Prop_IFG = 0.05, Prop_IMG = 0.039, Prop_IFSW = 0.266, Prop_IMC = 0.06,
            Prop_F2 = 0.125, Prop_M2 = 0.125, Prop_IF2 = 0.06, Prop_IM2 = 0.06,
            Tau = 0.242, mu = 0.001,
            mu_I = 0.06, mu_T = 0.01, mu_A = 0.0935, omega = 0.1, rho = 0.08,
            con_F2 = 0.274, con_FSW = 0.632, con_FG = 0.576, con_MG = 0.587, con_M2 = 0.213, con_MC = 0.45,con_eff = 0.875,
            circum = 0.54, circum_eff = 0.6, sup.2.rate = 0.08, perc.2nd.line = 0.1, unsup.1.rate = 0.05,
            unsup.2.rate = 0.05, fall.off.2 = 0.06, birth = 0.015
)

### no intervention #####
FSW.ode<- function(pars,times=seq(1990, 2033,by=1)){

  derivs <- function(t,state,pars)

{ # returns rate of change
  with(as.list(c(state, pars)),
       {
         #circumcision
         #psi=1-circum
         #if(t<2002){
        #   phi=0
         #}else if(t>=2002&t<=2006){
          # phi=0.04
         #}else {
        #   phi=0.01
        # }

         #  #Treatment rates

         if(t<2015){
           sigma = sigma.1
         }
         else if(t>=2015){
           sigma = sigma.2
         }

         if(t<2015){
           con_FG = con_FG.1
           con_MG = con_MG.1
           con_MC = con_MC.1
           con_FSW = con_FSW.1
           con_F2 = con_F2.1
           con_M2 = con_M2.1
           circum = circum.1
           mu_A = mu_A.1
           mu_T = mu_T.1
           mu_I = mu_I.1
         }
         else if(t>=2015){
           con_FG = con_FG.2
           con_MG = con_MG.2
           con_MC = con_MC.2
           con_FSW = con_FSW.2
           con_F2 = con_F2.2
           con_M2 = con_M2.2
           circum = circum.2
           mu_A = mu_A.2
           mu_T = mu_T.2
           mu_I = mu_I.2
         }

         if(t<2013){
           sup.1.rate = sup.1.rate.pre2013
           fall.off.1 = fall.off.1.pre2013
         }
         else if(t>=2013 & t<2015){
           sup.1.rate = sup.1.rate.pre2015
           fall.off.1 = fall.off.1.pre2015
         }
         else if(t>=2015){
           sup.1.rate = sup.1.rate.post2015
           fall.off.1 = fall.off.1.post2015
         }

         # else {
         #      tau=0.14
         #    }

         #Duration
         #t_gamma=1/t_du #treatment rate
         #alpha_BB=1/df_BB
         #alpha_TS=1/df_TS
         #delta_BB=1/dm_BB
         #delta_TS=1/dm_TS

         #Population sizes
         #initial sizes
         N0_FG=NF_0 #+ 1000
         N0_MG=NM_0 #+ 1500

         #variable population sizes
         N_FG=S_FG+IV_FG+I_FG+T1_FG+T2_FG+T1S_FG+T2S_FG+A_FG
         N_MG=S_MG+IV_MG+I_MG+T1_MG+T2_FG+T1S_MG+T2S_MG+A_MG

         N_FSW = S_FSW+IV_FSW+I_FSW+T1_FSW+T2_FSW+T1S_FSW+T2S_FSW+A_FSW
         N_MC = S_MC+IV_MC+I_MC+T1_MC+T2_MC+T1S_MC+T2S_MC+A_MC

         N_F2=S_F2+IV_F2+I_F2+T1_F2+T2_F2+T1S_F2+T2S_F2+A_F2
         N_M2=S_M2+IV_M2+I_M2+T1_M2+T2_M2+T1S_M2+T2S_M2+A_M2

         N_F=N_FG + N_FSW + N_F2
         N_M=N_MG + N_MC + N_M2

         N_D = D
         N_DHIV = D_HIV
         N_Gen = D - D_HIV
         #contact constraints
         #CFBB_CBB=CCBB_FBB*N_FBB/N_CBB
         #CFTS_CTS=(1-Z)*(CCTS_FTS*N_FTS)/N_CTS
         #CFTS_CBB=Z*CCTS_FTS*N_FTS/N_CBB
         #CFGMG=P_FGPMGP*(CMGP_FGP*N_FGP)/(P_MGPFGP*N_MGP)
         #CMGP_FBB=P_MGPFBB*(CFGP_MGP*N_MGP)/(pmf*N_FBB)
         #CFGP_CBB=P_FGPCBB*(CMGP_FGP*N_FGP)/(pmc*N_CBB)
         C_MGFG = C_FGMG*N_FG/N_MG
         C_CFSW = C_FSWC*N_FSW/N_MC
         C_M2F2 = C_F2M2*N_F2/N_M2



           #General population
           beta_FGI = P_transmission*C_FGMG
           beta_FGV = P_transmission*HighV_factor*C_FGMG
           beta_FGT = P_transmission*C_FGMG
           beta_FGTS = P_transmission*T_factor*C_FGMG
           beta_FGA = P_transmission*C_FGMG
           beta_MGI = P_transmission*C_MGFG
           beta_MGV = P_transmission*HighV_factor*C_MGFG
           beta_MGT = P_transmission*C_MGFG
           beta_MGTS = P_transmission*T_factor*C_MGFG
           beta_MGA = P_transmission*C_MGFG

           #FSW and C
           beta_FSWI = P_transmission*C_FSWC
           beta_FSWV = P_transmission*HighV_factor*C_FSWC
           beta_FSWT = P_transmission*C_FSWC
           beta_FSWTS = P_transmission*T_factor*C_FSWC
           beta_FSWA = P_transmission*C_FSWC
           beta_MCI = P_transmission*C_CFSW
           beta_MCV = P_transmission*HighV_factor*C_CFSW
           beta_MCT = P_transmission*C_CFSW
           beta_MCTS = P_transmission*T_factor*C_CFSW
           beta_MCA = P_transmission*C_CFSW

           # F2 and M2
           beta_F2I = P_transmission*C_F2M2
           beta_F2V = P_transmission*HighV_factor*C_F2M2
           beta_F2T = P_transmission*C_F2M2
           beta_F2TS = P_transmission*T_factor*C_F2M2
           beta_F2A = P_transmission*C_F2M2
           beta_M2I = P_transmission*C_M2F2
           beta_M2V = P_transmission*HighV_factor*C_M2F2
           beta_M2T = P_transmission*C_M2F2
           beta_M2TS = P_transmission*T_factor*C_M2F2
           beta_M2A = P_transmission*C_M2F2
           #Force of Infection
           #General Population
           lambda_MG=beta_MGI*(I_FG/N_FG)*((1-con_FG)+con_FG*(1-con_eff))*((1-circum) + circum*(1-circum_eff))+  # should this all be a single term?
             beta_MGV*(IV_FG/N_FG)*((1-con_FG)+con_FG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_MGT*(T1_FG/N_FG)*((1-con_FG)+con_FG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_MGTS*(T1S_FG/N_FG)*((1-con_FG)+con_FG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_MGT*(T2_FG/N_FG)*((1-con_FG)+con_FG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_MGTS*(T2S_FG/N_FG)*((1-con_FG)+con_FG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_MGA*(A_FG/N_FG)*((1-con_FG)+con_FG*(1-con_eff))*((1-circum) + circum*(1-circum_eff))
           lambda_FG=beta_FGI*(I_MG/N_MG)*((1-con_MG)+con_MG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_FGV*(IV_MG/N_MG)*((1-con_MG)+con_MG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_FGT*(T1_MG/N_MG)*((1-con_MG)+con_MG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_FGTS*(T1S_MG/N_MG)*((1-con_MG)+con_MG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_FGT*(T2_MG/N_MG)*((1-con_MG)+con_MG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_FGTS*(T2S_MG/N_MG)*((1-con_MG)+con_MG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_FGA*(A_MG/N_MG)*((1-con_MG)+con_MG*(1-con_eff))*((1-circum) + circum*(1-circum_eff))
           #FSW and C
           lambda_MC=beta_MCI*(I_FSW/N_FSW)*((1-con_FSW)+con_FSW*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_MCV*(IV_FSW/N_FSW)*((1-con_FSW)+con_FSW*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_MCT*(T1_FSW/N_FSW)*((1-con_FSW)+con_FSW*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_MCTS*(T1S_FSW/N_FSW)*((1-con_FSW)+con_FSW*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_MCT*(T2_FSW/N_FSW)*((1-con_FSW)+con_FSW*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_MCTS*(T2S_FSW/N_FSW)*((1-con_FSW)+con_FSW*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_MCA*(A_FSW/N_FSW)*((1-con_FSW)+con_FSW*(1-con_eff))*((1-circum) + circum*(1-circum_eff))
           ##### need to add condom usage values for clients ######
           lambda_FSW=beta_FSWI*(I_MC/N_MC)*((1-con_MC)+con_MC*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_FSWV*(IV_MC/N_MC)*((1-con_MC)+con_MC*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_FSWT*(T1_MC/N_MC)*((1-con_MC)+con_MC*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_FSWTS*(T1S_MC/N_MC)*((1-con_MC)+con_MC*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_FSWT*(T2_MC/N_MC)*((1-con_MC)+con_MC*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_FSWTS*(T2S_MC/N_MC)*((1-con_MC)+con_MC*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_FSWA*(A_MC/N_MC)*((1-con_MC)+con_MC*(1-con_eff))*((1-circum) + circum*(1-circum_eff))
           #General Population
           lambda_M2=beta_M2I*(I_F2/N_F2)*((1-con_F2)+con_F2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_M2V*(IV_F2/N_F2)*((1-con_F2)+con_F2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_M2T*(T1_F2/N_F2)*((1-con_F2)+con_F2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_M2TS*(T1S_F2/N_F2)*((1-con_F2)+con_F2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_M2T*(T2_F2/N_F2)*((1-con_F2)+con_F2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_M2TS*(T2S_F2/N_F2)*((1-con_F2)+con_F2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_M2A*(A_F2/N_F2)*((1-con_F2)+con_F2*(1-con_eff))*((1-circum) + circum*(1-circum_eff))
           lambda_F2=beta_F2I*(I_M2/N_M2)*((1-con_M2)+con_M2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_F2V*(IV_M2/N_M2)*((1-con_M2)+con_M2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_F2T*(T1_M2/N_M2)*((1-con_M2)+con_M2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_F2TS*(T1S_M2/N_M2)*((1-con_M2)+con_M2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_F2T*(T2_M2/N_M2)*((1-con_M2)+con_M2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_F2TS*(T2S_M2/N_M2)*((1-con_M2)+con_M2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
             beta_F2A*(A_M2/N_M2)*((1-con_M2)+con_M2*(1-con_eff))*((1-circum) + circum*(1-circum_eff))
         #clients
         #lambda_CBB=beta_FM*(CFBB_CBB*15*(1-ep*f_FBB)*(theta*I_FBV+I_FBB+eta*T_FBB)/N_FBB + CFTS_CBB*20*(1-ep*f_FTS)*(theta*I_FTV+I_FTS+eta*T_FTS)/N_FTS+CFGP_CBB*54*(1-ep*f_FGP)*(theta*I_FGV+I_FGP+eta*T_FGP)/N_FGP)
         #lambda_CTS=beta_FM*CFTS_CTS*20*(1-ep*f_FTS)*(theta*I_FTV+I_FTS+eta*T_FTS)/N_FTS

         #sex workers
         #lambda_FBB=beta_MF*(CCBB_FBB*15*(1-ep*f_FBB)*(theta*I_CBV+I_CBB+eta*T_CBB)/N_CBB+CMGP_FBB*54*(1-ep*f_FGP)*(theta*I_MGV+I_MGP+eta*T_MGP)/N_MGP)
         #lambda_FTS=beta_MF*CCTS_FTS*20*(1-ep*f_FTS)*(Z*(theta*I_CBV+I_CBB+eta*T_CBB)/N_CBB+(1-Z)*(theta*I_CTV+I_CTS+eta*T_CTS)/N_CTS)

         #Tinci=(lambda_FGP*S_FGP+lambda_MGP*S_MGP+psi*lambda_MGP*S_MGC+lambda_CTS*S_CTS+psi*lambda_CTS*S_CTC+lambda_CBB*S_CBB+psi*lambda_CBB*S_CBC+lambda_FTS*S_FTS+lambda_FBB*S_FBB)*100
         #General population
         dD = mu*(S_FG + IV_FG + I_FG + T1_FG + T1S_FG+ T2_FG + T2S_FG+ A_FG +
                    S_MG + IV_MG + I_MG + T1_MG + T1S_MG + T2_MG + T2S_MG + A_MG +
                    S_FSW + IV_FSW + I_FSW + A_FSW + T1_FSW + T1S_FSW + T2_FSW + T2S_FSW +
                    S_MC + IV_MC + I_MC + T1_MC + T1S_MC + T2_MC + T2S_MC + A_MC +
                    S_F2 + IV_F2 + I_F2 + T1_F2 + T1S_F2 + T2_F2 + T2S_F2 + A_F2 +
                    S_M2 + IV_M2 + I_M2 + T1_M2 + T1S_M2 + T2_M2 + T2S_M2 + A_M2) +

           mu_I*(I_FG + I_MG  + I_FSW + I_MC + I_F2 + I_M2 +
                   T1_FG + T1_MG + T1_FSW +T1_MC + T1_F2 + T1_M2 +
                   T2_FG + T2_MG + T2_FSW +T2_MC + T2_F2 + T2_M2) +
           mu_T*(T1S_FG + T1S_MG + T1S_FSW + T1S_MC + T1S_F2 + T1S_M2 +
                   T2S_FG + T2S_MG + T2S_FSW + T2S_MC + T2S_F2 + T2S_M2) +
           mu_A*(A_FG + A_MG + A_FSW + A_MC + A_F2 + A_M2)

         dD_HIV =  mu_I*(I_FG + I_MG  + I_FSW + I_MC + I_F2 + I_M2 +
                           T1_FG + T1_MG + T1_FSW +T1_MC + T1_F2 + T1_M2 +
                           T2_FG + T2_MG + T2_FSW +T2_MC + T2_F2 + T2_M2) +
           mu_T*(T1S_FG + T1S_MG + T1S_FSW + T1S_MC + T1S_F2 + T1S_M2 +
                   T2S_FG + T2S_MG + T2S_FSW + T2S_MC + T2S_F2 + T2S_M2) +
           mu_A*(A_FG + A_MG + A_FSW + A_MC + A_F2 + A_M2)
         #Females
         dS_FG = -lambda_FG*S_FG-(mu)*S_FG + birth*N_FG#+ mu*NF_0*(1-Prop_F2 - Prop_FSW)
         dIV_FG = lambda_FG*S_FG-mu*IV_FG -Tau*IV_FG
         dI_FG= Tau * IV_FG-(mu_I+mu)*I_FG -sigma*I_FG - omega*I_FG +rho*A_FG + fall.off.1*T1_FG + fall.off.2*T2_FG
         dT1_FG= sigma*I_FG - (mu+mu_I)*T1_FG - perc.2nd.line *T1_FG - sup.1.rate*(T1_FG) + unsup.1.rate*T1S_FG
         dT1S_FG = sup.1.rate*T1_FG - unsup.1.rate*T1S_FG -(mu+mu_T)*T1S_FG
         dT2_FG = perc.2nd.line*T1_FG-sup.2.rate*T2_FG - (mu+mu_I)*T2_FG - fall.off.2*T2_FG
         dT2S_FG = sup.2.rate*T2_FG - unsup.2.rate*T2S_FG -(mu +mu_T)*T2S_FG
         dA_FG= omega*I_FG -(mu+ mu_A)*A_FG - rho*A_FG
         #Males
         dS_MG = -lambda_MG*S_MG-(mu)*S_MG +birth*N_MG#+ mu*NM_0*(1-Prop_M2 - Prop_MC)
         dIV_MG = lambda_MG*S_MG-mu*IV_MG -Tau*IV_MG
         dI_MG = Tau * IV_MG-(mu_I+mu)*I_MG -sigma*I_MG - omega*I_MG + rho*A_MG + fall.off.1*T1_MG + fall.off.2*T2_MG
         dT1_MG= sigma*I_MG - (mu+mu_I)*T1_MG - perc.2nd.line *T1_MG - sup.1.rate*(T1_MG) + unsup.1.rate*T1S_MG
         dT1S_MG = sup.1.rate*T1_MG - unsup.1.rate*T1S_MG -(mu+mu_T)*T1S_MG
         dT2_MG = perc.2nd.line*T1_MG-sup.2.rate*T2_MG - (mu+mu_I)*T2_MG - fall.off.2*T2_MG
         dT2S_MG = sup.2.rate*T2_MG - unsup.2.rate*T2S_MG -(mu +mu_T)*T2S_MG
         dA_MG = omega*I_MG -(mu+mu_A)*A_MG - rho*A_MG

         #High risk males

         #Male 2+
         dS_M2 = -lambda_M2*S_M2-(mu)*S_M2 + birth*N_M2#+ mu*NM_0*Prop_M2
         dIV_M2 = lambda_M2*S_M2-mu*IV_M2 - Tau*IV_M2
         dI_M2 = Tau * IV_M2-(mu_I+mu)*I_M2 -sigma*I_M2 - omega*I_M2 + rho*A_M2 + fall.off.1*T1_M2 + fall.off.2*T2_M2
         dT1_M2= sigma*I_M2 - (mu+mu_I)*T1_M2 - perc.2nd.line *T1_M2 - sup.1.rate*(T1_M2) + unsup.1.rate*T1S_M2
         dT1S_M2 = sup.1.rate*T1_M2 - unsup.1.rate*T1S_M2 -(mu+mu_T)*T1S_M2
         dT2_M2 = perc.2nd.line*T1_M2-sup.2.rate*T2_M2 - (mu+mu_I)*T2_M2 - fall.off.2*T2_M2
         dT2S_M2 = sup.2.rate*T2_M2 - unsup.2.rate*T2S_M2 -(mu +mu_T)*T2S_M2
         dA_M2 = omega*I_M2 - rho*A_M2 -(mu+mu_A)*A_M2
         #FSW Clients
         dS_MC = -lambda_MC*S_MC-(mu)*S_MC + birth*N_MC#+ mu*NM_0*Prop_MC
         dIV_MC = lambda_MC*S_MC-mu*IV_MC - Tau*IV_MC
         dI_MC = Tau*IV_MC-(mu_I+mu)*I_MC - sigma*I_MC - omega*I_MC +rho*A_MC + fall.off.1*T1_MC + fall.off.2*T2_MC
         dT1_MC= sigma*I_MC - (mu+mu_I)*T1_MC - perc.2nd.line *T1_MC - sup.1.rate*(T1_MC) + unsup.1.rate*T1S_MC
         dT1S_MC = sup.1.rate*T1_MC - unsup.1.rate*T1S_MC -(mu+mu_T)*T1S_MC
         dT2_MC = perc.2nd.line*T1_MC-sup.2.rate*T2_MC - (mu+mu_I)*T2_MC - fall.off.2*T2_MC
         dT2S_MC = sup.2.rate*T2_MC - unsup.2.rate*T2S_MC -(mu +mu_T)*T2S_MC
         dA_MC = omega*I_MC - rho*A_MC -(mu+mu_A)*A_MC

         #High risk females

         #Female 2+
         dS_F2 = -lambda_F2*S_F2-(mu)*S_F2 + birth*N_F2#+ mu*NF_0*Prop_F2
         dIV_F2 = lambda_F2*S_F2-mu*IV_F2 - Tau*(IV_F2)
         dI_F2= Tau * IV_F2-(mu_I+mu)*I_F2 - sigma*I_F2 - omega*I_F2 + rho*I_F2 + fall.off.1*T1_F2 + fall.off.2*T2_F2
         dT1_F2= sigma*I_F2 - (mu+mu_I)*T1_F2 - perc.2nd.line *T1_F2 - sup.1.rate*(T1_F2) + unsup.1.rate*T1S_F2
         dT1S_F2 = sup.1.rate*T1_F2 - unsup.1.rate*T1S_F2 -(mu+mu_T)*T1S_F2
         dT2_F2 = perc.2nd.line*T1_F2-sup.2.rate*T2_F2 - (mu+mu_I)*T2_F2 - fall.off.2*T2_F2
         dT2S_F2 = sup.2.rate*T2_F2 - unsup.2.rate*T2S_F2 -(mu +mu_T)*T2S_F2
         dA_F2= omega*I_F2 - rho*A_F2 -(mu+mu_A)*A_F2
         #FSWs
         dS_FSW = -lambda_FSW*S_FSW-(mu)*S_FSW + birth*N_FSW#+ mu*Prop_FSW*NF_0
         dIV_FSW = lambda_FSW*S_FSW-mu*IV_FSW - Tau*IV_FSW
         dI_FSW= Tau * IV_FSW-(mu_I+mu)*I_FSW +rho*A_FSW - sigma*I_FSW - omega*I_FSW + fall.off.1*T1_FSW + fall.off.2*T2_FSW
         dT1_FSW= sigma*I_FSW - (mu+mu_I)*T1_FSW - perc.2nd.line *T1_FSW - sup.1.rate*(T1_FSW) + unsup.1.rate*T1S_FSW
         dT1S_FSW = sup.1.rate*T1_FSW - unsup.1.rate*T1S_FSW -(mu+mu_T)*T1S_FSW
         dT2_FSW = perc.2nd.line*T1_FSW-sup.2.rate*T2_FSW - (mu+mu_I)*T2_FSW - fall.off.2*T2_FSW
         dT2S_FSW = sup.2.rate*T2_FSW - unsup.2.rate*T2S_FSW -(mu +mu_T)*T2S_FSW
         dA_FSW= omega*I_FSW - rho*A_FSW -(mu+mu_A)*A_FSW


         return(list(c(dS_FG,dIV_FG,dI_FG,dT1_FG,dT1S_FG,dT2_FG,dT2S_FG,dA_FG,
                       dS_MG,dIV_MG,dI_MG,dT1_MG,dT1S_MG,dT2_MG,dT2S_MG,dA_MG,
                       dS_FSW,dIV_FSW,dI_FSW,dT1_FSW,dT1S_FSW,dT2_FSW,dT2S_FSW,dA_FSW,
                       dS_MC,dIV_MC,dI_MC,dT1_MC,dT1S_MC,dT2_MC,dT2S_MC,dA_MC,
                       dS_F2,dIV_F2,dI_F2,dT1_F2,dT1S_F2,dT2_F2,dT2S_F2,dA_F2,
                       dS_M2,dIV_M2,dI_M2,dT1_M2,dT1S_M2,dT2_M2,dT2S_M2,dA_M2,dD, dD_HIV),
                     N_FG=S_FG+IV_FG+I_FG+T1_FG+T2_FG+T1S_FG+T2S_FG+A_FG,
                     N_MG=S_MG+IV_MG+I_MG+T1_MG+T2_MG+T1S_MG+T2S_MG+A_MG,
                     N_FSW=S_FSW+IV_FSW+I_FSW+T1_FSW+T2_FSW+T1S_FSW+T2S_FSW+A_FSW,
                     N_MC=S_MC+IV_MC+I_MC+T1_MC+T2_MC+T1S_MC+T2S_MC+A_MC,
                     N_F2=S_F2+IV_F2+I_F2+T1_F2+T1S_F2+T2_F2+T2S_F2+A_F2,
                     N_M2=S_M2+IV_M2+I_M2+T1_M2+T2_M2+T1S_M2+T2S_M2+A_M2,
                     N_F = N_FG + N_FSW + N_F2, N_M = N_MG + N_MC + N_M2, N_D=D, N_HIV = D_HIV,
                     Prev_FG=((I_FG+IV_FG+T1_FG+T2_FG+T1S_FG+T2S_FG+A_FG)/N_FG)*100,
                     Prev_MG=((I_MG+IV_MG+T1_MG+T2_MG+T1S_MG+T2S_MG+A_MG)/N_MG)*100,
                     Prev_F2=((I_F2+IV_F2+T1_F2+T1S_F2+T2_F2+T2S_F2+A_F2)/N_F2)*100,
                     Prev_M2=((I_MG+IV_MG+T1_M2+T2_M2+T1S_M2+T2S_M2+A_MG)/N_M2)*100,
                     Prev_FSW=((I_FSW+IV_FSW+T1_FSW+T2_FSW+T1S_FSW+T2S_FSW+A_FSW)/N_FSW)*100,
                     Prev_MC=((I_MC+IV_MC+T1_MC+T2_MC+T1S_MC+T2S_MC+A_MC)/N_MC)*100,
                     Prev_F=((I_FG+IV_FG+T1_FG+T2_FG+T1S_FG+T2S_FG+A_FG+
                                I_FSW+IV_FSW+T1_FSW+T2_FSW+T1S_FSW+T2S_FSW+A_FSW+
                                I_F2+IV_F2+T1_F2+T1S_F2+T2_F2+T2S_F2+A_F2)/N_F)*100,
                     Inf_Total = I_FG+IV_FG+T1_FG+T2_FG+T1S_FG+T2S_FG+A_FG+
                       I_FSW+IV_FSW+T1_FSW+T2_FSW+T1S_FSW+T2S_FSW+A_FSW+
                       I_F2+IV_F2+T1_F2+T1S_F2+T2_F2+T2S_F2+A_F2+D_HIV+
                       I_MG+IV_MG+T1_MG+T2_MG+T1S_MG+T2S_MG+A_MG+
                       I_MC+IV_MC+T1_MC+T2_MC+T1S_MC+T2S_MC+A_MC+
                       I_M2+IV_M2+T1_M2+T2_M2+T1S_M2+T2S_M2+A_M2,
                     Inf_noT = I_FG+IV_FG +A_FG+ I_MG+IV_MG+A_MG+I_F2+IV_F2+ A_F2+
                       I_M2+IV_M2 +A_M2+ I_FSW+IV_FSW +A_FSW+ I_MC+ IV_MC +A_MC+ D_HIV,
                     Inf_T = I_FG+IV_FG+T1_FG+T2_FG+T1S_FG+T2S_FG+A_FG+
                       I_FSW+IV_FSW+T1_FSW+T2_FSW+T1S_FSW+T2S_FSW+A_FSW+
                       I_F2+IV_F2+T1_F2+T1S_F2+T2_F2+T2S_F2+A_F2+D_HIV+
                       I_MG+IV_MG+T1_MG+T2_MG+T1S_MG+T2S_MG+A_MG+
                       I_MC+IV_MC+T1_MC+T2_MC+T1S_MC+T2S_MC+A_MC+
                       I_M2+IV_M2+T1_M2+T2_M2+T1S_M2+T2S_M2+A_M2 -
                       (I_FG+IV_FG +A_FG+ I_MG+IV_MG+A_MG+I_F2+IV_F2+ A_F2+
                          I_M2+IV_M2 +A_M2+ I_FSW+IV_FSW +A_FSW+ I_MC+ IV_MC +A_MC+ D_HIV),
                     Prev = ((I_MG+IV_MG+T1_MG+T2_MG+T1S_MG+T2S_MG+A_MG+
                               I_MC+IV_MC+T1_MC+T2_MC+T1S_MC+T2S_MC+A_MC+
                               I_M2+IV_M2+T1_M2+T2_M2+T1S_M2+T2S_M2+A_M2) +
                               (I_FG+IV_FG+T1_FG+T2_FG+T1S_FG+T2S_FG+A_FG+
                                I_FSW+IV_FSW+T1_FSW+T2_FSW+T1S_FSW+T2S_FSW+A_FSW+
                                  I_F2+IV_F2+T1_F2+T1S_F2+T2_F2+T2S_F2+A_F2))*100/(N_F + N_M),
                     Prev_M=((I_MG+IV_MG+T1_MG+T2_MG+T1S_MG+T2S_MG+A_MG+
                                I_MC+IV_MC+T1_MC+T2_MC+T1S_MC+T2S_MC+A_MC+
                                I_M2+IV_M2+T1_M2+T2_M2+T1S_M2+T2S_M2+A_M2)/N_M)*100,
                     Tinci=(lambda_FG*S_FG+lambda_MG*S_MG+lambda_MC*S_MC+lambda_FSW*S_FSW+lambda_F2*S_F2+lambda_M2*S_M2)*
                       100/(S_FG+S_MG+S_FSW+S_MC+S_F2+S_M2)))})
}
#initial conditions
# state  = c(S_FG=pars$NF_0*(1-pars$Prop_IFG), IV_FG=0 ,I_FG=pars$NF_0*pars$Prop_IFG,T1_FG=0, T1S_FG=0,T2_FG=0, T2S_FG=0,A_FG = 0,
#          S_MG=pars$NM_0*(1-pars$Prop_IMG), IV_MG=0, I_MG=pars$NM_0*pars$Prop_IMG, T1_MG=0,T1S_MG=0, T2_MG=0,T2S_MG=0,A_MG = 0,
#          S_FSW=pars$NF_0*pars$Prop_FSW, IV_FSW=0, I_FSW=pars$NF_0*pars$Prop_FSW*pars$Prop_IFSW, T1_FSW=0,T1S_FSW=0,T2_FSW=0,T2S_FSW=0, A_FSW=0,
#          S_MC=pars$NM_0*pars$Prop_MC, IV_MC=0, I_MC=pars$NM_0*pars$Prop_MC*pars$Prop_IMC, T1_MC=0,T1S_MC=0,T2_MC=0,T2S_MC=0, A_MC=0,
#          S_F2=pars$NF_0*pars$Prop_F2, IV_F2=0 ,I_F2=pars$NF_0*pars$Prop_F2*pars$Prop_IF2,T1_F2=0,T1S_F2=0,T2_F2=0,T2S_F2=0, A_F2 = 0,
#          S_M2=pars$NM_0*pars$Prop_M2, IV_M2=0, I_M2=pars$NM_0*pars$Prop_IM2*pars$Prop_M2, T1_M2=0,T1S_M2=0,T2_M2=0,T2S_M2=0, A_M2 = 0,D=0, D_HIV = 0)
  pars <- as.list(pars)
  NF_0 = 12864738; NM_0 = 12594866
  state  = c(S_FG=NF_0*(1-pars$Prop_IFG), IV_FG=0 ,I_FG=NF_0*pars$Prop_IFG,T1_FG=0, T1S_FG=0,T2_FG=0, T2S_FG=0,A_FG = 0,
             S_MG=NM_0*(1-pars$Prop_IMG), IV_MG=0, I_MG=NM_0*pars$Prop_IMG, T1_MG=0,T1S_MG=0, T2_MG=0,T2S_MG=0,A_MG = 0,
             S_FSW=NF_0*pars$Prop_FSW, IV_FSW=0, I_FSW=NF_0*pars$Prop_FSW*pars$Prop_IFSW, T1_FSW=0,T1S_FSW=0,T2_FSW=0,T2S_FSW=0, A_FSW=0,
             S_MC=NM_0*pars$Prop_MC, IV_MC=0, I_MC=NM_0*pars$Prop_MC*pars$Prop_IMC, T1_MC=0,T1S_MC=0,T2_MC=0,T2S_MC=0, A_MC=0,
             S_F2=NF_0*pars$Prop_F2, IV_F2=0 ,I_F2=NF_0*pars$Prop_F2*pars$Prop_IF2,T1_F2=0,T1S_F2=0,T2_F2=0,T2S_F2=0, A_F2 = 0,
             S_M2=NM_0*pars$Prop_M2, IV_M2=0, I_M2=NM_0*pars$Prop_IM2*pars$Prop_M2, T1_M2=0,T1S_M2=0,T2_M2=0,T2S_M2=0, A_M2 = 0,D=0, D_HIV = 0)


# ode solves the model by integration
return(as.data.frame(ode(y = state,times = times, func = derivs, parms = pars)))
}

#### with intervention ####
FSW.ode.int<- function(pars,times=seq(1990, 2033,by=1))
{derivs <- function(t,state,pars)

{ # returns rate of change
  with(as.list(c(state, pars)),
       {
         #circumcision
         #psi=1-circum
         #if(t<2002){
         #   phi=0
         #}else if(t>=2002&t<=2006){
         # phi=0.04
         #}else {
         #   phi=0.01
         # }

         #  #Treatment rates

         if(t<2015){
           sigma = sigma.1
         }
         else if(t>=2015){
           sigma = sigma.2
         }

         if(t<2015){
           con_FG = con_FG.1
           con_MG = con_MG.1
           con_MC = con_MC.1
           con_FSW = con_FSW.1
           con_F2 = con_F2.1
           con_M2 = con_M2.1
           circum = circum.1
           mu_A = mu_A.1
           mu_T = mu_T.1
           mu_I = mu_I.1
         }
         else if(t>=2015){
           con_FG = con_FG.2
           con_MG = con_MG.2
           con_MC = con_MC.2
           con_FSW = con_FSW.2
           con_F2 = con_F2.2
           con_M2 = con_M2.2
           circum = circum.2
           mu_A = mu_A.2
           mu_T = mu_T.2
           mu_I = mu_I.2
         }

         if(t<2013){
           sup.1.rate = sup.1.rate.pre2013
           fall.off.1 = fall.off.1.pre2013
         }
         else if(t>=2013 & t<2015){
           sup.1.rate = sup.1.rate.pre2015int
           fall.off.1 = fall.off.1.pre2015int
         }
         else if(t>=2015){
           sup.1.rate = sup.1.rate.post2015int
           fall.off.1 = fall.off.1.post2015int
         }

         # else {
         #      tau=0.14
         #    }

         #Duration
         #t_gamma=1/t_du #treatment rate
         #alpha_BB=1/df_BB
         #alpha_TS=1/df_TS
         #delta_BB=1/dm_BB
         #delta_TS=1/dm_TS

         #Population sizes
         #initial sizes
         N0_FG=NF_0 #+ 1000
         N0_MG=NM_0 #+ 1500

         #variable population sizes
         N_FG=S_FG+IV_FG+I_FG+T1_FG+T2_FG+T1S_FG+T2S_FG+A_FG
         N_MG=S_MG+IV_MG+I_MG+T1_MG+T2_FG+T1S_MG+T2S_MG+A_MG

         N_FSW = S_FSW+IV_FSW+I_FSW+T1_FSW+T2_FSW+T1S_FSW+T2S_FSW+A_FSW
         N_MC = S_MC+IV_MC+I_MC+T1_MC+T2_MC+T1S_MC+T2S_MC+A_MC

         N_F2=S_F2+IV_F2+I_F2+T1_F2+T2_F2+T1S_F2+T2S_F2+A_F2
         N_M2=S_M2+IV_M2+I_M2+T1_M2+T2_M2+T1S_M2+T2S_M2+A_M2

         N_F=N_FG + N_FSW + N_F2
         N_M=N_MG + N_MC + N_M2

         N_D = D
         N_DHIV = D_HIV
         N_Gen = D - D_HIV
         #contact constraints
         #CFBB_CBB=CCBB_FBB*N_FBB/N_CBB
         #CFTS_CTS=(1-Z)*(CCTS_FTS*N_FTS)/N_CTS
         #CFTS_CBB=Z*CCTS_FTS*N_FTS/N_CBB
         #CFGMG=P_FGPMGP*(CMGP_FGP*N_FGP)/(P_MGPFGP*N_MGP)
         #CMGP_FBB=P_MGPFBB*(CFGP_MGP*N_MGP)/(pmf*N_FBB)
         #CFGP_CBB=P_FGPCBB*(CMGP_FGP*N_FGP)/(pmc*N_CBB)
         C_MGFG = C_FGMG*N_FG/N_MG
         C_CFSW = C_FSWC*N_FSW/N_MC
         C_M2F2 = C_F2M2*N_F2/N_M2



         #General population
         beta_FGI = P_transmission*C_FGMG
         beta_FGV = P_transmission*HighV_factor*C_FGMG
         beta_FGT = P_transmission*C_FGMG
         beta_FGTS = P_transmission*T_factor*C_FGMG
         beta_FGA = P_transmission*C_FGMG
         beta_MGI = P_transmission*C_MGFG
         beta_MGV = P_transmission*HighV_factor*C_MGFG
         beta_MGT = P_transmission*C_MGFG
         beta_MGTS = P_transmission*T_factor*C_MGFG
         beta_MGA = P_transmission*C_MGFG

         #FSW and C
         beta_FSWI = P_transmission*C_FSWC
         beta_FSWV = P_transmission*HighV_factor*C_FSWC
         beta_FSWT = P_transmission*C_FSWC
         beta_FSWTS = P_transmission*T_factor*C_FSWC
         beta_FSWA = P_transmission*C_FSWC
         beta_MCI = P_transmission*C_CFSW
         beta_MCV = P_transmission*HighV_factor*C_CFSW
         beta_MCT = P_transmission*C_CFSW
         beta_MCTS = P_transmission*T_factor*C_CFSW
         beta_MCA = P_transmission*C_CFSW

         # F2 and M2
         beta_F2I = P_transmission*C_F2M2
         beta_F2V = P_transmission*HighV_factor*C_F2M2
         beta_F2T = P_transmission*C_F2M2
         beta_F2TS = P_transmission*T_factor*C_F2M2
         beta_F2A = P_transmission*C_F2M2
         beta_M2I = P_transmission*C_M2F2
         beta_M2V = P_transmission*HighV_factor*C_M2F2
         beta_M2T = P_transmission*C_M2F2
         beta_M2TS = P_transmission*T_factor*C_M2F2
         beta_M2A = P_transmission*C_M2F2
         #Force of Infection
         #General Population
         lambda_MG=beta_MGI*(I_FG/N_FG)*((1-con_FG)+con_FG*(1-con_eff))*((1-circum) + circum*(1-circum_eff))+
           beta_MGV*(IV_FG/N_FG)*((1-con_FG)+con_FG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_MGT*(T1_FG/N_FG)*((1-con_FG)+con_FG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_MGTS*(T1S_FG/N_FG)*((1-con_FG)+con_FG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_MGT*(T2_FG/N_FG)*((1-con_FG)+con_FG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_MGTS*(T2S_FG/N_FG)*((1-con_FG)+con_FG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_MGA*(A_FG/N_FG)*((1-con_FG)+con_FG*(1-con_eff))*((1-circum) + circum*(1-circum_eff))
         lambda_FG=beta_FGI*(I_MG/N_MG)*((1-con_MG)+con_MG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_FGV*(IV_MG/N_MG)*((1-con_MG)+con_MG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_FGT*(T1_MG/N_MG)*((1-con_MG)+con_MG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_FGTS*(T1S_MG/N_MG)*((1-con_MG)+con_MG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_FGT*(T2_MG/N_MG)*((1-con_MG)+con_MG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_FGTS*(T2S_MG/N_MG)*((1-con_MG)+con_MG*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_FGA*(A_MG/N_MG)*((1-con_MG)+con_MG*(1-con_eff))*((1-circum) + circum*(1-circum_eff))
         #FSW and C
         lambda_MC=beta_MCI*(I_FSW/N_FSW)*((1-con_FSW)+con_FSW*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_MCV*(IV_FSW/N_FSW)*((1-con_FSW)+con_FSW*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_MCT*(T1_FSW/N_FSW)*((1-con_FSW)+con_FSW*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_MCTS*(T1S_FSW/N_FSW)*((1-con_FSW)+con_FSW*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_MCT*(T2_FSW/N_FSW)*((1-con_FSW)+con_FSW*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_MCTS*(T2S_FSW/N_FSW)*((1-con_FSW)+con_FSW*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_MCA*(A_FSW/N_FSW)*((1-con_FSW)+con_FSW*(1-con_eff))*((1-circum) + circum*(1-circum_eff))
         ##### need to add condom usage values for clients ######
         lambda_FSW=beta_FSWI*(I_MC/N_MC)*((1-con_MC)+con_MC*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_FSWV*(IV_MC/N_MC)*((1-con_MC)+con_MC*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_FSWT*(T1_MC/N_MC)*((1-con_MC)+con_MC*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_FSWTS*(T1S_MC/N_MC)*((1-con_MC)+con_MC*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_FSWT*(T2_MC/N_MC)*((1-con_MC)+con_MC*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_FSWTS*(T2S_MC/N_MC)*((1-con_MC)+con_MC*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_FSWA*(A_MC/N_MC)*((1-con_MC)+con_MC*(1-con_eff))*((1-circum) + circum*(1-circum_eff))
         #General Population
         lambda_M2=beta_M2I*(I_F2/N_F2)*((1-con_F2)+con_F2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_M2V*(IV_F2/N_F2)*((1-con_F2)+con_F2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_M2T*(T1_F2/N_F2)*((1-con_F2)+con_F2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_M2TS*(T1S_F2/N_F2)*((1-con_F2)+con_F2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_M2T*(T2_F2/N_F2)*((1-con_F2)+con_F2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_M2TS*(T2S_F2/N_F2)*((1-con_F2)+con_F2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_M2A*(A_F2/N_F2)*((1-con_F2)+con_F2*(1-con_eff))*((1-circum) + circum*(1-circum_eff))
         lambda_F2=beta_F2I*(I_M2/N_M2)*((1-con_M2)+con_M2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_F2V*(IV_M2/N_M2)*((1-con_M2)+con_M2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_F2T*(T1_M2/N_M2)*((1-con_M2)+con_M2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_F2TS*(T1S_M2/N_M2)*((1-con_M2)+con_M2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_F2T*(T2_M2/N_M2)*((1-con_M2)+con_M2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_F2TS*(T2S_M2/N_M2)*((1-con_M2)+con_M2*(1-con_eff))*((1-circum) + circum*(1-circum_eff)) +
           beta_F2A*(A_M2/N_M2)*((1-con_M2)+con_M2*(1-con_eff))*((1-circum) + circum*(1-circum_eff))
         #clients
         #lambda_CBB=beta_FM*(CFBB_CBB*15*(1-ep*f_FBB)*(theta*I_FBV+I_FBB+eta*T_FBB)/N_FBB + CFTS_CBB*20*(1-ep*f_FTS)*(theta*I_FTV+I_FTS+eta*T_FTS)/N_FTS+CFGP_CBB*54*(1-ep*f_FGP)*(theta*I_FGV+I_FGP+eta*T_FGP)/N_FGP)
         #lambda_CTS=beta_FM*CFTS_CTS*20*(1-ep*f_FTS)*(theta*I_FTV+I_FTS+eta*T_FTS)/N_FTS

         #sex workers
         #lambda_FBB=beta_MF*(CCBB_FBB*15*(1-ep*f_FBB)*(theta*I_CBV+I_CBB+eta*T_CBB)/N_CBB+CMGP_FBB*54*(1-ep*f_FGP)*(theta*I_MGV+I_MGP+eta*T_MGP)/N_MGP)
         #lambda_FTS=beta_MF*CCTS_FTS*20*(1-ep*f_FTS)*(Z*(theta*I_CBV+I_CBB+eta*T_CBB)/N_CBB+(1-Z)*(theta*I_CTV+I_CTS+eta*T_CTS)/N_CTS)

         #Tinci=(lambda_FGP*S_FGP+lambda_MGP*S_MGP+psi*lambda_MGP*S_MGC+lambda_CTS*S_CTS+psi*lambda_CTS*S_CTC+lambda_CBB*S_CBB+psi*lambda_CBB*S_CBC+lambda_FTS*S_FTS+lambda_FBB*S_FBB)*100
         #General population
         dD = mu*(S_FG + IV_FG + I_FG + T1_FG + T1S_FG+ T2_FG + T2S_FG+ A_FG +
                    S_MG + IV_MG + I_MG + T1_MG + T1S_MG + T2_MG + T2S_MG + A_MG +
                    S_FSW + IV_FSW + I_FSW + A_FSW + T1_FSW + T1S_FSW + T2_FSW + T2S_FSW +
                    S_MC + IV_MC + I_MC + T1_MC + T1S_MC + T2_MC + T2S_MC + A_MC +
                    S_F2 + IV_F2 + I_F2 + T1_F2 + T1S_F2 + T2_F2 + T2S_F2 + A_F2 +
                    S_M2 + IV_M2 + I_M2 + T1_M2 + T1S_M2 + T2_M2 + T2S_M2 + A_M2) +

           mu_I*(I_FG + I_MG  + I_FSW + I_MC + I_F2 + I_M2 +
                   T1_FG + T1_MG + T1_FSW +T1_MC + T1_F2 + T1_M2 +
                   T2_FG + T2_MG + T2_FSW +T2_MC + T2_F2 + T2_M2) +
           mu_T*(T1S_FG + T1S_MG + T1S_FSW + T1S_MC + T1S_F2 + T1S_M2 +
                   T2S_FG + T2S_MG + T2S_FSW + T2S_MC + T2S_F2 + T2S_M2) +
           mu_A*(A_FG + A_MG + A_FSW + A_MC + A_F2 + A_M2)

         dD_HIV =  mu_I*(I_FG + I_MG  + I_FSW + I_MC + I_F2 + I_M2 +
                           T1_FG + T1_MG + T1_FSW +T1_MC + T1_F2 + T1_M2 +
                           T2_FG + T2_MG + T2_FSW +T2_MC + T2_F2 + T2_M2) +
           mu_T*(T1S_FG + T1S_MG + T1S_FSW + T1S_MC + T1S_F2 + T1S_M2 +
                   T2S_FG + T2S_MG + T2S_FSW + T2S_MC + T2S_F2 + T2S_M2) +
           mu_A*(A_FG + A_MG + A_FSW + A_MC + A_F2 + A_M2)
         #Females
         dS_FG = -lambda_FG*S_FG-(mu)*S_FG + birth*N_FG#+ mu*NF_0*(1-Prop_F2 - Prop_FSW)
         dIV_FG = lambda_FG*S_FG-mu*IV_FG -Tau*IV_FG
         dI_FG= Tau * IV_FG-(mu_I+mu)*I_FG -sigma*I_FG - omega*I_FG +rho*A_FG + fall.off.1*T1_FG + fall.off.2*T2_FG
         dT1_FG= sigma*I_FG - (mu+mu_I)*T1_FG - perc.2nd.line *T1_FG - sup.1.rate*(T1_FG) + unsup.1.rate*T1S_FG
         dT1S_FG = sup.1.rate*T1_FG - unsup.1.rate*T1S_FG -(mu+mu_T)*T1S_FG
         dT2_FG = perc.2nd.line*T1_FG-sup.2.rate*T2_FG - (mu+mu_I)*T2_FG - fall.off.2*T2_FG
         dT2S_FG = sup.2.rate*T2_FG - unsup.2.rate*T2S_FG -(mu +mu_T)*T2S_FG
         dA_FG= omega*I_FG -(mu+ mu_A)*A_FG - rho*A_FG
         #Males
         dS_MG = -lambda_MG*S_MG-(mu)*S_MG +birth*N_MG#+ mu*NM_0*(1-Prop_M2 - Prop_MC)
         dIV_MG = lambda_MG*S_MG-mu*IV_MG -Tau*IV_MG
         dI_MG = Tau * IV_MG-(mu_I+mu)*I_MG -sigma*I_MG - omega*I_MG + rho*A_MG + fall.off.1*T1_MG + fall.off.2*T2_MG
         dT1_MG= sigma*I_MG - (mu+mu_I)*T1_MG - perc.2nd.line *T1_MG - sup.1.rate*(T1_MG) + unsup.1.rate*T1S_MG
         dT1S_MG = sup.1.rate*T1_MG - unsup.1.rate*T1S_MG -(mu+mu_T)*T1S_MG
         dT2_MG = perc.2nd.line*T1_MG-sup.2.rate*T2_MG - (mu+mu_I)*T2_MG - fall.off.2*T2_MG
         dT2S_MG = sup.2.rate*T2_MG - unsup.2.rate*T2S_MG -(mu +mu_T)*T2S_MG
         dA_MG = omega*I_MG -(mu+mu_A)*A_MG - rho*A_MG

         #High risk males

         #Male 2+
         dS_M2 = -lambda_M2*S_M2-(mu)*S_M2 + birth*N_M2#+ mu*NM_0*Prop_M2
         dIV_M2 = lambda_M2*S_M2-mu*IV_M2 - Tau*IV_M2
         dI_M2 = Tau * IV_M2-(mu_I+mu)*I_M2 -sigma*I_M2 - omega*I_M2 + rho*A_M2 + fall.off.1*T1_M2 + fall.off.2*T2_M2
         dT1_M2= sigma*I_M2 - (mu+mu_I)*T1_M2 - perc.2nd.line *T1_M2 - sup.1.rate*(T1_M2) + unsup.1.rate*T1S_M2
         dT1S_M2 = sup.1.rate*T1_M2 - unsup.1.rate*T1S_M2 -(mu+mu_T)*T1S_M2
         dT2_M2 = perc.2nd.line*T1_M2-sup.2.rate*T2_M2 - (mu+mu_I)*T2_M2 - fall.off.2*T2_M2
         dT2S_M2 = sup.2.rate*T2_M2 - unsup.2.rate*T2S_M2 -(mu +mu_T)*T2S_M2
         dA_M2 = omega*I_M2 - rho*A_M2 -(mu+mu_A)*A_M2
         #FSW Clients
         dS_MC = -lambda_MC*S_MC-(mu)*S_MC + birth*N_MC#+ mu*NM_0*Prop_MC
         dIV_MC = lambda_MC*S_MC-mu*IV_MC - Tau*IV_MC
         dI_MC = Tau*IV_MC-(mu_I+mu)*I_MC - sigma*I_MC - omega*I_MC +rho*A_MC + fall.off.1*T1_MC + fall.off.2*T2_MC
         dT1_MC= sigma*I_MC - (mu+mu_I)*T1_MC - perc.2nd.line *T1_MC - sup.1.rate*(T1_MC) + unsup.1.rate*T1S_MC
         dT1S_MC = sup.1.rate*T1_MC - unsup.1.rate*T1S_MC -(mu+mu_T)*T1S_MC
         dT2_MC = perc.2nd.line*T1_MC-sup.2.rate*T2_MC - (mu+mu_I)*T2_MC - fall.off.2*T2_MC
         dT2S_MC = sup.2.rate*T2_MC - unsup.2.rate*T2S_MC -(mu +mu_T)*T2S_MC
         dA_MC = omega*I_MC - rho*A_MC -(mu+mu_A)*A_MC

         #High risk females

         #Female 2+
         dS_F2 = -lambda_F2*S_F2-(mu)*S_F2 + birth*N_F2#+ mu*NF_0*Prop_F2
         dIV_F2 = lambda_F2*S_F2-mu*IV_F2 - Tau*(IV_F2)
         dI_F2= Tau * IV_F2-(mu_I+mu)*I_F2 - sigma*I_F2 - omega*I_F2 + rho*I_F2 + fall.off.1*T1_F2 + fall.off.2*T2_F2
         dT1_F2= sigma*I_F2 - (mu+mu_I)*T1_F2 - perc.2nd.line *T1_F2 - sup.1.rate*(T1_F2) + unsup.1.rate*T1S_F2
         dT1S_F2 = sup.1.rate*T1_F2 - unsup.1.rate*T1S_F2 -(mu+mu_T)*T1S_F2
         dT2_F2 = perc.2nd.line*T1_F2-sup.2.rate*T2_F2 - (mu+mu_I)*T2_F2 - fall.off.2*T2_F2
         dT2S_F2 = sup.2.rate*T2_F2 - unsup.2.rate*T2S_F2 -(mu +mu_T)*T2S_F2
         dA_F2= omega*I_F2 - rho*A_F2 -(mu+mu_A)*A_F2
         #FSWs
         dS_FSW = -lambda_FSW*S_FSW-(mu)*S_FSW + birth*N_FSW#+ mu*Prop_FSW*NF_0
         dIV_FSW = lambda_FSW*S_FSW-mu*IV_FSW - Tau*IV_FSW
         dI_FSW= Tau * IV_FSW-(mu_I+mu)*I_FSW +rho*A_FSW - sigma*I_FSW - omega*I_FSW + fall.off.1*T1_FSW + fall.off.2*T2_FSW
         dT1_FSW= sigma*I_FSW - (mu+mu_I)*T1_FSW - perc.2nd.line *T1_FSW - sup.1.rate*(T1_FSW) + unsup.1.rate*T1S_FSW
         dT1S_FSW = sup.1.rate*T1_FSW - unsup.1.rate*T1S_FSW -(mu+mu_T)*T1S_FSW
         dT2_FSW = perc.2nd.line*T1_FSW-sup.2.rate*T2_FSW - (mu+mu_I)*T2_FSW - fall.off.2*T2_FSW
         dT2S_FSW = sup.2.rate*T2_FSW - unsup.2.rate*T2S_FSW -(mu +mu_T)*T2S_FSW
         dA_FSW= omega*I_FSW - rho*A_FSW -(mu+mu_A)*A_FSW


         return(list(c(dS_FG,dIV_FG,dI_FG,dT1_FG,dT1S_FG,dT2_FG,dT2S_FG,dA_FG,
                       dS_MG,dIV_MG,dI_MG,dT1_MG,dT1S_MG,dT2_MG,dT2S_MG,dA_MG,
                       dS_FSW,dIV_FSW,dI_FSW,dT1_FSW,dT1S_FSW,dT2_FSW,dT2S_FSW,dA_FSW,
                       dS_MC,dIV_MC,dI_MC,dT1_MC,dT1S_MC,dT2_MC,dT2S_MC,dA_MC,
                       dS_F2,dIV_F2,dI_F2,dT1_F2,dT1S_F2,dT2_F2,dT2S_F2,dA_F2,
                       dS_M2,dIV_M2,dI_M2,dT1_M2,dT1S_M2,dT2_M2,dT2S_M2,dA_M2,dD, dD_HIV),
                     N_FG=S_FG+IV_FG+I_FG+T1_FG+T2_FG+T1S_FG+T2S_FG+A_FG,
                     N_MG=S_MG+IV_MG+I_MG+T1_MG+T2_MG+T1S_MG+T2S_MG+A_MG,
                     N_FSW=S_FSW+IV_FSW+I_FSW+T1_FSW+T2_FSW+T1S_FSW+T2S_FSW+A_FSW,
                     N_MC=S_MC+IV_MC+I_MC+T1_MC+T2_MC+T1S_MC+T2S_MC+A_MC,
                     N_F2=S_F2+IV_F2+I_F2+T1_F2+T1S_F2+T2_F2+T2S_F2+A_F2,
                     N_M2=S_M2+IV_M2+I_M2+T1_M2+T2_M2+T1S_M2+T2S_M2+A_M2,
                     N_F = N_FG + N_FSW + N_F2, N_M = N_MG + N_MC + N_M2, N_D=D, N_HIV = D_HIV,
                     Prev_FG=((I_FG+IV_FG+T1_FG+T2_FG+T1S_FG+T2S_FG+A_FG)/N_FG)*100,
                     Prev_MG=((I_MG+IV_MG+T1_MG+T2_MG+T1S_MG+T2S_MG+A_MG)/N_MG)*100,
                     Prev_F2=((I_F2+IV_F2+T1_F2+T1S_F2+T2_F2+T2S_F2+A_F2)/N_F2)*100,
                     Prev_M2=((I_MG+IV_MG+T1_M2+T2_M2+T1S_M2+T2S_M2+A_MG)/N_M2)*100,
                     Prev_FSW=((I_FSW+IV_FSW+T1_FSW+T2_FSW+T1S_FSW+T2S_FSW+A_FSW)/N_FSW)*100,
                     Prev_MC=((I_MC+IV_MC+T1_MC+T2_MC+T1S_MC+T2S_MC+A_MC)/N_MC)*100,
                     Prev_F=((I_FG+IV_FG+T1_FG+T2_FG+T1S_FG+T2S_FG+A_FG+
                                I_FSW+IV_FSW+T1_FSW+T2_FSW+T1S_FSW+T2S_FSW+A_FSW+
                                I_F2+IV_F2+T1_F2+T1S_F2+T2_F2+T2S_F2+A_F2)/N_F)*100,
                     Inf_Total = I_FG+IV_FG+T1_FG+T2_FG+T1S_FG+T2S_FG+A_FG+
                       I_FSW+IV_FSW+T1_FSW+T2_FSW+T1S_FSW+T2S_FSW+A_FSW+
                       I_F2+IV_F2+T1_F2+T1S_F2+T2_F2+T2S_F2+A_F2+D_HIV+
                       I_MG+IV_MG+T1_MG+T2_MG+T1S_MG+T2S_MG+A_MG+
                       I_MC+IV_MC+T1_MC+T2_MC+T1S_MC+T2S_MC+A_MC+
                       I_M2+IV_M2+T1_M2+T2_M2+T1S_M2+T2S_M2+A_M2,
                     Inf_noT = I_FG+IV_FG +A_FG+ I_MG+IV_MG+A_MG+I_F2+IV_F2+ A_F2+
                       I_M2+IV_M2 +A_M2+ I_FSW+IV_FSW +A_FSW+ I_MC+ IV_MC +A_MC+ D_HIV,
                     Inf_T = I_FG+IV_FG+T1_FG+T2_FG+T1S_FG+T2S_FG+A_FG+
                       I_FSW+IV_FSW+T1_FSW+T2_FSW+T1S_FSW+T2S_FSW+A_FSW+
                       I_F2+IV_F2+T1_F2+T1S_F2+T2_F2+T2S_F2+A_F2+D_HIV+
                       I_MG+IV_MG+T1_MG+T2_MG+T1S_MG+T2S_MG+A_MG+
                       I_MC+IV_MC+T1_MC+T2_MC+T1S_MC+T2S_MC+A_MC+
                       I_M2+IV_M2+T1_M2+T2_M2+T1S_M2+T2S_M2+A_M2 -
                       (I_FG+IV_FG +A_FG+ I_MG+IV_MG+A_MG+I_F2+IV_F2+ A_F2+
                          I_M2+IV_M2 +A_M2+ I_FSW+IV_FSW +A_FSW+ I_MC+ IV_MC +A_MC+ D_HIV),
                     Prev_M=((I_MG+IV_MG+T1_MG+T2_MG+T1S_MG+T2S_MG+A_MG+
                                I_MC+IV_MC+T1_MC+T2_MC+T1S_MC+T2S_MC+A_MC+
                                I_M2+IV_M2+T1_M2+T2_M2+T1S_M2+T2S_M2+A_M2)/N_M)*100,
                     Prev = ((I_MG+IV_MG+T1_MG+T2_MG+T1S_MG+T2S_MG+A_MG+
                                I_MC+IV_MC+T1_MC+T2_MC+T1S_MC+T2S_MC+A_MC+
                                I_M2+IV_M2+T1_M2+T2_M2+T1S_M2+T2S_M2+A_M2) +
                               (I_FG+IV_FG+T1_FG+T2_FG+T1S_FG+T2S_FG+A_FG+
                                  I_FSW+IV_FSW+T1_FSW+T2_FSW+T1S_FSW+T2S_FSW+A_FSW+
                                  I_F2+IV_F2+T1_F2+T1S_F2+T2_F2+T2S_F2+A_F2))*100/(N_F + N_M),
                     Tinci=(lambda_FG*S_FG+lambda_MG*S_MG+lambda_MC*S_MC+lambda_FSW*S_FSW+lambda_F2*S_F2+lambda_M2*S_M2)*
                       100/(S_FG+S_MG+S_FSW+S_MC+S_F2+S_M2)))})
}
#initial conditions
state  = c(S_FG=pars$NF_0*(1-pars$Prop_IFG), IV_FG=0 ,I_FG=pars$NF_0*pars$Prop_IFG,T1_FG=0, T1S_FG=0,T2_FG=0, T2S_FG=0,A_FG = 0,
           S_MG=pars$NM_0*(1-pars$Prop_IMG), IV_MG=0, I_MG=pars$NM_0*pars$Prop_IMG, T1_MG=0,T1S_MG=0, T2_MG=0,T2S_MG=0,A_MG = 0,
           S_FSW=pars$NF_0*pars$Prop_FSW, IV_FSW=0, I_FSW=pars$NF_0*pars$Prop_FSW*pars$Prop_IFSW, T1_FSW=0,T1S_FSW=0,T2_FSW=0,T2S_FSW=0, A_FSW=0,
           S_MC=pars$NM_0*pars$Prop_MC, IV_MC=0, I_MC=pars$NM_0*pars$Prop_MC*pars$Prop_IMC, T1_MC=0,T1S_MC=0,T2_MC=0,T2S_MC=0, A_MC=0,
           S_F2=pars$NF_0*pars$Prop_F2, IV_F2=0 ,I_F2=pars$NF_0*pars$Prop_F2*pars$Prop_IF2,T1_F2=0,T1S_F2=0,T2_F2=0,T2S_F2=0, A_F2 = 0,
           S_M2=pars$NM_0*pars$Prop_M2, IV_M2=0, I_M2=pars$NM_0*pars$Prop_IM2*pars$Prop_M2, T1_M2=0,T1S_M2=0,T2_M2=0,T2S_M2=0, A_M2 = 0,D=0, D_HIV = 0)

# ode solves the model by integration
return(as.data.frame(ode(y = state,times = times, func = derivs, parms = pars)))
}


out=FSW.ode(pars)
out = noInt
out.2 = FSW.ode.int(pars)
out.2 = Int
t = out[,1]
deaths = cbind(out$D, out.2$D)
deaths.hiv = cbind(out$D_HIV, out.2$D_HIV)
out$D_HIV - out.2$D_HIV
matplot(t,deaths.hiv,type="l", lty = c(1,2), lwd = c(2,2), col =c("blue","red"), xlab = "Time (years)", ylab = "Deaths", main = "HIV Deaths- int vs ctrl")
legend(x = "topleft", legend = c("control", "intervention"), col = c("blue", "red"), lty = c(1,2),
       bty = "n", xjust = 0)
out$D_HIV- out.2$D_HIV
dim(deaths)
# after 10 years closed system 2.27729% reductions in deaths after 1 year int
# after 5 years closed system 1.1578811% reductions in deaths after 1 year int

# after 5 years 2.0% reductions in deaths after 2 year int
# after 10 years closed system 2.678348% reductions in deaths after 2 year int


### observed prevalence ### doesnt change each run
prev.mat = matrix(c(4.1, 3.9, 3.7, 3.6, 6.1, 6, 5.9, 5.8, 5.1, 5, 4.8, 4.7), nrow = 4)
years = c(2013, 2014, 2015, 2016)
matplot(years, prev.mat[,-3], xlim = c(2013, 2016),
        ylim = c(0, 7), type = "l", lwd = 2, col = c("blue", "red", "green"))
legend("bottomright", c("Prev_M", "Prev_F"), col = c("blue", "red"),
       lty = c(1,2), lwd = 2, cex = 0.5, bty = "n")


################## changes each run #######
prev = cbind(out$Prev_F, out$Prev_M, out.2$Prev_F, out.2$Prev_M)
matplot(t,prev,type="l",lty = c(1,1,2,2),lwd = c(1,1,2,2), col =c("red","blue", "red", "blue"),
        xlab = "Time (years)", ylab = "Prev", main = "Prevalence by Gender")
legend("topright", c("Prev_F", "Prev_M", "Prev_FI", "Prev_MI"), col = c("red","blue", "red", "blue"), lty = c(1,1,2,2),
       bty = "n", cex = 0.6, lwd = c(1,1,2,2))
out.2$Prev_M-out$Prev_M
sim.prevs = prev[c(1,3,5,7),c(2,1)]
prevs = cbind(prev.mat[,-3], sim.prevs)
matplot(years, prevs, xlim = c(2013, 2016),
        ylim = c(0, 15), type = "l", lwd = 2, col = c("blue", "red"), lty = c(1,1,2,2),
        main = "Observed vs. Simulated: no intervention")
legend("bottomright", c("Male Observed", "Female Observed", "Male Simulated", "Female Simulated"), col = c("blue", "red"),
       lty = c(1,1,2,2), lwd = 2, cex = 0.4, bty = "n")

infections = cbind(out$Inf_Total+out$D_HIV, out.2$Inf_Total+out.2$D_HIV)
# after five years 0.9301303% reduction in infections with 1 year int.

matplot(t,infections,type="l",lty = c(1,2),lwd = c(1,2), col =c("red","blue"),
        xlab = "Time (years)", ylab = "Infections", main = "number of infections")
legend("bottomright", c("no intervention", "intervention"), col = c("red","blue"), lty = c(1,1,2,2),
       bty = "n", cex = 0.6, lwd = c(1,1,2,2))
(out$Inf_Total+out$D_HIV) -  (out.2$Inf_Total+out.2$D_HIV)
## after 5 years with 1 year int, 0.9301303% reduction in numbers infected
## after 10 years with 1 year int, 1.310972% reduction in numbers infected

## after 5 years with 2 year int, 1.390414% reduction in numbers infected
## after 10 years with 2 year int, 1.685699% reduction in numbers infected

### infections on treatment as opposed to not
inf.no.treat = cbind(out$Inf_noT/out$Inf_Total, out.2$Inf_noT/out.2$Inf_Total)
matplot(t,inf.no.treat,type="l",lty = c(1,2),lwd = c(1,2),
        col =c("red","blue"), xlab = "Time (years)", ylab = "Infections",
        main = "Percent Infected not on Treatment")
legend("topright", c("no intervention", "intervention"), col = c("red","blue", "red", "blue"), lty = c(1,1,2,2),
       bty = "n", cex = 0.6, lwd = c(1,1,2,2))
inf.treat = cbind(out$Inf_T/out$Inf_Total, out.2$Inf_T/out.2$Inf_Total)
matplot(t,inf.treat,type="l",lty = c(1,2),lwd = c(1,2), col =c("red","blue"),
        xlab = "Time (years)", ylab = "Infections", main = "Percent Infected")
legend("topright", c("no intervention", "intervention"), col = c("red","blue"), lty = c(1,2),
       bty = "n", cex = 0.6, lwd = c(1))
out$Inf_T/out$Inf_Total- out.2$Inf_T/out.2$Inf_Total

out$I_FSW

#### art initiates ###
art = cbind(out$T1_F2+out$T1_FG+out$T1_FSW+out$T1_M2+out$T1_MC+out$T1_MG, out.2$T1_F2+out.2$T1_FG+out.2$T1_FSW+out.2$T1_M2+out.2$T1_MC+out.2$T1_MG)
matplot(t,art,type="l",lty = c(1,2),lwd = c(2,2), col =c("red","blue"), xlab = "Time (years)", ylab = "ART #", main = "number on ART")
legend("topright", c("ART no intervention", "ART with intervention"), col = c("red","blue"), lty = c(1,2),
       bty = "n", cex = 0.6, lwd = c(2,2))
out.2$Prev_M-out$Prev_M

#### incidence ###
inc = cbind(out$Tinci, out.2$Tinci)
out$Tinci-out.2$Tinci
matplot(t,inc, type = "l", lty = c(1,2), lwd = c(2,2), col = c("red", "blue"),
        xlab = "time (years)",ylab = "Incidence", main = "Incidence")
legend(x = "topright", legend = c("no intervention", "intervention"), col = c("red", "blue"), lty = 1, bty = "n")
#Time series plots
R0<-cbind(out$N_F,out$N_M)
R1<-cbind(out$S_FG,out$IV_FG, out$I_FG, out$T_FG, out$A_FG,
          out$S_MG,out$IV_MG, out$I_MG, out$T_MG, out$A_MG)
R2 = cbind(out$S_FSW,out$IV_FSW,out$I_FSW,out$T_FSW,out$A_FSW,
           out$S_MC,out$IV_MC,out$I_MC,out$T_MC,out$A_MC)
#Prevalences
R4<-cbind(out$Prev_F,out$Prev_M,out$Prev_FG,out$Prev_MG,out$Prev_FSW,out$Prev_MC)

#High Viraemia
R7<-cbind(out$S_FG,out$IV_FG,out$S_MG,out$IV_MG)

t<-out[,1]
par(mfrow=c(1,1))
matplot(t,R0,type="l",lwd = c(2,2), col =c("blue","red"), xlab = "Time (years)", ylab = "Population size")
legend("topright", c("N_F","N_M"),lty = 1:2,bty="n",col =c("blue","red"), lwd = c(2,2),cex = 0.6)

#plotting general population
par(mfrow=c(1,1))
matplot(t,R1,type="l",lwd = c(1,1,1,1,1,1,1,1,1,1), col =c("blue","red","black","yellow", "green","blue","red","black","yellow", "green"), xlab = "Time (years)", ylab = "Population size")
legend("center", c("S_FG","IV_FG", "I_FG","T_FG","A_FG",
                   "S_MG","IV_MG","I_MG", "T_MG", "A_MG"),lty = 1:4,bty="n",col =c("blue","red","black","yellow", "green"), lwd = c(1,1,1,1,1),cex = 0.3)

#plotting FSW and MC
par(mfrow=c(1,1))
matplot(t,R2,type="l",lwd = c(1,1,1,1,1,1,1,1,1,1), col =c("red","red","red","red", "red","blue","blue","blue","blue", "blue"), xlab = "Time (years)", ylab = "Population size")
legend("center", c("S_FSW","IV_FSW", "I_FSW","T_FSW","A_FSW",
                   "S_MC","IV_MC","I_MC", "T_MC", "A_MC"),lty = 1:4,bty="n",col =c("red","red","red","red", "red","blue","blue","blue","blue", "blue"), lwd = c(1,1,1,1,1),cex = 0.4)
#Prevalences
par(mfrow=c(1,1))
matplot(t,R4,type="l",lty = c(1,2,3,1,2,3), lwd = c(2,2,2,2,2,2), col =c("blue","red"), xlab = "Time (years)", ylab = "Prevalence")
legend("right", c("Prev_F","Prev_M","Prev_FG","Prev_MG","Prev_FSW","Prev_MC"),lty = c(1,2,3,1,2,3), bty="n",col =c("blue","red"), lwd = c(2,2,2,2,2,2),cex = 0.6)

#########    plotting females   #######
#susceptibles
par(mfrow = c(1,1))
R.Sus.Fem = cbind(out$S_FG, out$S_F2, out$S_FSW)
matplot(t, R.Sus.Fem, type = "l", lty = c(2,2,2), lwd = c(2,2,2), col = c("red", "blue", "green"),
        xlab = "Time (years)", ylab = "Population Size")
legend("topright", c("S_FG", "S_F2", "S_FSW"), col = c("red", "blue", "green"), lty = c(2,2,2),
       bty = "n", cex = 0.6, lwd = c(2,2,2))
# high viremia
par(mfrow = c(1,1))
R.HV.Fem = cbind(out$IV_FG, out$IV_F2, out$IV_FSW)
matplot(t, R.HV.Fem, type = "l", lty = c(2,2,2), lwd = c(2,2,2), col = c("red", "blue", "green"),
        xlab = "Time (years)", ylab = "Population Size")
legend("topright", c("IV_FG", "IV_F2", "IV_FSW"), col = c("red", "blue", "green"), lty = c(2,2,2),
       bty = "n", cex = 0.6, lwd = c(2,2,2))
#infected
par(mfrow = c(1,1))
R.Inf.Fem = cbind(out$I_FG, out$I_F2, out$I_FSW)
matplot(t, R.Inf.Fem, type = "l", lty = c(2,2,2), lwd = c(2,2,2), col = c("red", "blue", "green"),
        xlab = "Time (years)", ylab = "Population Size")
legend("topright", c("I_FG", "I_F2", "I_FSW"), col = c("red", "blue", "green"), lty = c(2,2,2),
       bty = "n", cex = 0.6, lwd = c(2,2,2))

#AIDS
par(mfrow = c(1,1))
R.Aids.Fem = cbind(out$A_FG, out$A_F2, out$A_FSW)
matplot(t, R.Aids.Fem, type = "l", lty = c(2,2,2), lwd = c(2,2,2), col = c("red", "blue", "green"),
        xlab = "Time (years)", ylab = "Population Size")
legend("topright", c("A_FG", "A_F2", "A_FSW"), col = c("red", "blue", "green"), lty = c(2,2,2),
       bty = "n", cex = 0.6, lwd = c(2,2,2))

#treatment
par(mfrow = c(1,1))
R.T.Fem = cbind(out$T_FG, out$T_F2, out$T_FSW)
matplot(t, R.T.Fem, type = "l", lty = c(2,2,2), lwd = c(2,2,2), col = c("red", "blue", "green"),
        xlab = "Time (years)", ylab = "Population Size")
legend("topright", c("T_FG", "T_F2", "T_FSW"), col = c("red", "blue", "green"), lty = c(2,2,2),
       bty = "n", cex = 0.6, lwd = c(2,2,2))

#########    plotting females   #######
#susceptibles
par(mfrow = c(1,1))
R.Sus.Male = cbind(out$S_MG, out$S_M2, out$S_MC)
matplot(t, R.Sus.Male, type = "l", lty = c(2,2,2), lwd = c(2,2,2), col = c("red", "blue", "green"),
        xlab = "Time (years)", ylab = "Population Size")
legend("topright", c("S_MG", "S_M2", "S_MC"), col = c("red", "blue", "green"), lty = c(2,2,2),
       bty = "n", cex = 0.6, lwd = c(2,2,2))
# high viremia
par(mfrow = c(1,1))
R.HV.Male = cbind(out$IV_MG, out$IV_M2, out$IV_MC)
matplot(t, R.HV.Male, type = "l", lty = c(2,2,2), lwd = c(2,2,2), col = c("red", "blue", "green"),
        xlab = "Time (years)", ylab = "Population Size")
legend("topright", c("IV_MG", "IV_M2", "IV_MC"), col = c("red", "blue", "green"), lty = c(2,2,2),
       bty = "n", cex = 0.6, lwd = c(2,2,2))
#infected
par(mfrow = c(1,1))
R.Inf.Male = cbind(out$I_MG, out$I_M2, out$I_MC)
matplot(t, R.Inf.Male, type = "l", lty = c(2,2,2), lwd = c(2,2,2), col = c("red", "blue", "green"),
        xlab = "Time (years)", ylab = "Population Size")
legend("topright", c("I_MG", "I_M2", "I_MC"), col = c("red", "blue", "green"), lty = c(2,2,2),
       bty = "n", cex = 0.6, lwd = c(2,2,2))

#AIDS
par(mfrow = c(1,1))
R.Aids.Male = cbind(out$A_MG, out$A_M2, out$A_MC)
matplot(t, R.Aids.Male, type = "l", lty = c(2,2,2), lwd = c(2,2,2), col = c("red", "blue", "green"),
        xlab = "Time (years)", ylab = "Population Size")
legend("topright", c("A_MG", "A_M2", "A_MC"), col = c("red", "blue", "green"), lty = c(2,2,2),
       bty = "n", cex = 0.6, lwd = c(2,2,2))

#treatment
par(mfrow = c(1,1))
R.T.Male = cbind(out$T_MG, out$T_M2, out$T_MC)
matplot(t, R.T.Male, type = "l", lty = c(2,2,2), lwd = c(2,2,2), col = c("red", "blue", "green"),
        xlab = "Time (years)", ylab = "Population Size")
legend("topright", c("T_MG", "T_M2", "T_MC"), col = c("red", "blue", "green"), lty = c(2,2,2),
       bty = "n", cex = 0.6, lwd = c(2,2,2))

##### Female and Male ######
#Susceptible
par(mfrow = c(1,1))
R.Sus = cbind(R.Sus.Fem, R.Sus.Male)
matplot(t, R.Sus, type = "l", lty = c(2,2,2,1,1,1), lwd = c(2,2,2,2,2,2), col = c("red", "blue", "green", "red", "blue", "green"),
        xlab = "Time (years)", ylab = "Population Size")
legend("topright", c("S_FG", "S_F2", "S_FSW", "S_MG", "S_M2", "S_MC"), col = c("red", "blue", "green", "red", "blue", "green"), lty = c(2,2,2,1,1,1),
       bty = "n", cex = 0.4, lwd = c(2,2,2))

#High Viremia
par(mfrow = c(1,1))
R.HV = cbind(R.HV.Fem, R.HV.Male)
matplot(t, R.HV, type = "l", lty = c(2,2,2,1,1,1), lwd = c(2,2,2,2,2,2), col = c("red", "blue", "green", "red", "blue", "green"),
        xlab = "Time (years)", ylab = "Population Size")
legend("topright", c("IV_FG", "IV_F2", "IV_FSW", "IV_MG", "IV_M2", "IV_MC"), col = c("red", "blue", "green", "red", "blue", "green"), lty = c(2,2,2,1,1,1),
       bty = "n", cex = 0.4, lwd = c(2,2,2))
#Infected
par(mfrow = c(1,1))
R.Inf = cbind(R.Inf.Fem, R.Inf.Male)
matplot(t, R.Inf, type = "l", lty = c(2,2,2,1,1,1), lwd = c(2,2,2,2,2,2), col = c("red", "blue", "green", "red", "blue", "green"),
        xlab = "Time (years)", ylab = "Population Size")
legend("topright", c("I_FG", "I_F2", "I_FSW", "I_MG", "I_M2", "I_MC"), col = c("red", "blue", "green", "red", "blue", "green"), lty = c(2,2,2,1,1,1),
       bty = "n", cex = 0.4, lwd = c(2,2,2))
#AIDS
par(mfrow = c(1,1))
R.Aids = cbind(R.Aids.Fem, R.Aids.Male)
matplot(t, R.Aids, type = "l", lty = c(2,2,2,1,1,1), lwd = c(2,2,2,2,2,2), col = c("red", "blue", "green", "red", "blue", "green"),
        xlab = "Time (years)", ylab = "Population Size")
legend("topright", c("A_FG", "A_F2", "A_FSW", "A_MG", "A_M2", "A_MC"), col = c("red", "blue", "green", "red", "blue", "green"), lty = c(2,2,2,1,1,1),
       bty = "n", cex = 0.4, lwd = c(2,2,2))
#treatment
par(mfrow = c(1,1))
R.Treat = cbind(R.T.Fem, R.T.Male)
matplot(t, R.Treat, type = "l", lty = c(2,2,2,1,1,1), lwd = c(2,2,2,2,2,2), col = c("red", "blue", "green", "red", "blue", "green"),
        xlab = "Time (years)", ylab = "Population Size")
legend("topright", c("T_FG", "T_F2", "T_FSW", "T_MG", "T_M2", "T_MC"), col = c("red", "blue", "green", "red", "blue", "green"), lty = c(2,2,2,1,1,1),
       bty = "n", cex = 0.4, lwd = c(2,2,2))

########### new from uganda share ###############






#v     "P_transmission*, HighV_factor*, T_factor*, C_FGMG, C_FSWC, C_F2M2, Prop_FSW, Prop_MC, Prop_IFG, Prop_IMG, Prop_IFSW, Prop_IMC, Prop_F2, Prop_M2, Prop_IF2, Prop_IM2,  Tau,     mu,     mu_I.2,  mu_T.2,  mu_A.2, mu_I.1,  mu_T.1,  mu_A.1, omega, rho, con_F2.2, con_FSW.2, con_FG.2, con_MG.2, con_M2.2, con_MC.2,con_F2.1, con_FSW.1, con_FG.1, con_MG.1, con_M2.1, con_MC.1,con_eff*,circum.2 , circum.1, circum_eff*, sup.2.rate, perc.2nd.line, unsup.1.rate, unsup.2.rate, fall.off.2, birth, sigma.1, sigma.2, sup.1.rate.pre2013, sup.1.rate.pre2015, sup.1.rate.post2015, sup.1.rate.pre2015int, sup.1.rate.post2015.int, fall.1.off.pre2013, fall.1.off.pre2015, fall.1.off.post2015, fall.1.off.pre2015int, fall.1.off.post2015int"
min = c(0.006,           5,            0.02,         0.5,   70   ,  2     ,    0.08,   0.06,      0.035 ,   0.03  ,    0.18   ,   0.04,    0.1  ,  0.1 ,    0.03,     0.05    ,0.103,  0.0005,   0.03,   0.007,   0.09,   0.05,   0.009,    0.07,   0.09 , 0.09, 0.2  ,    0.5  ,     0.4  ,    0.4  ,    0.25 ,    0.25  ,    0.1,    0.02,      0.05,     0.05,     0.05,      0.04,    0.8  ,  0.45  ,   0.4,         0.5      ,  0.03     ,    0.06      ,   0.02      ,      0.02   ,     0.02  , 0.01,  0.05,    0.4,       0.46,                0.540,                0.46,                 0.808,                 0.7,                    0.1116,             0.0628,             0.1116,              0.0028,                 0.0438)
max = c(0.03,            15,           0.16,         1.5,   150  ,  7     ,    0.1,     0.08,      0.055 ,  0.05  ,    0.35   ,   0.1,     0.15 , 0.15,     0.08,      0.15   ,0.381,  0.0015,   0.06,   0.015,   0.18,   0.15,   0.05,      0.4,   0.15 , 0.15, 0.4 ,    0.75  ,     0.65 ,    0.65 ,    0.4  ,    0.6  ,     0.5,     0.5,       0.5,      0.4,      0.5,       0.4,    0.95 ,  0.95  ,   0.85,        0.6      ,  0.15     ,    0.15      ,   0.1       ,      0.1    ,     0.1   , 0.03,  0.5,     0.9 ,       0.664,               0.718,                0.718,                0.884,                 0.792,                  0.2587,             0.183,              0.2587,              0.0268,                 0.1)


### with ten parameter sets
#v     "P_transmission*, HighV_factor*, T_factor*, C_FGMG, C_FSWC, C_F2M2, Prop_FSW, Prop_MC, Prop_IFG, Prop_IMG, Prop_IFSW, Prop_IMC, Prop_F2, Prop_M2, Prop_IF2, Prop_IM2,  Tau,     mu,     mu_I.2,  mu_T.2,  mu_A.2, mu_I.1,  mu_T.1,  mu_A.1, omega, rho, con_F2.2, con_FSW.2, con_FG.2, con_MG.2, con_M2.2, con_MC.2,con_F2.1, con_FSW.1, con_FG.1, con_MG.1, con_M2.1, con_MC.1,con_eff*,circum.2 , circum.1, circum_eff*, sup.2.rate, perc.2nd.line, unsup.1.rate, unsup.2.rate, fall.off.2, birth, sigma.1, sigma.2, sup.1.rate.pre2013, sup.1.rate.pre2015, sup.1.rate.post2015, sup.1.rate.pre2015int, sup.1.rate.post2015.int, fall.1.off.pre2013, fall.1.off.pre2015, fall.1.off.post2015, fall.1.off.pre2015int, fall.1.off.post2015int, phiM, phiF"
minNew = c(0.006,           5.5,          0,            0.5,   70   ,  2     ,    0.06,   0.07,      0.04 ,    0.03  ,    0.18   ,   0.04,    0.1  ,  0.1 ,    0.03,     0.01     ,0.103,  0.0005,   0.03,    0.007,   0.07,   0.03,   0.007,    0.07,   0.05 , 0.05, 0.2  ,    0.5  ,     0.4  ,    0.4  ,    0.15 ,    0.2  ,     0.02,    0.02,      0.02,     0.02,     0.02,      0.02,    0.8  ,  0.4   ,   0.4,         0.5      ,  0.03     ,    0.04      ,   0.02      ,      0.02   ,     0.02  , 0.01,  0,       0.05,       0.46,                0.540,                0.46,                 0.808,                 0.7,                    0.1116,             0.0628,             0.1116,              0.0028,                 0.0438, 0, 0)
maxNew = c(0.06,            24,           0.16,         1.5,   150  ,  7     ,    0.1,     0.1,      0.06 ,    0.05  ,    0.35   ,   0.1,     0.15 , 0.15,     0.1,      0.2     ,0.381,    0.0015,   0.07,   0.015,   0.18,   0.15,   0.1,      0.4,    0.15 , 0.15, 0.4 ,    0.75  ,     0.65 ,    0.65 ,    0.4  ,    0.6  ,     0.5,     0.5,       0.5,      0.5,      0.5,       0.5,     0.95 ,  0.95  ,   0.85,        0.6      ,  0.15     ,    0.15      ,   0.1       ,      0.1    ,     0.1   , 0.03,  0.5,     0.9 ,       0.664,               0.718,                0.718,                0.884,                 0.792,                  0.2587,             0.183,              0.2587,              0.0268,                 0.1, 1, 1)

parRanges<-cbind(min,max)
parRangesNew = cbind(minNew, maxNew)
rownames(parRanges)= c("P_transmission", "HighV_factor", "T_factor", "C_FGMG", "C_FSWC", "C_F2M2", "Prop_FSW", "Prop_MC", "Prop_IFG", "Prop_IMG", "Prop_IFSW", "Prop_IMC", "Prop_F2", "Prop_M2", "Prop_IF2", "Prop_IM2",  "Tau",     "mu",     "mu_I.2",  "mu_T.2",  "mu_A.2", "mu_I.1",  "mu_T.1",  "mu_A.1", "omega", "rho", "con_F2.2", "con_FSW.2", "con_FG.2", "con_MG.2", "con_M2.2", "con_MC.2","con_F2.1", "con_FSW.1", "con_FG.1", "con_MG.1", "con_M2.1", "con_MC.1","con_eff" ,"circum.2" , "circum.1", "circum_eff", "sup.2.rate", "perc.2nd.line", "unsup.1.rate", "unsup.2.rate", "fall.off.2", "birth", "sigma.1", "sigma.2", "sup.1.rate.pre2013", "sup.1.rate.pre2015", "sup.1.rate.post2015", "sup.1.rate.pre2015int", "sup.1.rate.post2015int", "fall.off.1.pre2013", "fall.off.1.pre2015", "fall.off.1.post2015", "fall.off.1.pre2015int", "fall.off.1.post2015int", "phiM", "phiF")
rownames(parRangesNew) = c("P_transmission", "HighV_factor", "T_factor", "C_FGMG", "C_FSWC", "C_F2M2", "Prop_FSW", "Prop_MC", "Prop_IFG", "Prop_IMG", "Prop_IFSW", "Prop_IMC", "Prop_F2", "Prop_M2", "Prop_IF2", "Prop_IM2",  "Tau",     "mu",     "mu_I.2",  "mu_T.2",  "mu_A.2", "mu_I.1",  "mu_T.1",  "mu_A.1", "omega", "rho", "con_F2.2", "con_FSW.2", "con_FG.2", "con_MG.2", "con_M2.2", "con_MC.2","con_F2.1", "con_FSW.1", "con_FG.1", "con_MG.1", "con_M2.1", "con_MC.1","con_eff" ,"circum.2" , "circum.1", "circum_eff", "sup.2.rate", "perc.2nd.line", "unsup.1.rate", "unsup.2.rate", "fall.off.2", "birth", "sigma.1", "sigma.2", "sup.1.rate.pre2013", "sup.1.rate.pre2015", "sup.1.rate.post2015", "sup.1.rate.pre2015int", "sup.1.rate.post2015int", "fall.off.1.pre2013", "fall.off.1.pre2015", "fall.off.1.post2015", "fall.off.1.pre2015int", "fall.off.1.post2015int", "phiM", "phiF")
parMeans = as.data.frame(t(colMeans(d45)))

mat = matrix(c(12864738, 12594866), nrow = 2)
rownames(mat) = c("NF_0", "NM_0")
dat = as.data.frame(mat)
parMeansNew = rbind(parMeans, dat)
nrow(parMeans)
parRanges
dim(parRanges)
# Sensitivity analysis of all parameters (LHS)
SA0=sensRange(func=FSW.ode, parms = pars, dist = "latin", sensvar = c("N_F","N_M","Prev_F","Prev_M", "Prev"), parRange = parRanges, num=80000)
SA0Means=sensRange(func=FSW.ode, parms = c(parMeans, NF_0 = 12864738, NM_0 = 12594866), dist = "latin", sensvar = c("N_F","N_M","Prev_F","Prev_M", "Prev"), parRange = parRangesNew, num=80000)

write.table(SA0,"SA04new.txt",sep="\t")
SA1 <- summary(SA0)
SA1Means = summary(SA0Means)
par(mfrow=c(1,1))
plot(SA1, xlab = "time (years)",ylab = "HIV prevalence", mfrow = NULL,quant = F, col = c("lightblue","darkblue"), legpos = "topright")



mtext(outer = TRUE, line = -1.5, side = 1, "",cex = 1.25)


#selecting output parameters
paramFits = d45
save(paramFits, file = "paramFits.RData")
d1 = read.table("SA04new.txt", sep="\t",header = T)
prev.minM = c(4.4, 5.1, 5.7, 6.2, 6.6, 6.9, 7, 7.1, 7, 6.9, 6.7, 6.5, 6.1, 5.7, 5.3, 5, 4.7, 4.4, 4.2, 4, 3.9, 3.7, 3.6, 3.5, 3.3, 3.2, 3.1, 2.7)
prev.maxM = c(5.9, 6.8, 7.5, 8.1, 8.7, 9, 9.2, 9.3, 9.2, 8.9, 8.6, 8.3, 7.9, 7.5, 7.1, 6.7, 6.4, 6, 5.8, 5.5, 5.3, 5.1, 4.8, 4.6, 4.5, 4.3, 4.1, 3.6)
prev.minF = c(5, 5.8, 6.7, 7.5, 8.1, 8.6, 8.9, 9.1, 9.2, 9.1, 8.9, 8.6, 8.2, 7.7, 7.3, 6.8, 6.4, 6.1, 5.9, 5.7, 5.5, 5.4, 5.3, 5.2, 5.1, 5, 4.9, 5.8)
prev.maxF = c(6.6, 7.8, 8.9, 9.8, 10.7, 11.3, 11.7, 11.9, 12, 11.8, 11.4, 11.1, 10.7, 10.2, 9.7, 9.2, 8.8, 8.4, 8.1, 7.8, 7.5, 7.3, 7.1, 7, 6.9, 6.7, 6.6, 6.9)
prev.min = c(4.7, 5.4, 6.2, 6.8, 7.3, 7.7, 7.9, 8.1, 8.1, 8, 7.8, 7.5, 7.1, 6.7, 6.3, 5.9,5.5, 5.3, 5.1, 4.9, 4.7, 4.6, 4.4, 4.3, 4.2, 4.1, 4, 4.4)
prev.max = c(6.3, 7.3, 8.2, 8.9, 9.6, 10.1, 10.4, 10.6, 10.5, 10.3, 10, 9.7, 9.3, 8.8, 8.4, 8, 7.6, 7.2, 6.9, 6.6, 6.4, 6.2, 6, 5.8, 5.7, 5.5, 5.4, 5.1)

prevM = c(5.1, 5.9, 6.5, 7.1, 7.5, 7.9, 8.1, 8.1, 8.1, 7.9, 7.7, 7.4, 7.1, 6.7,
          6.3, 6, 5.6, 5.3, 5.1, 4.8, 4.6, 4.4, 4.2, 4.1, 3.9, 3.7, 3.6, 3.1)
prevF = c(5.7, 6.7, 7.7, 8.6, 9.3, 9.9, 10.3, 10.5, 10.5, 10.4, 10.2, 9.9, 9.5,
          9.1, 8.6, 8.2, 7.7, 7.3, 7.1, 6.8, 6.6, 6.4, 6.2, 6.1, 6, 5.9, 5.8, 6.4)
matplot(cbind(prevM,prevF), type = "l")
NM.5 = c(12594866,14818600,16914713,19441222,22749655,26627475,27474896,28339804,29232512)
NF.5 = c(12864738,15142176,17263329,19969323,23348936,27252482,28097305,28970215,29858880)
Nyears = c(1990, 1995, 2000, 2005, 2010, 2015, 2016, 2017, 2018)
NM = c(12594866, 14818600, 16914713, 19441222, 22749655 ,26627475,27474896,28339804,29232512)
NF = c(12864738, 15142176, 17263329, 19969323, 23348936, 27252482,28097305, 28970215,29858880)
d3 = subset(d1,(Prev_F1990>=prev.minF[1] & Prev_F1990 <= prev.maxF[1])&(Prev_F1991>=prev.minF[2] & Prev_F1991 <= prev.maxF[2])&(Prev_F1992>=prev.minF[3] & Prev_F1992 <= prev.maxF[3])&(Prev_F1993>=prev.minF[4] & Prev_F1993 <= prev.maxF[4])&(Prev_F1994>=prev.minF[5] & Prev_F1994 <= prev.maxF[5])&(Prev_F1995>=prev.minF[6] & Prev_F1995 <= prev.maxF[6])&(Prev_F1996>=prev.minF[7] & Prev_F1996 <= prev.maxF[7])&(Prev_F1997>=prev.minF[8] & Prev_F1997 <= prev.maxF[8])&(Prev_F1998>=prev.minF[9] & Prev_F1998 <= prev.maxF[9])&(Prev_F1999>=prev.minF[10] & Prev_F1999 <= prev.maxF[10])&(Prev_F2000>=prev.minF[11] & Prev_F2000 <= prev.maxF[11])&(Prev_F2001>=prev.minF[12] & Prev_F2001 <= prev.maxF[12])&(Prev_F2002>=prev.minF[13] & Prev_F2002 <= prev.maxF[13])&(Prev_F2003>=prev.minF[14] & Prev_F2003 <= prev.maxF[14])&(Prev_F2004>=prev.minF[15] & Prev_F2004 <= prev.maxF[15])&(Prev_F2005>=prev.minF[16] & Prev_F2005 <= prev.maxF[16])&(Prev_F2006>=prev.minF[17] & Prev_F2006 <= prev.maxF[17])&(Prev_F2007>=prev.minF[18] & Prev_F2007 <= prev.maxF[18])&(Prev_F2008>=prev.minF[19] & Prev_F2008 <= prev.maxF[19])&(Prev_F2009>=prev.minF[20] & Prev_F2009 <= prev.maxF[20])&(Prev_F2010>=prev.minF[21] & Prev_F2010 <= prev.maxF[21])&(Prev_F2011>=prev.minF[22] & Prev_F2011 <= prev.maxF[22])&(Prev_F2012>=prev.minF[23] & Prev_F2012 <= prev.maxF[23])&(Prev_F2013>=prev.minF[24] & Prev_F2013 <= prev.maxF[24])&(Prev_F2014>=prev.minF[25] & Prev_F2014 <= prev.maxF[25])&(Prev_F2015>=prev.minF[26] & Prev_F2015 <= prev.maxF[26])&(Prev_F2016>=prev.minF[27] & Prev_F2016 <= prev.maxF[27])&(Prev_F2017>=prev.minF[28] & Prev_F2017 <= prev.maxF[28])&
              (Prev_M1990>=prev.minM[1] & Prev_M1990 <= prev.maxM[1])&(Prev_M1991>=prev.minM[2] & Prev_M1991 <= prev.maxM[2])&(Prev_M1992>=prev.minM[3] & Prev_M1992 <= prev.maxM[3])&(Prev_M1993>=prev.minM[4] & Prev_M1993 <= prev.maxM[4])&(Prev_M1994>=prev.minM[5] & Prev_M1994 <= prev.maxM[5])&(Prev_M1995>=prev.minM[6] & Prev_M1995 <= prev.maxM[6])&(Prev_M1996>=prev.minM[7] & Prev_M1996 <= prev.maxM[7])&(Prev_M1997>=prev.minM[8] & Prev_M1997 <= prev.maxM[8])&(Prev_M1998>=prev.minM[9] & Prev_M1998 <= prev.maxM[9])&(Prev_M1999>=prev.minM[10] & Prev_M1999 <= prev.maxM[10])&(Prev_M2000>=prev.minM[11] & Prev_M2000 <= prev.maxM[11])&(Prev_M2001>=prev.minM[12] & Prev_M2001 <= prev.maxM[12])&(Prev_M2002>=prev.minM[13] & Prev_M2002 <= prev.maxM[13])&(Prev_M2003>=prev.minM[14] & Prev_M2003 <= prev.maxM[14])&(Prev_M2004>=prev.minM[15] & Prev_M2004 <= prev.maxM[15])&(Prev_M2005>=prev.minM[16] & Prev_M2005 <= prev.maxM[16])&(Prev_M2006>=prev.minM[17] & Prev_M2006 <= prev.maxM[17])&(Prev_M2007>=prev.minM[18] & Prev_M2007 <= prev.maxM[18])&(Prev_M2008>=prev.minM[19] & Prev_M2008 <= prev.maxM[19])&(Prev_M2009>=prev.minM[20] & Prev_M2009 <= prev.maxM[20])&(Prev_M2010>=prev.minM[21] & Prev_M2010 <= prev.maxM[21])&(Prev_M2011>=prev.minM[22] & Prev_M2011 <= prev.maxM[22])&(Prev_M2012>=prev.minM[23] & Prev_M2012 <= prev.maxM[23])&(Prev_M2013>=prev.minM[24] & Prev_M2013 <= prev.maxM[24])&(Prev_M2014>=prev.minM[25] & Prev_M2014 <= prev.maxM[25])&(Prev_M2015>=prev.minM[26] & Prev_M2015 <= prev.maxM[26])&(Prev_M2016>=prev.minM[27] & Prev_M2016 <= prev.maxM[27])&(Prev_M2017>=prev.minM[28] & Prev_M2017 <= prev.maxM[28]))
d3 = subset(d1,(Prev_F1990>=prev.minF[1] & Prev_F1990 <= prev.maxF[1])&(Prev_F1995>=prev.minF[6] & Prev_F1995 <= prev.maxF[6])&(Prev_F2000>=prev.minF[11] & Prev_F2000 <= prev.maxF[11])&(Prev_F2005>=prev.minF[16] & Prev_F2005 <= prev.maxF[16])&(Prev_F2010>=prev.minF[21] & Prev_F2010 <= prev.maxF[21])&(Prev_F2015>=prev.minF[26] & Prev_F2015 <= prev.maxF[26])&
              (Prev_M1990>=prev.minM[1] & Prev_M1990 <= prev.maxM[1])&(Prev_M1995>=prev.minM[6] & Prev_M1995 <= prev.maxM[6])&(Prev_M2000>=prev.minM[11] & Prev_M2000 <= prev.maxM[11])&(Prev_M2005>=prev.minM[16] & Prev_M2005 <= prev.maxM[16])&(Prev_M2010>=prev.minM[21] & Prev_M2010 <= prev.maxM[21])&(Prev_M2015>=prev.minM[26] & Prev_M2015 <= prev.maxM[26]))
d3 = subset(d1, (Prev1990<= prev.max[1] & Prev1990>= prev.min))
sum(Prev1990 <= prev.max[1])
#d3=subset(d1,(Prev_F2005>=6&Prev_F2005<=14)&(Prev_M2005>=6&Prev_M2005<=14))
#d3=subset(d1,(Prev_F2005>=6&Prev_F2005<=14)&(Prev_M2005>=6&Prev_M2005<=14)&(Tinci2005>=1.02&Tinci2005<=1.68)&(Tinci2006>=0.82&Tinci2006<=1.33)&(Tinci2008>=1.02&Tinci2008<=1.54)&(Tinci2009>=1.03&Tinci2009<=1.5)&(Tinci2011>=0.67&Tinci2011<=1.07))
#d3=subset(d1,(Prev_F1992>=11&Prev_F1992<=14)&(Prev_M1992>=8&Prev_M1992<=12)&(Prev_F2001>=8&Prev_F2001<=11)&(Prev_M2001>=6&Prev_M2001<=9)&(Prev_F2005>=6&Prev_F2005<=14)&(Prev_M2005>=6&Prev_M2005<=14)&(Tinci2005>=0.62&Tinci2005<=1.22)&(Tinci2006>=0.79&Tinci2006<=1.37)&(Tinci2008>=0.57&Tinci2008<=0.98)&(Tinci2009>=0.65&Tinci2009<=1.04)&(Tinci2011>=0.72&Tinci2011<=1.16))
#d3=subset(d,(Prev_F2005>=0&Prev_F2005<=14)&(Prev_M2005>=0&Prev_M2005<=14)&(Tinci2005>=1.02&Tinci2005<=1.68)&(Tinci2006>=0.82&Tinci2006<=1.68)&(Tinci2008>=0.92&Tinci2008<=1.54)&(Tinci2009>=0.9&Tinci2009<=1.68)&(Tinci2011>=0.67&Tinci2011<=1.68))
#d3=subset(d,(Prev_F2005>=0&Prev_F2005<=140)&(Prev_M2005>=0&Prev_M2005<=140)&(Tinci2005>=0&Tinci2005<=150)&(Tinci2006>=0&Tinci2006<=150)&(Tinci2008>=0 &Tinci2008<=150)&(Tinci2009>=0 &Tinci2009<=150)&(Tinci2011>=0 &Tinci2011<=150))
#d3 = subset(d1, (Prev2012 >= 0.01&Prev2005 <= 0.2))
load("d3R.RData")
d3
write.table(d3,"Fits.txt",sep="\t")
d45 = subset(d3, select = c(P_transmission, HighV_factor, T_factor, C_FGMG, C_FSWC, C_F2M2, Prop_FSW, Prop_MC, Prop_IFG, Prop_IMG, Prop_IFSW, Prop_IMC, Prop_F2, Prop_M2, Prop_IF2, Prop_IM2,  Tau,     mu,     mu_I.2,  mu_T.2,  mu_A.2, mu_I.1,  mu_T.1,  mu_A.1, omega, rho, con_F2.2, con_FSW.2, con_FG.2, con_MG.2, con_M2.2, con_MC.2,con_F2.1, con_FSW.1, con_FG.1, con_MG.1, con_M2.1, con_MC.1,con_eff,circum.2 , circum.1, circum_eff, sup.2.rate, perc.2nd.line, unsup.1.rate, unsup.2.rate, fall.off.2, birth, sigma.1, sigma.2, sup.1.rate.pre2013, sup.1.rate.pre2015, sup.1.rate.post2015, sup.1.rate.pre2015int, sup.1.rate.post2015int, fall.off.1.pre2013, fall.off.1.pre2015, fall.off.1.post2015, fall.off.1.pre2015int, fall.off.1.post2015int))
write.table(d45,"parameterz.txt",sep="\t")
years <- c(1990, 1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013, 2014, 2015, 2016, 2017)
d5 <- rbind(years, d4[1:nrow(d4), ])
length(d3[1,])
plot(years,d45[1,1:28], type = "l", ylim = c(0,100),  xlab = "Time (years)", ylab = "HIV prevalence",col = 1)
dim(d3)
SS<-(d3$Prev_F1995-9.9)^2+(d3$Prev_F2000-10.2)^2+(d3$Prev_F2005-8.2)^2+(d3$Prev_F2010-6.6)^2+(d3$Prev_F2015-5.9)^2+
  (d3$Prev_M1995-7.9)^2+(d3$Prev_M2000-7.7)^2+(d3$Prev_M2005-6)^2+(d3$Prev_M2010-4.6)^2+(d3$Prev_M2015-3.7)^2
attach(d3)
d3$sums<- SS
detach(d3)
d4=d3[order(SS), ]
dim(d3)
length(SS)
d4p1=subset(d4,select=c(Prev_F1990, Prev_F1991,Prev_F1992,Prev_F1993,Prev_F1994,Prev_F1995,Prev_F1996,Prev_F1997,Prev_F1998,Prev_F1999,Prev_F2000,Prev_F2001,Prev_F2002,Prev_F2003,Prev_F2004,Prev_F2005,Prev_F2006,Prev_F2007,Prev_F2008,Prev_F2009,Prev_F2010,Prev_F2011,Prev_F2012,Prev_F2013,Prev_F2014,Prev_F2015,Prev_F2016,Prev_F2017))
d4p2=subset(d4,select=c(Prev_M1990, Prev_M1991,Prev_M1992,Prev_M1993,Prev_M1994,Prev_M1995,Prev_M1996,Prev_M1997,Prev_M1998,Prev_M1999,Prev_M2000,Prev_M2001,Prev_M2002,Prev_M2003,Prev_M2004,Prev_M2005,Prev_M2006,Prev_M2007,Prev_M2008,Prev_M2009,Prev_M2010,Prev_M2011,Prev_M2012,Prev_M2013,Prev_M2014,Prev_M2015,Prev_M2016,Prev_M2017))
par(mfrow = c(1,1))
plot(years,d4p2[1,], type = "l", ylim = c(0,30),  xlab = "Time (years)", ylab = "HIV Prev M",col = 1)
dim(d4p1)
for (i in 1:nrow(d4p2)){
  lines(years,d4p2[i,], col = topo.colors(nrow(d4p1))[i])
}

test = FSW.ode(c(d4, NF_0 = 12864738, NM_0 = 12594866))
plot(t, test$Prev_F2)
colnames(d45)
dim(d3)
head(d3)
test = FSW.ode(c(d45[1,], NF_0 = 100000, NM_0 = 100000))
out = FSW.ode(c(d45[1,], NF_0 = 12864738, NM_0 = 12594866))
out.2 = FSW.ode.int(c(d45[1,], NF_0 = 12864738, NM_0 = 12594866))
dim(d4)
head(d4)

mat = matrix(cbind(out$Tinci,out$D_HIV, art[,1], out$Inf_Total + out$D_HIV,(out$Inf_Total+out$D_HIV)/(out$N_F + out$N_M + out$D),out$Inf_T/out$Inf_Total,
        out.2$Tinci,out.2$D_HIV, art[,2], (out.2$Inf_Total+out.2$D_HIV)/(out.2$N_F + out.2$N_M + out.2$D),out.2$Inf_Total + out.2$D_HIV,out.2$Inf_T/out.2$Inf_Total), ncol = 12)
rownames(mat) = c(1990:2033)
colnames(mat) =c("Incidence", "HIV Deaths", "Number on ART", "Total Infection #", "Percent Infected","Percent on Treatment",
         "IncidenceINT", "HIV DeathsINT", "Number on ARTINT", "Total Infection #INT", "Percent InfectedINT","Percent on TreatmentINT")
write.csv(mat, file = "output.csv", row.names = TRUE)
?write.csv
tail(out$Inf_Total + out$D_HIV - out.2$Inf_Total - out.2$D_HIV)[6]

### BayesianTools Approach ###
library(BayesianTools)
pop = cbind(NF, NM)
priors = createUniformPrior(lower = minNew, upper = maxNew)

estBetaParams <- function(mu, phi) {
  var <- (mu*(1-mu))*phi
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}


lfunc = function(theta){
  names(theta) <- thetaNames
  fit = FSW.ode(theta)
  prev.mat = matrix(0, nrow = nrow(prev), ncol = 2)
  pop.vec = c(1,6,11,16,21,26,27,28,29)
  pop.mat = matrix(0, nrow = length(pop.vec),ncol = 2)

  for(i in 1:nrow(prev)){
    betaparamsM = estBetaParams(mu = fit$Prev_M[i]*1e-2, phi = theta["phiM"])
    betaparamsF = estBetaParams(mu = fit$Prev_F[i]*1e-2, phi = theta["phiF"])
    prev.mat[i,] =  c(dbeta(prev[i,1]*1e-2,betaparamsF$alpha, betaparamsF$beta, log = T), dbeta(prev[i,2]*1e-2,betaparamsM$alpha, betaparamsM$beta, log = T))
  }

  for(j in 1:nrow(pop)){
    pop.mat[j,] = c(dpois(pop[j,1], fit$N_F[pop.vec[j]], log = T), dpois(pop[j,2], fit$N_M[pop.vec[j]], log = T))
  }

  return(sum(prev.mat) + sum(pop.mat))
}

setup = createBayesianSetup(likelihood = lfunc, prior  = priors, names = rownames(parRangesNew))
settings = list(iterations = 1000000, message = T)
DEzs_out = runMCMC(bayesianSetup = setup, sampler = "DEzs", settings = settings)
save(DEzs_out, file = "AfyaMCMCPars.RData")
load("AfyaMCMCPars.RData")
summary(DEzs_out)
plot(DEzs_out)
lfunc(theta = c(unlist(paramFits[1,]),0.1,0.1))
pars.mcmc = MAP(DEzs_out)
NF_0 = 12864738
NM_0 = 12594866
pars.mcmc2 = c(NF_0, NM_0,pars.mcmc$parametersMAP[1:60])
length(pars.mcmc2)
out = FSW.ode(pars.mcmc2)
plot(run)
