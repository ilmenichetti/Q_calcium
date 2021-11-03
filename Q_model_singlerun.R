# inputs=time_series_input_interp$`Littertrap C`*0.01
# beta=6.72
# eta_11=0.34
# e0=0.29
# fc=0.48
# u0=1.13
# q0=1.1
# tmax=1.06 
# spinup=100 
# delay=0


#Q decomposition model, fraction remaining
Q_decomposition<-function(inputs, #input vector
                          beta=6.72,
                          eta_11=0.34,
                          e0=0.29,
                          fc=0.48,
                          u0=1.13,
                          q0=1.1,
                          tmax=1.06, 
                          spinup=100, 
                          delay){
  
  
  inputs_spun=c(rep(mean(inputs, na.rm=T),spinup),inputs)
  
  OM_matrix<-matrix(NA, nrow = length(inputs_spun), ncol = length(inputs_spun))#mat.or.vec(length(inputs_spun), length(inputs_spun))
  RESP_matrix<-matrix(NA, nrow = length(inputs_spun), ncol = length(inputs_spun))#mat.or.vec(length(inputs_spun), length(inputs_spun))
  
  colnames(OM_matrix)<-seq(1:length(inputs_spun))
  rownames(OM_matrix)<-seq(1:length(inputs_spun))
  
  for(t in 1:(length(inputs_spun)-delay)){
    
    zeta = (1-e0)/(beta*eta_11*e0)
    alpha = fc*beta*eta_11*u0*q0^beta;
    dec_vec=c()
    resp_vec=c()
    
    for(i in 1:length(inputs_spun)){
      #for all the times when time<tmax
      if(i<tmax){
        dec_vec[i] <- ((2/(tmax))*(1/(alpha*(1-zeta)))*((1+alpha*i)^(1-zeta)-
                                                          (1-(i/(tmax))))+
                         ((2/(tmax)^2)*(1/(alpha^2*(1-zeta)*(2-zeta)))*(1-(1+alpha*i)^(2-zeta)))+
                         (1-(i/(tmax)))^2)
      }else{
        #for all the times when time>=tmax
        dec_vec[i] <-  (2/(tmax))*(1/(alpha*(1-zeta)))*(1+alpha*i)^(1-zeta)+
          ((2/((tmax)^2))*(1/(alpha^2*(1-zeta)*(2-zeta)))*(((1+alpha*(i-(tmax)))^(2-zeta))-((1+alpha*i)^(2-zeta))))
      }#closing the if
      resp_vec[i]=1-dec_vec[i]
      
    }#closing the for (i)
  
    OM_matrix[(t+delay):length(inputs_spun),t] = (inputs_spun[t]*dec_vec)[1:(length(inputs_spun)-(t+delay)+1)]
    RESP_matrix[(t+delay):length(inputs_spun),t] = c(NA,inputs_spun[t]*diff(resp_vec))[1:(length(inputs_spun)-(t+delay)+1)]
    
  }#closing the for (t)
  
  OM_vec=rowSums(OM_matrix, na.rm=T)
  RESP_vec=rowSums(RESP_matrix, na.rm=T)
  OM_vec[(spinup+1):length(OM_vec)][1]
  
  (list(Litter=OM_vec[(spinup+1):length(OM_vec)], Respiration=RESP_vec[(spinup+1):length(RESP_vec)]))
  
} #closing the function


# 
# dec_litter_run<-Q_decomposition(inputs=I_F_C, spinup = 2000)
# dec_litter_run1<-Q_decomposition(inputs=I_F_C, spinup = 1000)
# dec_litter_run2<-Q_decomposition(inputs=I_F_C, spinup = 100)
# dec_CWD_run<-Q_decomposition(inputs=I_D_C, spinup = 2000)
# dec_litter=dec_litter_run$OM_decomposition
# dec_litter1=dec_litter_run1$OM_decomposition
# dec_litter2dec_litter_run2$OM_decomposition
# dec_CWD=Q_decomposition(inputs=I_D_C)$OM_decomposition
# resp_litter=dec_litter_run1$Respiration
# resp_CWD=dec_CWD_run$Respiration
# 
# plot(dec_litter, type="l", ylim=c(0,max(dec_litter)*1.1))
# lines(dec_litter1, lty=2)
# lines(dec_litter2, lty=3)
# 
# plot(resp_litter, type="l")
# 
# plot(dec_CWD, type="l")
# plot(resp_CWD, type="l")
# 

# biomass=time_series_input_interp$Biomass #biomass vector, annual values
# I_F_C=time_series_input_interp$`Littertrap C`*0.01#annual C inputs, litterfacll, in g/m2 converted to Mg/ha
# I_D_C=time_series_input_interp$CWD#annual C inputs, CWD, 
# ## parameters
# rho_alive=1/8.3
# rho_dead=1/13
# P_Ca=meq2Mg(time_series_input_interp$`Ca++ through`) #net calcium deposition every year
# spinup=20
# e_tmax=0.6
# delay=0


Q_Ca<-function(#variables:
  biomass, 
  I_F_C,#annual C inputs, litterfacll, 
  I_D_C,#annual C inputs, CWD, 
  r, #climatic scaling
  ## parameters
  rho_alive,
  rho_dead,
  P_Ca, #net calcium deposition every year
  spinup, 
  e_tmax,
  delay){
  
  
  ###### Annual biomass increments
  DeltaB<-c(mean(diff(biomass[1:5])), diff(biomass))
  
  ###### CARBON WITH Q
  Litter_run=Q_decomposition(inputs=I_F_C, spinup = spinup, tmax = 1.06 * e_tmax, delay=1)
  CWD_run=Q_decomposition(inputs=I_D_C, spinup = spinup, tmax = 1.06 * e_tmax, delay=delay)
  
  dec_litter=Litter_run$Litter
  dec_CWD=CWD_run$Litter
  
  resp_litter=Litter_run$Respiration
  resp_CWD=CWD_run$Respiration
  
  
  ###### CALCIUM
  
  #calcium in the litter
  F_Ca = dec_litter*rho_dead
  
  #calcium from the decomposing organic matter (litter and CWD)
  xi_Ca = resp_litter*rho_dead + resp_CWD*rho_dead
  
  #Calcium into the vegetation
  #Veg_Ca = pmax(DeltaB*0.45, 0)*rho_alive

  #total Calcium outflux from the catchment
  S_ca = xi_Ca + P_Ca #- Veg_Ca
  

return(list(outflux_Ca=S_ca, litter_Ca=F_Ca, decomp_Ca=xi_Ca, litter_C=(dec_litter+dec_CWD))) #the functions output is the simulated Ca in the catchment outflux

}

# plot(Litter_run$Litter, type="l", ylim=c(0,1.5))
# lines(Litter_run$Respiration, col="blue")
# 
# plot(CWD_run$Litter, type="l", ylim=c(0,18.5))
# lines(CWD_run$Respiration, col="blue")
# 
# CWD_run$Litter[18]+Litter_run$Litter[18]
# 
# rownames(time_series_input_interp)[18]

