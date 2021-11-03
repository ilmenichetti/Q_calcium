
library(readxl)
library(forecast)
library(zoo)
library(dfoptim)
library(hydroGOF)
library(matrixStats)

## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

time_series<-as.data.frame(read_excel("./Data/Collected_data_from_reports_interpolated.xlsx", sheet=1))
dim(time_series)


png("./Figures/Balance_data.png", height=2000, width=3000, res=300)
par(mfrow=c(2,3))
plot(colnames(time_series[i,3:31]),time_series[1,3:31]*0.45*1000, ylab = time_series[6,2], xlab="Years", main= time_series[i,1])
lines(na.approx(ts(as.numeric(time_series[1,3:31]*0.45*1000), start=1991, end=2019)), lty=2)
for(i in 2:6){
  plot(colnames(time_series[i,3:31]),time_series[i,3:31], ylab = time_series[i,2], xlab="Years", main= time_series[i,1])
  lines(na.approx(ts(as.numeric(time_series[i,3:31]), start=1991, end=2019)), lty=2)
}
dev.off()


plot(colnames(time_series[i,3:31]),time_series[1,3:31]*0.45*1000, ylab = time_series[6,2], xlab="Years", main= time_series[i,1])
lines(na.approx(ts(as.numeric(time_series[1,3:31]*0.45*1000), start=1991, end=2019)), lty=2)


par(mar=c(4,5,4,5))
plot(colnames(time_series[1,3:31]),time_series[13,3:31]*time_series[12,3:31],ylab = "ug m2", xlab="Years", main="Ca++ in the litterfall", type="l", lwd=2)
par(new = TRUE)
plot(colnames(time_series[1,3:31]),time_series[13,3:31]*time_series[12,3:31], type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", pch=NA, ylim=c(0, max(as.numeric(time_series[4,3:31]), na.rm=T)))
lines(names(na.approx(time_series[4,3:31])),na.approx(as.numeric(time_series[4,3:31])),  col="red", lwd=2, lty=2)
axis(side=4, at = pretty(range(na.approx(as.numeric(time_series[4,3:31])))), col="red")
mtext(time_series[4,2], side=4, line=3, col="red")
lines(colnames(time_series[1,3:31]),time_series[4,3:31])





### Fit an exponential to the decay, utilizing ICBM young pool as reference for decomposition

#exponential function (in percent of materail decayed)
exp_decay_percent<-function(time, r){
  k=0.400297610547104 #ICBM k1 recalibration
  100*exp(-k*time*r)
}

#cost function
exp_decay_fit<-function(x){
  ts_sim<-exp_decay_percent(seq(1:3), x)
  return(rmse(ts, ts_sim))
}

#fitting the cost function to the vector
r_vec<-c()
for(i in 9:31){
  ts<-time_series[7:9,i]
  r_vec[i]<-optimize(exp_decay_fit, interval=c(0,10))$minimum
}

png("./Figures/Decomposition_scaling.png", height=1500, width=2000, res=300)
plot(colnames(time_series[i,3:31]),r_vec[3:31], xlab="Years",  ylab="Decomposition scaling", pch=17, main="r scaling factor (assuming ICBM2022 kinetic of decay)")
dev.off()

time_series<-rbind(time_series, r_vec)
time_series[14,1]<-"r_scaling"


#creating the input table for running the model (with a bit of interpolation)
time_series_input<-rbind(time_series[14,-c(1,2)],
                         time_series[1,-c(1,2)],
                         time_series[2,-c(1,2)],
                         time_series[3,-c(1,2)],
                         time_series[4,-c(1,2)],
                         time_series[5,-c(1,2)],
                         time_series[6,-c(1,2)],
                         time_series[10,-c(1,2)],
                         time_series[13,-c(1,2)])

time_series_input=as.data.frame(t(time_series_input))
time_series_input
colnames(time_series_input)<-time_series[,1][c(14,1:6,10,13)]

time_series_input$CWD[time_series_input$CWD<0]=0

time_series_input$`Ca++ Litter`<-as.numeric(time_series[13,-c(1,2)]) #substituting the CA++ in the litterfall with Stefan's values

time_series_input_interp<-time_series_input
time_series_input_interp$Biomass[1:26]<-na.approx(time_series_input_interp$Biomass)
time_series_input_interp$Biomass[27:29]=time_series_input_interp$Biomass[26]

#adding climate scaling before 1996
time_series_input_interp$r_scaling[1:6]<-mean(time_series_input_interp$r_scaling[7:10])

#interpolating deposition based on averages
time_series_input_interp$`Ca++ through`[is.na(time_series_input_interp$`Ca++ through`)]=mean(time_series_input_interp$`Ca++ through`[!is.na(time_series_input_interp$`Ca++ through`)])

# #assuming litterfall time series
# C_litterfall_interp<-time_series_input_interp$`C litterfall`
# C_litterfall_interp[which(!is.na(C_litterfall_interp))[1]:length(C_litterfall_interp)]<-na.approx(C_litterfall_interp)
# C_litterfall_interp[is.na(C_litterfall_interp)]<-mean(time_series_input_interp$`C litterfall`, na.rm=T)
# time_series_input_interp<-cbind(time_series_input_interp, C_litterfall_interp)

#assuming CWD time series
C_lCWD_interp<-diff(time_series_input_interp$Biomass)
C_lCWD_interp<-c(C_lCWD_interp[1],C_lCWD_interp)
C_lCWD_interp[!C_lCWD_interp<0]=NA
C_lCWD_interp<-(-C_lCWD_interp)
C_lCWD_interp[is.na(C_lCWD_interp)]<-0.2 #assuming a minimal CWD input of 0.2 t ha y-1

#convert mEq m2 y-1 to Mg ha y-1
#mEq to mg =20
#mg to kg
#m2 to ha
meq2Mg<-function(x){
  res=(x*40.078*0.5)*10^-9*10^4
  return(res)
}



time_series_input_interp$`Ca++ through`
meq2Mg(x=time_series_input_interp$`Ca++ through`)

((time_series_input_interp$`C litterfall`)/1000)/

plot(meq2Mg(time_series_input_interp$`Ca++ stream`))
unique(meq2Mg(time_series_input_interp$`Ca++ stream`))




#### using Q, assuming tmax

plot(time_series_input_interp$`Littertrap C`)

plot(meq2Mg(time_series_input_interp$`Ca++ Litter`), type="l", ylim=c(0,0.08))
lines(meq2Mg(time_series_input_interp$`Ca++ stream`))

par(mfrow=c(2,1))
plot(rownames(time_series_input_interp),time_series_input_interp$Biomass, type="l", ylab = "C Mg ha-1 y-1", xlab="year", col="green", main="Biomass Carbon", ylim=c(0,225))
plot(rownames(time_series_input_interp),time_series_input_interp$CWD, type="l", ylab = "C Mg ha-1 y-1", xlab="year", col="red", main="Litter and CWD Carbon", lty=2)
lines(rownames(time_series_input_interp),time_series_input_interp$`Littertrap C`*0.01, col="blue", lty=3)
legend("topleft", c("CWD", "Litter"), col=c("red", "blue"), lty=c(2,3), bty="n")

source("Q_model_singlerun.R")


####!!! here I am assuming the Ca/C ratio in the litter is the measured one, it is very far from Hyvönen paper!!!
mean_litter_ratio<-mean(as.numeric(time_series[14,-c(1,2)]), na.rm=T)

sim_Q<-Q_Ca(
  #variables:
  biomass=time_series_input_interp$Biomass, #biomass vector, annual values
  I_F_C=time_series_input_interp$`Littertrap C`*0.01,#annual C inputs, litterfacll, in g/m2 converted to Mg/ha
  I_D_C=time_series_input_interp$CWD,#annual C inputs, CWD, 
  ## parameters
  rho_alive=1/(8.3*10^3),
  rho_dead=1/mean_litter_ratio,
  P_Ca=meq2Mg(time_series_input_interp$`Ca++ through`), #net calcium deposition every year
  spinup=2000,
  e_tmax=0.6,
  delay=0) 



par(mfrow=c(4,1))
plot(sim_Q$litter_Ca, type="l", ylim=c(0, max(c(meq2Mg(time_series_input_interp$`Ca++ Litter`), sim_Q$litter_Ca), na.rm=T)), ylab="Ca flux in the litterfall (Mg ha-1 y-1)" , xlab="Year")
points(meq2Mg(time_series_input_interp$`Ca++ Litter`))


plot(time_series_input_interp$CWD, type="l", col="green", yaxt="n", ann=F, lty=2)
legend("topleft", c("Q simulation", "measurements", "Assumed CWD input"), lty=c(1,NA,2), pch=c(NA,1,NA), col=c("cadetblue3", "black", "green"), bty="n")

plot(sim_Q$litter_C, type="l", ylab="C stored in the litter (litterfall+CWD) (Mg ha-1)", col="brown", xlab="Year")

plot(time_series_input_interp$Biomass, type="l", ylab="Aboveground biomass  (Mg ha-1)", col="brown", xlab="Year", lty=2)
points(as.numeric(time_series[1,3:29]))


png("./Figures/Simulation_draft.png", height = 2000, width = 2800, res=300)
plot(sim_Q$litter_Ca, type="l", ylab="Potential Ca flux in the litterfall (Mg ha-1 y-1)" , xlab="Year", col="firebrick", lwd=2, ylim=c(0,0.06))
points(meq2Mg(time_series_input_interp$`Ca++ Litter`), pch=16, col="red")
lines(sim_Q$decomp_Ca, type="l", ylim=c(0, max(c(meq2Mg(time_series_input_interp$`Ca++ Litter`), sim_Q$litter_Ca), na.rm=T)), col="darkorchid", lwd=2, lty=2)
legend("topleft", c("Litterfall (predictions)", "Litterfall (measured)", "CWD (predictions)"), lty=c(1,NA,2), pch=c(NA, 16, NA), bty="n", lwd=c(2,NA,2), col=c("firebrick", "red", "darkorchid"))
dev.off()





##### Extending the simulation with uncertainty

general_parms<-as.data.frame(read.csv("./General_parameters.csv"))
local_parms<-as.data.frame(read.csv("./Local_parameters.csv"))

posteriors<-as.data.frame(read.csv("./Data/posteriors_general_model.csv"))
dim(posteriors)




beta_min=local_parms[1,]$Min
beta_max=local_parms[1,]$Max
eta_11_min=local_parms[2,]$Min
eta_11_max=local_parms[2,]$Max
e0_min=local_parms[7,]$Min
e0_max=local_parms[7,]$Max
fc_min=local_parms[8,]$Max
fc_max=local_parms[8,]$Max
u0_min=(0.0855+0.0157*(50.6-0.768*57.12067))*general_parms[8,]$min
u0_max=(0.0855+0.0157*(50.6-0.768*57.12067))*general_parms[8,]$max
q0_min=min(local_parms[c(3:6),]$Min)
q0_max=max(local_parms[c(3:6),]$Max)
tmax_min= 1.06*general_parms[7,]$min #tmax comes from the average of Mäkinen dataset
tmax_max= 1.06*general_parms[7,]$max
delay_min=min(local_parms[c(9:10),]$Min)
delay_max=max(local_parms[c(9:10),]$Max)

colnames(posteriors)



###multirun here

runs=10
name_vec<-c()
for(i in 1:runs){name_vec[i]<-paste("run", i)}

params_table<-mat.or.vec(8, runs)
rownames(params_table)<-c("beta","eta11", "e0", "fc", "u0", "q0", "tmax", "delay")
colnames(params_table)<-name_vec
params_table[1,]<-runif(runs, min=beta_min, max=beta_max)
params_table[2,]<-runif(runs, min=eta_11_min, max=eta_11_max)
params_table[3,]<-runif(runs, min=e0_min, max=e0_max)
params_table[4,]<-runif(runs, min=fc_min, max=fc_max)
params_table[5,]<-runif(runs, min=u0_min, max=u0_max)
params_table[6,]<-runif(runs, min=q0_min, max=q0_max)
params_table[7,]<-runif(runs, min=tmax_min, max=tmax_max)
params_table[8,]<-runif(runs, min=delay_min, max=delay_max)



runs=dim(posteriors)[1]

params_table<-mat.or.vec(8, dim(posteriors)[1])
rownames(params_table)<-c("beta","eta11", "e0", "fc", "u0", "q0", "tmax", "delay")
params_table[1,]<-posteriors$beta
params_table[2,]<-posteriors$eta_11
params_table[3,]<-posteriors$e0
params_table[4,]<-posteriors$fc
params_table[5,]<-1.13*posteriors$u0_error
params_table[6,]<-posteriors$q0
params_table[7,]<-1.06*posteriors$tmax_error
params_table[8,]<-posteriors$delay

colMeans(posteriors)


outflux_Ca_multirun<-mat.or.vec(length(sim_Q$outflux_Ca), runs)
litter_Ca_multirun<-mat.or.vec(length(sim_Q$litter_Ca), runs)
decomp_Ca_multirun<-mat.or.vec(length(sim_Q$decomp_Ca), runs)
litter_C_multirun<-mat.or.vec(length(sim_Q$litter_C), runs)

source("Q_model_multirun.R")


for (k in 1:runs){

  # sim_Q_single<-Q_Ca(
  # #variables:
  # biomass=time_series_input_interp$Biomass, #biomass vector, annual values
  # I_F_C=time_series_input_interp$`Littertrap C`*0.01,#annual C inputs, litterfacll, in g/m2 converted to Mg/ha
  # I_D_C=time_series_input_interp$CWD,#annual C inputs, CWD, 
  # ## parameters
  # rho_alive=1/(8.3*10^3),
  # rho_dead=1/mean_litter_ratio,
  # P_Ca=meq2Mg(time_series_input_interp$`Ca++ through`), #net calcium deposition every year
  # spinup=100,
  # e_tmax=0.6,
  # delay=0) 

  sim_Q_single<-Q_Ca_multi(
    #variables:
    biomass=time_series_input_interp$Biomass, #biomass vector, annual values
    I_F_C=time_series_input_interp$`Littertrap C`*0.01,#annual C inputs, litterfacll, in g/m2 converted to Mg/ha
    I_D_C=time_series_input_interp$CWD,#annual C inputs, CWD, 
    ## parameters
    rho_alive=1/(8.3*10^3),
    rho_dead=1/mean_litter_ratio,
    P_Ca=meq2Mg(time_series_input_interp$`Ca++ through`), #net calcium deposition every year
    spinup=100, 
    beta=params_table[1,k],
    eta_11=params_table[2,k],
    e0=params_table[3,k],
    fc=params_table[4,k],
    u0=params_table[5,k],
    q0=params_table[6,k],
    tmax=params_table[7,k],
    delay=params_table[8,k]) 

  outflux_Ca_multirun[,k]<-sim_Q_single$outflux_Ca
  litter_Ca_multirun[,k]<-sim_Q_single$litter_Ca
  decomp_Ca_multirun[,k]<-sim_Q_single$decomp_Ca
  litter_C_multirun[,k]<-sim_Q_single$litter_C
  
}




colnames(outflux_Ca_multirun)<-name_vec
colnames(litter_Ca_multirun)<-name_vec
colnames(decomp_Ca_multirun)<-name_vec
colnames(litter_C_multirun)<-name_vec

years_sim<-seq(1:length(sim_Q$outflux_Ca))

rowMeans(litter_Ca_multirun)

rowQuantiles(litter_Ca_multirun, probs=c(0.025, 0.975))

library(modeest)

png("./Figures/Simulation_range.png", height = 3500, width = 2800, res=300)
par(mfrow=c(2,1))
plot(years_sim, apply(litter_Ca_multirun, 1, FUN=mlv),  type="l", ylab="Potential Ca flux in the litterfall (Mg ha-1 y-1)" , xlab="Year", col="firebrick", lwd=2, ylim=c(0,max(litter_Ca_multirun)*1.1))
points(meq2Mg(time_series_input_interp$`Ca++ Litter`), pch=16, col="red")
polygon(c(years_sim, rev(years_sim)), c(rowQuantiles(litter_Ca_multirun, probs=c(0.025, 0.975))[,1], rev(rowQuantiles(litter_Ca_multirun, probs=c(0.025, 0.975))[,2])), col=add.alpha("darkorange", 0.2))
legend("topleft", c("Litterfall (mode of predictions)", "Litterfall (measured)"), lty=c(1,NA,2), pch=c(NA, 16, NA), bty="n", lwd=c(2,NA,2), col=c("firebrick", "red", "darkorchid"))
legend("topright", "uncertainty", pch=15, col=add.alpha("darkorange", 0.2), bty="n")

plot(years_sim, apply(decomp_Ca_multirun, 1, FUN=mlv),  type="l", ylab="Potential Ca flux in the CWD (Mg ha-1 y-1)" , xlab="Year", col="darkgreen", lwd=2, ylim=c(0,max(decomp_Ca_multirun)*1.1))
polygon(c(years_sim, rev(years_sim)), c(rowQuantiles(decomp_Ca_multirun, probs=c(0.025, 0.975))[,1], rev(rowQuantiles(decomp_Ca_multirun, probs=c(0.025, 0.975))[,2])), col=add.alpha("chartreuse3", 0.2))
legend("topleft", c("CWD (mode of predictions)"), lty=c(1), pch=c(NA), bty="n", lwd=c(2), col=c("darkgreen"))
legend("topright", "uncertainty", pch=15, col=add.alpha("chartreuse3", 0.2), bty="n")

dev.off()






########### Temperature and precipitation future scenarios
#https://www.smhi.se/en/climate/future-climate/climate-scenarios/sweden/nation/rcp45/year/temperature

#function from Berg, B., Johansson, M.B., Liu, C., Faituri, M., Sanborn, P., Vesterdal, L., Ni, X., Hansen, K., Ukonmaanaho, L., 2017. Calcium in decomposing foliar litter – A synthesis for boreal and temperate coniferous forests. For. Ecol. Manage. 403, 137–144. https://doi.org/10.1016/j.foreco.2017.08.022
spruce_CA<-function(MAP){
  Ca=60.816*exp(-0.003*MAP)
  return(Ca)
}
#the paper contains functions for many other species too!

spruce_CA(750)


