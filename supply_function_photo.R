library(plantecophys)

# seq_Ci <- seq(5,600, length=101)
# aci <- AciC4(seq_Ci)
# #plot(seq_Ci, aci$ALEAF)
# with(aci, plot(Ci, ALEAF, type='l', ylim=c(0,max(ALEAF))))

p50 <- 1.8
p_crit <- 10 # as in Sperry

# k_leaf <- function(p50) {((61.922 + (7.758 * p50)) * (18 * 0.001)) / 1000} from aDGVM2 (kg m^-1 s^-1 MPa^-1)
k_leaf <- function(p50) {(61.922 + (7.758 * p50 * -1))/0.2} # (mmol m^-2 s^-1 MPa^-1) (0.2 assumes leaf lenfth of 20cm)

K_max <- k_leaf(p50)
res <- 1/K_max
conductance1 <- function(x) (K_max*(1 - (1 / (1 + exp(2.0*(p50 + x)))))) # (mmol m^-2 s^-1 MPa^-1) I don't get why sperry doesn't have the difference between soil and leaf----
test <- seq(0, 8, length=1000)

gs_seq_test <- seq( ((0.000001 * (18 * 0.001))/(1000*0.2))*43200, ((0.1* (18 * 0.001))/(1000*0.2))*43200, length=1000) # stomatal conductance to co2 (mol m^-2 s^-1) (Medlyn 2011)

gs_seq <- seq(0.0001, 0.5, length=1000) # stomatal conductance to co2 (mol m^-2 s^-1) (Medlyn 2011)

gctogw <- 1.57  # conversion
Ca <- 400
gc <- gs / gctogw  # stomatal conductance to CO2

Evapo <- rep(0, length=1000) # this should be (mmol m^-2 s^-1) (ref. plantecophys package)
psi_leaf <- rep(0, length=1000)
psi_leaf_x <- rep(0, length=1000)
psi_leaf_delta <- rep(0, length=1000)
psi_soil <- 0
As <- rep(0, length=1000)
GS_out <- rep(0, length=1000)
slope_As <- rep(0, length=1000)
slope_cost <- rep(0, length=1000)  
kak_plc <- rep(0, length=1000)  
kak_plcx <- rep(0, length=1000)  

psi_leafxx <- seq(0.0, 8, length=1000)
cum_can_transportx <- matrix(0,0,nrow=1000, ncol=1) # this is the supply function
predawn_soil_mat_pot <- psi_soil # this assumes initial plant matric potential is the same as the soil matric potential
conductance1x <- function(x) (K_max*(1 - (1 / (1 + exp(2.0*(p50 - x)))))) # (mmol m^-2 s^-1 MPa^-1) I don't get why sperry doesn't have the difference between soil and leaf----

  for(i in 1:1000)
  {
    ffx <- integrate(conductance1x, predawn_soil_mat_pot, psi_leafxx[i] )
    cum_can_transportx[i,] <- pmax(0, ffx$value)
  }


for(i in  1:1000)
{
  
  photo1 <- Photosyn( GS=gs_seq[i], Ca=Ca )
  #aci_1 <- Aci(Ci=photo1$Ci, Ca=Ca) # A-ci calculates transpiration
  Evapo[i] <- photo1$ELEAF # (mmol m^-2 s^-1) (mmol H2O m^-2 s^-1) (1000*GS*VPD/Patm) (https://bitbucket.org/remkoduursma/plantecophys/src/0a72ce787d8962d07203f4cd1427d92ff11833df/R/photosyn.R?at=master&fileviewer=file-view-default)
  As[i] <- photo1$ALEAF
  GS_out[i] <- photo1$GS  # mol m^-2 s^-1 (Medlyn et al 2011)
  
  # leaf matric potential is something like
  if(i == 1)
  {
  psi_leaf_delta[i] <-  psi_soil - (Evapo[i]*(1/conductance1(psi_leaf[i])) - psi_leaf[i]) # psi_leaf starts at 0 MPa (i.e. psi_soil)
  psi_leaf_x[i] <-  psi_soil - (Evapo[i]*(1/conductance1(psi_soil))) # psi_leaf starts at 0 MPa (i.e. psi_soil)
  }
  
  print("---------------------------------")
  print("i")
  print(i)
  print("psi_leaf[i-1]")
  print(psi_leaf[i-1])
  print(psi_leaf_x[i-1])
  
  if(i > 1)
  {
    
  print("here")  
  psi_leaf_delta[i] <-  psi_soil - (Evapo[i]*(1/conductance1( psi_leaf[i-1])) - psi_leaf[i-1]) 
  psi_leaf_x[i] <-  psi_soil - (Evapo[i]*(1/conductance1( psi_leaf_x[i-1]))) 
  print("here2")  
  print("psi_leaf_x")
  print(psi_leaf_x[i])
  
  psi_leaf[i] <- psi_leaf[i-1] + psi_leaf_delta[i]
  slope_As[i] <- ( As[i] - As[i-1] ) / ( GS_out[i] - GS_out[i-1] )
  kak_plc[i] <- 1 - (conductance1(psi_leaf[i]) / conductance1(0)) 
  kak_plcx[i] <- 1 - (conductance1(psi_leaf_x[i]) / conductance1(0)) 
  slope_cost[i] <- ( kak_plc[i] - kak_plc[i-1] ) / ( GS_out[i] - GS_out[i-1] )
#  slope_cost[i] <- 1 - ( (conductance1(p50, 0) - conductance1(p50, psi_leaf[i]))/conductance1(p50, 0) ) # wrong
  }

#  slope_cost[i] <- 1 - ( (conductance1(p50, 0) - conductance1(p50, psi_leaf[i]))/conductance1(p50, 0) ) # wrong
#  slope_cost[i] <- ( (conductance1(p50, 0) - conductance1(p50, psi_leaf[i]))/conductance1(p50, 0) ) # wrong
  #  slope_cost[i] <- conductance1(p50, psi_leaf[i])
#  kak_plc[i] <- 1 - conductance1(p50, psi_leaf[i]) /conductance1(p50, 0) 
#  slope_cost[i] <- kak_plc[]
  
}





# NOTE: I can't find an easy way to calculate leaf temperaute easily. We have it in ADGVM though. 
# plot(photo1$GS, photo1$ALEAF)
plot(GS_out, Evapo) # Evapo (mmol m^-2 s^-1), GS (mol m^-2 s^-1) 
plot(Evapo,GS_out) # Evapo (mmol m^-2 s^-1), GS (mol m^-2 s^-1) 
plot(GS_out[1:1000], psi_leaf[1:1000], ylab=c("psi_leaf (MPa)"), xlab=("Stomatal conductance (check units)"))
plot(Evapo, psi_leaf[1:1000], ylab=c("psi_leaf (MPa)"), xlab=("Transpiraiton"))
plot(Evapo, psi_leaf_delta[1:1000], ylab=c("psi_leafdelta (MPa)"), xlab=("Transpiraiton"))
plot(Evapo, psi_leaf_x[1:1000], ylab=c("psi_leaf (MPa)"), xlab=("Transpiraiton"))
plot(GS_out[1:1000], psi_leaf_x[1:1000], ylab=c("psi_leaf (MPa)"), xlab=("Stomatal conductance"))

plot(GS_out, Evapo, ylab=c("Evapo"), xlab=("GS"))

plot(GS_out, As, ylab=c("As"), xlab=("GS"))

plot(Evapo, As, ylab=c("As"), xlab=("Evapo"))

plot(GS_out, kak_plc, ylab=c("loss of conductance "), xlab=("GS"))

plot(GS_out, kak_plcx, ylab=c("loss of conductance "), xlab=("GS"))

##compare integral transpiration to interative transpirtation for a given psi_leaf
plot(psi_leafxx, cum_can_transportx)
plot(psi_leaf_x[1:1000], Evapo, xlab=c("psi_leaf (MPa)"), ylab=("Transpiraiton"))
## VERY CLOSE - for iterative method with psi leaf = -0.032189 transp=7.5, for integral with psi leaf = 0.032 transp=7.47 

# par(yaxs="i")
# #with(p, plot(Ci, ALEAF, type='l', ylim=c(0,max(ALEAF))))
# with(photo1, points(Ci, ALEAF, pch=19, col="red"))
# abline(gc * Ca, -gc, lty=5)
# 
# 

# gs <- 0.2  # stomatal conductance to H2O
# Ca <- 400  # ambient CO2
# gctogw <- 1.57  # conversion
# gc <- gs / gctogw  # stomatal conductance to CO2
# # Demand curve (Farquhar model)
# p <- Aci(seq(60,500,length=101), Ca=400)
# # Provide stomatal conductance as input, gives intersection point.
# g <- Photosyn(GS=gs, Ca=Ca)
# # Intersection point visualized
# par(yaxs="i")
# with(p, plot(Ci, ALEAF, type='l', ylim=c(0,max(ALEAF))))
# with(g, points(Ci, ALEAF, pch=19, col="red"))
# abline(gc * Ca, -gc, lty=5)
# legend("topleft", c(expression("Demand:"~~A==f(C[i])),
#                     expression("Supply:"~~A==g[c]*(C[a]-C[i])),
#                     "Operating point"),
#        lty=c(1,5,-1),pch=c(-1,-1,19),
#        col=c("black","black","red"),
#        bty='n', cex=0.9)









  
e_can <- function(psi_leaf, psi_soil, res, gi) { ( psi_leaf - psi_soil )*gi / res }
gi <- function(psi_leaf, p50) { 1 - (1 / (1 + exp(5.0*(p50 - psi_leaf))))  }

p50 <- 3
res <- 10
psi_leaf <- seq(0.01, 7, length=1000)
psi_soil <- 0.0

sum_e <- (e_can(psi_leaf, psi_soil, res, gi(psi_leaf, p50) ))
new <- rep(0,length=length(sum_e))

for(i in 1:length(sum_e))
  {
    new[i] <- sum(sum_e[1:i])
  }

#e_can <- function(x, res, gi) { ( psi_leaf )*gi / res }
#gi <- function(psi_leaf, p50) { 1 - (1 / (1 + exp(2.0*(p50 - psi_leaf))))  }

## conductance here is assuming the soil matric potential is 0. 
conductance <- function(x) (8*(1 - (1 / (1 + exp(2.0*(3 - x))))))
conductance1 <- function(x) (1*(1 - (1 / (1 + exp(2.0*(3 - x))))))

plot(psi_leaf, pmax(0.0, (conductance(0) - conductance(psi_leaf))/conductance(0)) )
lines(psi_leaf, pmax(0.0, (conductance1(0) - conductance1(psi_leaf))/conductance1(0)), col="red", lwd=2 )

p50 <- 0.5
K_max <- 8
res <- 1/K_max
psi_leaf <- seq(0.0, 8, length=1000)
psi_soil <- 0.0
cum_can_transport <- rep(0, length=1000)

## transpiration rate assuming the soil matric potential is 0. 
e_can <- function(psi_leaf) { ( psi_leaf - psi_soil )*((1 - (1 / (1 + exp(2.0*(3 - psi_leaf)))))) / res }
sum_e <- e_can(psi_leaf)

##------------------------------------------------------------------------------------------------------------------------------------------
## Sperry's version of supply-demand-loss theory 
## trying to work out why Sperry doensn't have the potential difference between soil and leaf included when calculating evaporation

sperry_cond <- function(psi_leaf) { ((1 - (1 / (1 + exp(3.0*(2.5 - psi_leaf)))))) / res }
#cum_can_transportx <- rep(0, length=1000)
cum_can_transportx <- matrix(0,0,nrow=1000, ncol=5) # this is the supply function
predawn_soil_mat_pot <- seq(0,2, length=5) # this assumes initial plant matric potential is the same as the soil matric potential

for(j in 1:length(predawn_soil_mat_pot))
  {
    for(i in 1:1000)
      {
        ffx <- integrate(sperry_cond, predawn_soil_mat_pot[j], psi_leaf[i] )
        cum_can_transportx[i,j] <- pmax(0, ffx$value)
      }
  }

for(l in 1:length(predawn_soil_mat_pot)) # remove the zeros from the 
{
  kop <- which(cum_can_transportx[,l] == 0)
  cum_can_transportx[kop,l] <- NA
}

cond_max_slope_sperry <- rep(0, lenght=5)
  
for(i in 1:length(predawn_soil_mat_pot))
  {
    cond_max_slope_sperry[i] <-  sperry_cond(predawn_soil_mat_pot[i]) # for Sperry the maximum conductance is always the pre-dawn matric potential/soil matric potential 
  }

##------------------------------------------------------------------------------------------------------------------------------------------
## make transpiration demand
Gmax <- 12563.1 # (Sperry 2016, 2130 kg h^-1 m^-2) NOTE it should be (kg h^-1 MPa^-1 m^-2) (12563.1 in Excel doc)
G <- rep(Gmax, length=5)
VPD <- 1.5*0.001   # (Sperry 2016, leaf-to-air vapor pressure deficit 1 kPa) (0.001 converts to MPa) NOTE in Sperry 2016 this is 1 which is crazy, in the Excel script it is 100
#evap_demand <- G*VPD

# need to define an Pcrit, i.e. a matric potential we choose wheere we decide conductance is effectively zero. We use this to get Ecrit, i.e. maximum transpiration beyond which leads to runaway cavitation  
P_crit <- 4 # MPa - this is arbitrary and could be a plant trait 
# get maximum transpiration possible based on Pcrit and soil matric potential 
E_crit <- rep(0, lenght=5)
evap_demand <- rep(0, lenght=5)

for(j in 1:5)
  {
    ffx <- integrate(sperry_cond, predawn_soil_mat_pot[j], P_crit )
    E_crit[j] <- pmax(0, ffx$value)
    evap_demand[j] <- G[j]*VPD
    if(evap_demand[j] > E_crit[j]) evap_demand[j] <- E_crit[j] # demand can't be greater than maximum supply
  }

## quick and dirty method to find the psi where demand is met. 
demand_met_at_sperry <- rep(0, length=5)
psi_demand_met_at_sperry <- rep(0, length=5)
demand_met_at_slope_sperry <- rep(0, length=5)
loss_function_sperry <- rep(0, length=5)
regulated_transpiration <- rep(0, length=5)
regulated_leaf_psi <- rep(0, length=5)
regulated_leaf_psi1 <- rep(0, length=5)

for(i in 1:5)
  {
    demand_met_at_sperry[i] <- which(cum_can_transportx[,i] >= evap_demand[i])
    psi_demand_met_at_sperry[i] <- psi_leaf[demand_met_at_sperry[i]]
    demand_met_at_slope_sperry[i] <- sperry_cond(psi_demand_met_at_sperry[i])
    loss_function_sperry[i] <- demand_met_at_slope_sperry[i] / cond_max_slope_sperry 
    regulated_leaf_psi[i] <- predawn_soil_mat_pot[i] + ((psi_demand_met_at_sperry[i] - predawn_soil_mat_pot[i])*loss_function_sperry[i])
    regulated_leaf_psi1[i] <- ((psi_demand_met_at_sperry[i] - predawn_soil_mat_pot[i])*loss_function_sperry[i])
    
    ffx <- integrate(sperry_cond, predawn_soil_mat_pot[i], regulated_leaf_psi[i] )
    regulated_transpiration[i] <- pmax(0, ffx$value)
    G[i] <- G[i]*loss_function_sperry[i]
    
  }


demand_met_at_sperry1 <- rep(0, length=5)

for(i in 1:5)
{
  demand_met_at_sperry1[i] <- which(cum_can_transportx[,i] >= regulated_transpiration[i])
  demand_met_at_sperry1[i] <- psi_leaf[demand_met_at_sperry1[i]]
}

##--------------------------------------------------------------------------------
# Calculation of regulated stomatal conductance - these two are equivalent
#G_new <- demand_new/VPD
#G <- G*loss_function

##------------------------------------------------------------------------------------------------------------------------------------------

par(mfrow=c(2,2))
plot(psi_leaf, sperry_cond(psi_leaf), cex=0.1, main="Sperry", ylab="Conductance") # Sperry's conductance rate (Sperry and Love 2015 Fig. 1)
plot(psi_leaf, conductance(psi_leaf), cex=0.1, main="Liam", ylab="Conductance") # My conductance rate assuming both transpiraiton and conductance are affected by the matric potential difference

#graphics.off()
#pdf( file="Supply_demand_loss_Sperrys_way.pdf", width=14,height=8 )
#par( mfcol=c(2,2), mar=c(5,5,4,1))

par(mfrow=c(1,1))
plot(psi_leaf, cum_can_transportx[,1], type="l", lwd=2, cex=0.1, main="Sperry", ylab="Transport Rate") # Sperry's conductance rate (Sperry and Love 2015 Fig. 1)
points(psi_demand_met_at_sperry[1], evap_demand[1], type="p", cex=0.9, pch=19, col="black")
points(regulated_leaf_psi[1], regulated_transpiration[1], type="p", cex=1.6, pch=21, col="red", lwd=2)
lines(psi_leaf, cum_can_transportx[,2], lwd=2)
points(psi_demand_met_at_sperry[2], evap_demand[2], type="p", cex=0.9, pch=19, col="black")
points(regulated_leaf_psi[2], regulated_transpiration[2], type="p", cex=1.6, pch=21, col="red", lwd=2)
lines(psi_leaf, cum_can_transportx[,3], lwd=2)
points(psi_demand_met_at_sperry[3], evap_demand[3], type="p", cex=0.9, pch=19, col="black")
points(regulated_leaf_psi[3], regulated_transpiration[3], type="p", cex=1.6, pch=21, col="red", lwd=2)
lines(psi_leaf, cum_can_transportx[,4], lwd=2)
points(psi_demand_met_at_sperry[4], evap_demand[4], type="p", cex=0.9, pch=19, col="black")
points(regulated_leaf_psi[4], regulated_transpiration[4], type="p", cex=1.6, pch=21, col="red", lwd=2)
lines(psi_leaf, cum_can_transportx[,5], lwd=2)
points(psi_demand_met_at_sperry[5], evap_demand[5], type="p", cex=0.9, pch=19, col="black")
points(regulated_leaf_psi[5], regulated_transpiration[5], type="p", cex=1.6, pch=21, col="red", lwd=2)
legend("topleft", c("unregulated leaf psi","regulated leaf psi","p50=2.5"), col=c("black", "red", "white"), pch = 1, cex=1.5, lwd=2)
#graphics.off()


plot(psi_leaf, conductance(psi_leaf), cex=0.1, main="Liam", ylab="Conductance") # My conductance rate assuming both transpiraiton and conductance are affected by the matric potential difference

##------------------------------------------------------------------------------------------------------------------------------------------
#sperry_cond <- function(psi_leaf) { ( psi_leaf )*((1 - (1 / (1 + exp(2.0*(3 - psi_leaf)))))) / res }
G <- matrix(Gmax, Gmax, ,nrow=1000, ncol=5)

#demandx <- rep(0, length=1000)
psi_1 <- matrix(0,0,nrow=1000, ncol=5)#
slope_demand <- matrix(0,0,nrow=1000, ncol=5)
loss_fun_sp_gs <- matrix(0,0,nrow=1000, ncol=5)
reg_leaf_psi <- matrix(0,0,nrow=1000, ncol=5)

regulated_trans <- matrix(0,0,nrow=1000, ncol=5)
regulated_Gs <- matrix(0,0,nrow=1000, ncol=5)

regulated_trans_trap <- matrix(0,0,nrow=1000, ncol=5)
regulated_Gs_trap <- matrix(0,0,nrow=1000, ncol=5)

non_regulated_trans <- matrix(0,0,nrow=1000, ncol=5)
non_regulated_Gs <- matrix(0,0,nrow=1000, ncol=5)

#max_slo_sp <- sperry_cond(0) 

max_slo_sp <- rep(0, length=5)

for(j in 1:5) 
{
  max_slo_sp[j] <- sperry_cond(predawn_soil_mat_pot[j])
}

for(j in 1:5)
{
  for(i in 1:1000)
  {

  #  demandx[i] <- cum_can_transportx[i,1]
   psi_1[,j] <- seq(0, max(psi_leaf), length=1000) #psi_leaf[i]
   slope_demand[i,j] <- sperry_cond(psi_1[i,j])
  
   loss_fun_sp_gs[i,j] <- slope_demand[i,j] / max_slo_sp[j] 
   #reg_leaf_psi[i,j] <- predawn_soil_mat_pot[j] + ((psi_1[i,j] - predawn_soil_mat_pot[j])*loss_fun_sp_gs[i,j])
   if(i==1) reg_leaf_psi[i,j] <- pmax(0, ((psi_1[i,j] - predawn_soil_mat_pot[j])*loss_fun_sp_gs[i,j]))
      
   if(i>1)
   {
     reg_leaf_psi[i,j] <- pmax(0, pmax(reg_leaf_psi[i-1,j], ((psi_1[i,j] - predawn_soil_mat_pot[j])*loss_fun_sp_gs[i,j])))
   }
    
   ffx <- integrate(sperry_cond, predawn_soil_mat_pot[j], psi_1[i,j] )
   non_regulated_trans[i,j] <- pmax(0, ffx$value)
   non_regulated_Gs[i,j] <- non_regulated_trans[i,j]/0.001 # E = G*VPD ---- G = E/VPD (VPD=1, 0.001 transforms to MPa)
  # 
   ffx <- integrate(sperry_cond, predawn_soil_mat_pot[j], predawn_soil_mat_pot[j] + reg_leaf_psi[i,j] )
   regulated_trans[i,j] <- pmax(0, ffx$value)
   regulated_Gs[i,j] <- regulated_trans[i,j]/0.001 # E = G*VPD ---- G = E/VPD (VPD=1, 0.001 transforms to MPa)
   G[i,j] <- G[i,j]*loss_fun_sp_gs[i,j]
  
  }
}

holder_4_max <- rep(0, length=5)

for(i in 1:5)
{
  holder_4_max[i] <- which(reg_leaf_psi[,i]==max(reg_leaf_psi[,i]))
}

kss <- which(reg_leaf_psi==max(reg_leaf_psi))
reg_leaf_psi_new <- reg_leaf_psi
reg_leaf_psi_new[kss:length(reg_leaf_psi_new)] <- max(reg_leaf_psi)

for(i in 1:1000)
{
  
ffx <- integrate(sperry_cond, 0, reg_leaf_psi_new[i] )
regulated_trans_trap[i] <- pmax(0, ffx$value)
regulated_Gs_trap[i] <- regulated_trans_trap[i]/0.001 # E = G*VPD ---- G = E/VPD (VPD=1, 0.001 transforms to MPa)
# G[i] <- G[i]*loss_function_sperry[i]
}

 plot(psi_1, (1-slope_demand/slope_demand[1])*100, type ="l", lwd=1)
 lines(psi_1, (1-regulated_Gs/non_regulated_Gs)*100, col="red")
 lines(psi_1, (1-regulated_Gs_trap/non_regulated_Gs)*100, col="red")
 lines(psi_1, (1-regulated_Gs/Gmax)*100, col="red")
 
 plot(psi_1, (1-slope_demand/slope_demand[1])*100, type ="l", lwd=1)
 lines(psi_1, (1-regulated_Gs_trap/non_regulated_Gs)*100, col="red")
 lines(psi_1, (1-regulated_Gs_trap/Gmax)*100, col="red")
 
### 
cum_can_transport_Gs <- rep(0, length=1000)
slope_sperry1 <- rep(0, length=1000)
loss_sperry1 <- rep(0, length=1000)
slope_max <- sperry_cond(0)
  

for(i in 1:1000)
{
  ff <- integrate(sperry_cond, 0, psi_leaf[i] )
  cum_can_transport_Gs[i] <- ff$value
  slope_sperry1 <- sperry_cond(psi_leaf)
  loss_sperry1 <- slope_sperry1/slope_max
}

## quick and dirty method to find the psi where demand is met. 
demand_met_at_sperry <- rep(0, length=5)
psi_demand_met_at_sperry <- rep(0, length=5)
demand_met_at_slope_sperry <- rep(0, length=5)
loss_function_sperry <- rep(0, length=5)
regulated_transpiration <- rep(0, length=5)
regulated_leaf_psi <- rep(0, length=5)

for(i in 1:5)
{
  demand_met_at_sperry[i] <- which(cum_can_transportx[,i] >= evap_demand[i])
  psi_demand_met_at_sperry[i] <- psi_leaf[demand_met_at_sperry[i]]
  demand_met_at_slope_sperry[i] <- sperry_cond(psi_demand_met_at_sperry[i])
  loss_function_sperry[i] <- demand_met_at_slope_sperry[i] / cond_max_slope_sperry 
  regulated_leaf_psi[i] <- predawn_soil_mat_pot[i] + ((psi_demand_met_at_sperry[i] - predawn_soil_mat_pot[i])*loss_function_sperry[i])
  
  ffx <- integrate(sperry_cond, predawn_soil_mat_pot[i], regulated_leaf_psi[i] )
  regulated_transpiration[i] <- pmax(0, ffx$value)
  G[i] <- G[i]*loss_function_sperry[i]
  
}




##------------------------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------------------
## Can't work out how Sperry's Gs declines faster than conductance.......

cum_can_transport <- rep(0, length=1000)

## surrogate to get derivitive of the transpiration rate
e_can2 <- expression(( x - psi_soil )*((1 - (1 / (1 + exp(2.0*(3 - x)))))) / res )
D(D(e_can2, 'x'),'x') # need to manually copy this derivitive - cant find a way to automate it.
d_de_can2 <- function(x) {(((1 - (1/(1 + exp(2 * (3 - x)))))) - (x - psi_soil) * (exp(2 * (3 - x)) * 2/(1 + exp(2 * (3 - x)))^2))/res}
##dd_ecan2 <- function(x) {-((exp(2 * (3 - x)) * 2/(1 + exp(2 * (3 - x)))^2 + ((exp(2 * (3 - x)) * 2/(1 + exp(2 * (3 - x)))^2) - (x - psi_soil) * (exp(2 * (3 - x)) * 2 * 2/(1 + exp(2 * (3 - x)))^2 - exp(2 * (3 - x)) * 2 * (2 * (exp(2 * (3 - x)) * 2 * (1 + exp(2 * (3 - x)))))/((1 + exp(2 * (3 - x)))^2)^2)))/res)}
psi_max_slope <- uniroot(d_de_can2, c(0, 5))$root # find the leaf matric potential which yields the highest transpiraiton rate

for(i in 1:1000)
{
  ff <- integrate(e_can, 0, psi_leaf[i] )
  cum_can_transport[i] <- ff$value
}

cond_max_slope <- max(conductance(psi_leaf)) #(Sperry's )

Gmax <- 12563.1 # (Sperry 2016, 2130 kg h^-1 m^-2) NOTE it should be (kg h^-1 MPa^-1 m^-2) (12563.1 in Excel doc)
G <- Gmax
VPD <- 2*0.001   # (Sperry 2016, leaf-to-air vapor pressure deficit 1 kPa) (0.001 converts to MPa) NOTE in Sperry 2016 this is 1 which is crazy, in the Excel script it is 100
evap_demand <- G*VPD
## quick and dirty method to find the psi where demand is met.
demand_met_at <- which(cum_can_transport > evap_demand)
demand_met_at <- psi_leaf[demand_met_at[1]]
demand_met_at_slope <- conductance(demand_met_at)
loss_function <- demand_met_at_slope / cond_max_slope

demand_new <- evap_demand * loss_function

##--------------------------------------------------------------------------------
# Calculation of regulated stomatal conductance - these two are equivalent
#G_new <- demand_new/VPD
G <- G*loss_function

##------------------------------------------------------------------------------------------------------------------------------------------
library(plotrix)
conductance <- function(x) (K_max*(1 - (1 / (1 + exp(2.0*(3 - x))))))
y_int_cum <- function(x){ integrate(e_can, 0, x )$value - (e_can(x)*x)} # the y-intersept for the tangent
#abline(y_int_cum(2), e_can(2)) # the slope of the tangent is the transpiration rate (e_can)

plot(psi_leaf, conductance(psi_leaf), cex=0.1)
plot(psi_leaf, e_can(psi_leaf), cex=0.1)
plot(psi_leaf, cum_can_transport, cex=0.1)
ablineclip(y_int_cum(2), e_can(2), x1=1,x2=3, lty=2, col="red", lwd=2)

par(mfrow=c(1,1))
plot(psi_leaf, cum_can_transport, cex=0.1)
ablineclip(y_int_cum(1), e_can(1), x1=0.5,x2=1.5, lty=3, col="red", lwd=4)
ablineclip(y_int_cum(psi_max_slope), e_can(psi_max_slope), x1=psi_max_slope*0.85,x2=psi_max_slope*1.15, lty=3, col="red", lwd=4)
ablineclip(y_int_cum(4.0), e_can(4.0), x1=3.4,x2=4.5, lty=3, col="red", lwd=4)

##------------------------------------------------------------------------------------------------------------------------------------------
graphics.off()
pdf( file="Supply_demand_Sperrys_way.pdf", width=14,height=8 )
par( mfcol=c(2,2), mar=c(5,5,4,1))
plot(psi_leaf, sperry_cond(psi_leaf), cex=0.1, main="Sperry", ylab="Conductance") # Sperry's conductance rate (Sperry and Love 2015 Fig. 1)
text(7, 7, "(A)")
plot(psi_leaf, cum_can_transportx, cex=0.1, main="Sperry", ylab="Transport Rate") # Sperry's transport rate (Sperry and Love 2015 Fig. 1)
text(7, 20, "(B)")
graphics.off()




##------------------------------------------------------------------------------------------------------------------------------------------
## my way
cond_max_slope <- max(conductance(psi_leaf)) #(Sperry's )

Gmax <- 12563.1 # (Sperry 2016, 2130 kg h^-1 m^-2) NOTE it should be (kg h^-1 MPa^-1 m^-2) (12563.1 in Excel doc)
G <- Gmax
VPD <- 2*0.001   # (Sperry 2016, leaf-to-air vapor pressure deficit 1 kPa) (0.001 converts to MPa) NOTE in Sperry 2016 this is 1 which is crazy, in the Excel script it is 100
evap_demand <- G*VPD
## quick and dirty method to find the psi where demand is met.
demand_met_at <- which(cum_can_transport > evap_demand)
demand_met_at <- psi_leaf[demand_met_at[1]]
demand_met_at_slope <- conductance(demand_met_at)
loss_function <- demand_met_at_slope / cond_max_slope

demand_new <- demand * loss_function

##--------------------------------------------------------------------------------


graphics.off()
pdf( file="Supply_demand_Liams_way.pdf", width=14,height=8 )
par( mfcol=c(2,2), mar=c(5,5,4,1))
plot(psi_leaf, conductance(psi_leaf), cex=0.1, main="Liam", ylab="Conductance") # My conductance rate assuming both transpiraiton and conductance are affected by the matric potential difference
text(6, 7, "(C)")
plot(psi_leaf, cum_can_transport, cex=0.1, main="Liam", ylab="Transport Rate") # My transport rate assuming both transpiraiton and conductance are affected by the matric potential difference
text(6, 7, "(D)")
points(demand_met_at, demand, type="p", cex=1.5, pch=19, col="red")
ablineclip(y_int_cum(demand_met_at), e_can(demand_met_at), x1=demand_met_at*0.85,x2=demand_met_at*1.15, lty=2, col="red", lwd=2)
plot(psi_leaf, e_can(psi_leaf), cex=0.1, main="Liam", ylab="Conductance moderated psi diff and cavitation") # My conductance rate assuming both transpiraiton and conductance are affected by the matric potential difference
text(6, 7, "(E)")
graphics.off()

# (Sperry and Love) "The maximum dE/dPcanopy [max conductance and max slope of the line in Sperry's transport rate] is at the predawn start of the curve
# (Fig. 2, dashed dE/dPmax tangent) and equals the maximum soil–canopy k"
# While this is indeed the maximum conductance the maximum flow rate of water happens at a higher matric potential, see (E).
# This implies that the maximum slope of the transport function does not happen at 0MPa but at a higher MPa. I'm not sure what is more appropriate or correct.
# I'm obviously correct and also the best at being modest.
##------------------------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------------------


transpiration <- rep(0, length=1000)
transpiration1 <- rep(0, length=1000)
#psi_leaf <- seq(0.01, 6.0, length=1000)

for(i in 1:1000)
{

  if(i==1)  x <- 0
  ff <- integrate(conductance, 0, psi_leaf[i] )
  transpiration[i] <- ff$value
  x <- x + psi_leaf[i]
}



for(i in 1:100)
{
  if(i==1)  x <- 0
  ff <- conductance(x + psi_leaf[i] )
  transpiration1[i] <- ff
  x <- x + psi_leaf[i]
}


new1 <- rep(0,length=length(transpiration1))

for(i in 1:length(transpiration1))
{
  new1[i] <- sum(transpiration1[1:i])
}

new2 <- rep(0,length=length(new1))

for(i in 1:length(new1))
{
  if(i==1)  x <- 0
  new2[i] <- (new1[i])/(psi_leaf[i])
  x <- new1[i]
}



poo <- expression( (5*(1 - (1 / (1 + exp(2.0*(3 - x)))))))

pam <- function(x) {-(5 * (exp(2 * (3 - x)) * 2/(1 + exp(2 * (3 - x)))^2))}

psi_leaf <- seq(0.01, 7, length=1000)
psi_soil <- 0
e_can <- function(psi_leaf) { ( psi_leaf - psi_soil )*((1 - (1 / (1 + exp(2.0*(3 - psi_leaf)))))) / res }

suppx <- e_can(psi_leaf)

sum_suppx <- rep(0, length(suppx))

for(i in 1:length(suppx))
{
  if(i < 2) 
  {
    sum_suppx[i] <- suppx[i]
  }
  
  if(i > 1)
  {
    sum_suppx[i] <- sum_suppx[i-1] + abs(suppx[i] - suppx[i-1]) 
  }
}

supp1 <- rep(0, length=1000)

for(i in 1:1000)
{
  ff <- integrate(e_can, 0, psi_leaf[i] )
  supp1[i] <- ff$value
}


 plot(psi_leaf, crap(psi_leaf))
 plot(psi_leaf, e_can(psi_leaf))
 abline(v=psi_leaf[335])
 plot(psi_leaf, crap(psi_leaf))
 abline(v=psi_leaf[335])
 fuck <- crap(psi_leaf)
 which(fuck==min(fuck))
 abline(v=psi_leaf[365])
 plot(psi_leaf, e_can(psi_leaf))
 plot(psi_leaf, ball(psi_leaf))
 plot(psi_leaf, e_can(psi_leaf))
 plot(psi_leaf, ball(psi_leaf))
 abline(h=0)
 plot(psi_leaf, e_can(psi_leaf))

