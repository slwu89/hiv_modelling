################################################################################
# 
#   Testing a basic first-order upwind scheme for McKendrick-Von Foerster equation
#   see https://en.wikipedia.org/wiki/Von_Foerster_equation for details.
#   Using implementation described in:
#   Abia, L. M., Angulo, O., & Lòpez-Marcos, J. C. (2002). Numerical Solution of Stuctured Population Models, (May 2014).
# 
################################################################################

rm(list=ls());gc()

# maximum size of mesh (all time units in years)
T_max <- 100
A_max <- 80

# size of mesh cells for finite differences method
dt <- 7/365
da <- 7/365

# number of cells in mesh
# t_cells_n <- T_max/dt
# a_cells_n <- A_max/da

t_cells <- seq(from=0,to=T_max,by=dt)
t_cells_n <- length(t_cells)

a_cells <- seq(from=0,to=A_max,by=da)
a_cells_n <- length(a_cells)

# the solution mesh
U <- matrix(0,nrow = t_cells_n,ncol = a_cells_n,dimnames = list(
  as.character(round(t_cells,digits = 4)),
  as.character(round(a_cells,digits = 4))
))

# initial density
n <- 1e4
u0 <- dnorm(x = a_cells,mean = 20,sd = 5)
# u0 <- dexp(x = a_cells,rate = 1/52)
# u0 <- dlnorm(x = a_cells,meanlog = log(10),sdlog = log(5))
u0 <- (u0*n)/sum(u0)

U[1,] <- u0
v <- dt/da

# # fertility function
# beta <- function(U,t){
#   N <- sum(U[t,])
#   return((1/52)*N)
# }
# 
# # mortality function
# mu <- function(a,t){
#   return(1/52)
# }
# 
# # run simulation
# pb <- txtProgressBar(min = 1,max = (t_cells_n-1))
# for(t in 1:(t_cells_n-1)){
#   
#   a <- 1 # for births
#   U[t+1,a] <- beta(U,t)*dt
#   
#   for(a in 2:a_cells_n){
#     
#     U[t+1,a] <- U[t,a] - (v*(U[t,a] - U[t,(a-1)])) - (dt*mu(a,t)*U[t,a])
#     
#   }
#   
#   setTxtProgressBar(pb,t)
# }

Rcpp::sourceCpp('Desktop/git/hiv_modelling/age-structured/basic_upwind.cpp')

vonFoerster_upwind(U,dt,da)

library(reshape2)
library(ggplot2)
library(viridis)
U_gg <- melt(U,varnames = c("time","age"))

ggplot() +
  geom_raster(aes(x=time,y=age,fill=value),data=U_gg[U_gg$time<10,]) +
  scale_fill_viridis() +
  theme_bw() +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
  panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
  panel.grid=element_blank(),strip.background = element_blank())
