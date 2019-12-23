rm(list=ls())
require(deSolve)

#### Solver de Lisis determinista
X_theta_lisis<-function(theta,t=tiempos, X_ini){
  SEIRmod_Lisis <- function(t, x, theta) { 
    with(as.list(c(theta, x)), 
         { 
           ds <- alp*s*(1-s/s_m)-lam*s*s_v-mu_c*s 
           de <- lam*s*s_v-del*e-mu_c*e
           dr_i<- del*e
           dr <- mu_c*s+mu_c*e
           ds_v<- m*r_i-mu_v*s_v-lam*s*s_v
           dr_v<- mu_v*s_v+lam*s*s_v 
           res <- c(ds, de, dr_i, dr, ds_v, dr_v) 
           list(res) }
    ) }
  ## Solver
  out <- lsoda(X_ini, t, SEIRmod_Lisis, theta)
  out[which(out[,6]<0),6]<-0
  return(out)
}

ej_l<-X_theta_lisis(theta=c(alp=0.3,lam=0.002,del=0.24,m=100,mu_c=0.0395,
                    mu_v=3,s_m=1500),t=seq(0,100,by=.1),X_ini=c(s=300,e=0,r_i=0,r=0,s_v=10,r_v=0))

par(mfrow=c(2,1),mar=c(4,4,4,1))
matplot(ej_l[,1],ej_l[,2:5],main="Células: Modelo determinístico para Lisis",xlab="Tiempo (días)",ylab="Individuos",t="l",lty=c(1,1,1,1),col=c("turquoise1","limegreen","coral","darkorchid"),lwd=2)
legend("right",legend=c("S","E","R_I","R"), col=c("turquoise1","limegreen","coral","darkorchid"),lty=1,
       bg= ("white"), cex=0.4,horiz=F,bty = "o",lwd=2)

matplot(ej_l[,1],ej_l[,6:7],main="Virus: Modelo determinístico para Lisis",xlab="Tiempo (días)",ylab="Individuos",t="l",lty=c(1,1,1,1),col=c("mediumblue","khaki4"),lwd=2)
legend("right",legend=c("S_v","R_v"), col=c("mediumblue","khaki4"),lty=1,
       bg= ("white"),cex= 0.4,horiz=F,bty = "o",lwd=2)


#### Solver de Brote determinista
X_theta_brote<-function(theta,t=tiempos, X_ini){
  SEIRmod_Brote <- function(t, x, theta) { 
    with(as.list(c(theta, x)), 
         { 
           ds <- alp*s*(1-s/s_m)-lam*s*s_v-mu_c*s
           de <- lam*s*s_v-del*e-mu_c*e
           di <- del*e-mu_c*i
           dr <- mu_c*s+mu_c*e+mu_c*i
           ds_v <- gam*i-mu_v*s_v-lam*s*s_v
           dr_v <- mu_v*s_v+lam*s*s_v 
             
           res <- c(ds, de, di, dr, ds_v, dr_v) 
           list(res) }
    ) }
  ## Solver
  out <- lsoda(X_ini, t, SEIRmod_Brote, theta)
  out[which(out[,6]<0),6]<-0
  return(out)
}

ej_b<-X_theta_brote(theta=c(alp=0.3,lam=0.002,del=0.24,gam=3.95,mu_c=0.0395,
                            mu_v=3,s_m=100),t=seq(0,1500,by=.1),X_ini=c(s=300,e=0,i=0,r=0,s_v=10,r_v=0))

par(mfrow=c(2,1),mar=c(4,4,4,1))
matplot(ej_b[,1],ej_b[,2:5],main="Células: Modelo determinístico para Brote",xlab="Tiempo (días)",ylab="Individuos",t="l",lty=c(1,1,1,1),col=c("turquoise1","limegreen","coral","darkorchid"),lwd=2)
legend("right",legend=c("S","E","I","R"), col=c("turquoise1","limegreen","coral","darkorchid"),lty=1,
       bg= ("white"), cex= 0.4,horiz=F,bty = "o",lwd=2)

matplot(ej_b[,1],ej_b[,6:7],main="Virus: Modelo determinístico para Brote",xlab="Tiempo (días)",ylab="Individuos",t="l",lty=c(1,1,1,1),col=c("mediumblue","khaki4"),lwd=2)
legend("right", legend=c("S_v","R_v"), col=c("mediumblue","khaki4"),lty=1,
       bg= ("white"),cex = 0.4,horiz=F,bty = "o",lwd=2)
