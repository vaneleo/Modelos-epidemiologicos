# CODIGO PARA LAS FUNCIONES DE GILLESPIE EN REDES: LISIS Y BROTE

#ENTRADA:
# n_I es el numero de individuos de una poblacion (ya sean celulas o virus) al tiempo 0 o iniciales en el sistema.
# d es la medida del lado del espacio total E.
#SALIDA:
# Regresa una matriz Coord de dimensiones (n_I)x2 donde cada renglon k indica las coordenadas (x,y) del individuo k,
# las cuales se eligen aletoriamente dentro del espacio E 

####  YA  ####

Coloca<-function(n_I,d){
       Coord_x<-runif(n_I,0,d) 
       Coord_y<-runif(n_I,0,d)  
  return(Coord=cbind(Coord_x,Coord_y))  
}

#ENTRADA:
# Coord_v=(x,y) son las coordenadas de un virus 
# Coord_c es la matriz de coordenadas de las celulas, de dim |S|x2
# Si un virus esta a una distancia menor que r de una celula entonces la celula puede ser infectada por ese virus. 
#SALIDA:
# Regresa un vector de dim |S| cuyas entradas son TRUE o FALSE, dependiendo si el virus puede infectar a cada una de las celulas 

### YA ###
Distancia<-function(Coord_c,Coord_v,r){
  res<-0
  res<-sqrt((Coord_c[,1]-Coord_v[1])^2+(Coord_c[,2]-Coord_v[2])^2)
  res<-res<r
  return(res)  
}

#Ejemplo
#r<-Coloca(100,10)
#plot(r[,1],r[,2])
#v<-c(5,5)
#ff<-Distancia(r,v,2)
#lines(r[ff,1],r[ff,2],t="p",col="red",pch=19)
  

#ENTRADA:  
# Coord_S y Coord_S_v son matrices de dim |S|x2 y |S_v|x2 que contienen las coordenadas de S y S_v, respectivamente.
# Si un virus esta a una distancia menor que r de una celula entonces la celula puede ser infectada por ese virus.
#SALIDA:
# Regresa una matriz A de tamaño |S|x|S_v|, cuyas entradas A_ik son 0 o 1, dependiendo si la celula i esta al alcance del virus k.
# Llamaremos a cada alcance positivo conexion.
# NOTA: Tambien nos sera util calcular la suma de A que es el numero de conexiones para actualizar las poblaciones.
###  YA  ###

Mat_Ady<-function(Coord_S,Coord_S_v,r){
  s<-length(Coord_S[,1])
  s_v<-length(Coord_S_v[,1])
  A<-matrix(0,s,s_v)
    for(j in 1:s_v){
         A[,j]<-Distancia(Coord_S,Coord_S_v[j,],r) #¿Por que si la funcion Distancia es booleana A es numerica?
         }
  return(A)
}

#Ejemplo
#v<-Coloca(10,10)
AA<-Mat_Ady(r,v,2)
#numero de conexiones:
sum(AA)


#ENTRADA:
# El vector Propensidades 
#SALIDA:
#Elige el j que indica el indice de la propensidad (en el vector Propensidades) de la siguiente reaccion,
#el cual definira que tipo de reaccion sucedio y entre cuales individuos.
###   YA   ### 

Sig_Reac<-function(Propensidades){#seed=1){
  #set.seed(seed)
  r_2<-runif(1) #r_2 ~Unif(0,1)
  p<-cumsum(Propensidades)
  alp_0<-sum(Propensidades)
  res<-min(which(p>=r_2*alp_0))
  return(res) #(list(res,r_2))  
} 

#Ejemplo
#prop<-c(.1,.3,.5,.6)
#r2<-0.9
#Sig_Reac(prop)


#CONTEXTO:
# Aqui las poblaciones son: Poblaciones=(S,E,R_I,R,S_v,R_v) 
# Estamos suponiendo que cuando comenzamos el proceso el sistema tiene un numero de celulas expuestas E=0
# El vector de reacciones esta dado con el siguiente orden:
#Reac=(Nacimiento S, Muerte natural de S, Muerte natural de E, Exposicion, Infeccion E->R_I y Liberacion de S_v, Muerte de S_v) 
# tiempo es el vector que contiene (empezando por 0) los tiempos de espera entre reacciones, los cuales tienen distribucion exponencial.
#ENTRADA: 
# tm es el numero de observaciones que se haran al sistema, por lo que habra tm-1 reacciones en la historia
# m es el numero constante de virus que son liberados en el proceso de Lisis
# Coord_S y Coord_S_v son matrices de dim |S|x2 y |S_v|x2 que contienen las coordenadas de S y S_v, respectivamente.
# Si un virus esta a una distancia menor que r de una celula entonces la celula puede ser infectada por ese virus. 
# Las nuevas celulas y virus nacen a una distancia menor que r_p de la celula infectada

#SALIDA: ##########  REVISAR  ###########
#tiempo_hist es el vector que contiene (empezando por 0) los tiempos totales a los que ocurren las reacciones. 
# La diferencia entre cada par de entradas consecutivas de tiempo_hist tiene distribucion exponencial.

# La historia del proceso en dos matrices Hist_Sist_Cel que por renglones describe la historia de cada celula (que siempre empezara siendo sana) y Hist_Sist_V que por renglones describe la historia de cada virus (que siempre empezara vivo).
# Hist_Sist_Cel tiene dim numero de celulas totales en el proceso por 7 columnas: numero que identifica a la celula, su coord en x, su coord en y,S,E,R_I,R (S se usara para nacimiento de celulas
# en donde se indicaran los tiempos (que son entradas de tiempo_hist) en los cuales la celula se volvio parte de cada una de estas poblaciones, resp. Se usara un -1 si la celula nunca llega a una poblacion (E,R_I o R) y 0 si estaba en el sistema desde el principio.)
# Hist_Sist_V tiene dim numero de virus totales en el proceso por 5 columnas: numero que identifica al virus, su coord en x, su coord en y,S_v,R_v (S_v se usara para nacimiento de virus
# en donde se indicaran los tiempos (que son entradas de tiempo_hist) en los cuales el virus se volvio parte de cada una de estas poblaciones, resp. Se usara un -1 si el virus nunca llega a R_v y 0 si estaba en el sistema desde el principio.)

Gillespie_Redes_Lisis<-function(tm,m,Coord_S,Coord_S_v,r,r_p,alp,mu_c,lamb,del,mu_v){
  tiempo<-rep(0,tm)
  tiempo_hist<-rep(0,tm)
  P<-matrix(0,6,tm) # Tamaño de las poblaciones S, E, R_I, R (celulas), S_v, R_v (virus)
 
  Hist_Sist_Cel<-cbind(c(1:length(Coord_S[,1])),Coord_S[,1],Coord_S[,2],rep(0,length(Coord_S[,1])),
                       rep(0,length(Coord_S[,1])),rep(0,length(Coord_S[,1])),rep(0,length(Coord_S[,1]))) #Decimos que comenzamos en tiempo 0?
  Hist_Sist_V<-cbind(c(1:length(Coord_S_v[,1])),Coord_S_v[,1],Coord_S_v[,2],rep(0,length(Coord_S_v[,1])),
                                                                           rep(0,length(Coord_S_v[,1]))) #Decimos que comenzamos en tiempo 0?
  
   Reac<-c(1:6)# Las diferentes reacciones:Nacimiento S
              #                            Muerte natural de S
              #                            Muerte natural de E
              #                            Exposiciones de celulas
              #                            Muerte de celula infectada y liberacion de nuevos virus
              #                            Muerte de S_v
  Coord_E<-Coord_R_I<-matrix(0,1,2) ######## Es necesario guardar las coord de R_I?? Si sí agregar en el codigo el procedimiento#######
  A<-Mat_Ady(Coord_S,Coord_S_v,r) #Mat_Ady(Coord_S,Coord_S_v,r)
  P[,1]<-c(dim(Coord_S)[1],0,0,0,dim(Coord_S_v)[1],0)
  
  k<-c(0,3)
  k[1]<-P[1,1]
  k[2]<-P[2,1]
  
  v<-P[5,1]
  
  for(i in 2:tm){
    WA<-which(A==1,arr.ind = TRUE)
    n_e<-sum(A)
    
    #Revisar que estos casos vayan de acuerdo a Lisis porque estan de acuerdo a Brote:
    if(k[2]==0 && n_e==0){
      Propensidades<-c(rep(alp,k[1]),rep(mu_c,k[1]),rep(mu_v,v))
      Tipos_Reacciones<-c(rep(1,k[1]),rep(2,k[1]),rep(6,v))
    }
    else if(k[2]==0){
      Propensidades<-c(rep(alp,k[1]),rep(mu_c,k[1]),rep(lamb,n_e),rep(mu_v,v))
      Tipos_Reacciones<-c(rep(1,k[1]),rep(2,k[1]),rep(4,n_e),rep(6,v))
    }
    else if(n_e==0){
      Propensidades<-c(rep(alp,k[1]),rep(mu_c,k[1]),rep(mu_c,k[2]),rep(del,k[2]),rep(mu_v,v))
      Tipos_Reacciones<-c(rep(1,k[1]),rep(2,k[1]),rep(3,k[2]),rep(5,k[2]),rep(6,v))
    }
    
    Propensidades<-c(rep(alp,k[1]),rep(mu_c,k[1]),rep(mu_c,k[2]),rep(lamb,n_e),rep(del,k[2]),rep(mu_v,v))
    Tipos_Reacciones<-c(rep(1,k[1]),rep(2,k[1]),rep(3,k[2]),rep(4,n_e),rep(5,k[2]),rep(6,v))
    tiempo[i]<-rexp(1,sum(Propensidades)) #Tiempo exponencial de parametro alp_0=sum(Propensidades)
    tiempo_hist[i]<-tiempo_hist[i-1]+tiempo[i]
    sig<-Sig_Reac(Propensidades)
    tipo_sig<-Tipos_Reacciones[sig]
    if(tipo_sig==1){# Nacimiento de S
      nueva_cel<-c(runif(1,Coord_S[sig,1]-r_p,Coord_S[sig,1]+r_p),runif(1,Coord_S[sig,2]-r_p,Coord_S[sig,2]+r_p))
      Coord_S<-rbind(Coord_S,nueva_cel)
      P[,i]<-P[,i-1]+c(1,0,0,0,0,0)
      Hist_Sist_Cel<-rbind(Hist_Sist_Cel,c(length(Hist_Sist_Cel[,1]+1),nueva_cel[1],nueva_cel[2],tiempo_hist[i],0,0,0))
      
    }else if(tipo_sig==2){# Muerte de S
      Coord_S<-Coord_S[-(sig-k[1]),]
      A<-A[-(sig-k[1]),]   #Aqui marca error
      P[,i]<-P[,i-1]+c(-1,0,0,0,0,0)
      Hist_Sist_Cel[sig-k[1],]<-Hist_Sist_Cel[sig-k[1],]+c(0,0,0,0,0,0,tiempo_hist[i])
    }else if(tipo_sig==3){# Muerte de E
      Coord_S<-Coord_S[-(sig-2*k[1]),]
      P[,i]<-P[,i-1]+c(0,-1,0,0,0,0)
      Hist_Sist_Cel[sig-2*k[1],]<-Hist_Sist_Cel[sig-2*k[1],]+c(0,0,0,0,0,0,tiempo_hist[i])
      
    }else if(tipo_sig==4){# Exposicion
      Coord_E<-rbind(Coord_E,Coord_S[WA[sig,1],])
      Coord_S<-Coord_S[-WA[sig,1],]
      Coord_S_v<-Coord_S_v[-WA[sig,2],]
      A<-A[-WA[sig,1],-WA[sig,2]]
      P[,i]<-P[,i-1]+c(-1,1,0,0,-1,1)
      Hist_Sist_Cel[sig-(2*k[1]+k[2]),]<-Hist_Sist_Cel[sig-(2*k[1]+k[2]),]+c(0,0,0,0,tiempo_hist[i],0,0)
      
      }else if(tipo_sig==5){
      coord_cel_daluz<-Coord_S[-(sig-(2*k[1]+k[2]+n_e)),]
      Coord_S_v<-rbind(Coord_S_v,cbind(runif(m,coord_cel_daluz[1]-r_p,coord_cel_daluz[1]+r_p),
                                       runif(m,coord_cel_daluz[2]-r_p,coord_cel_daluz[2]+r_p)))
      Coord_E<-Coord_E[-(sig-(2*k[1]+k[2]+n_e)),]
      P[,i]<-P[,i-1]+c(0,-1,1,0,m,0) 
      Hist_Sist_Cel[sig-(2*k[1]+k[2]+n_e),]+c(0,0,0,0,0,tiempo_hist[i],0)
      Hist_Sist_V<-rbind(Hist_Sist_V,c(length(Hist_Sist_V[,1])+1,Coord_S_v[length(Coord_S_v)][1],
                                                Coord_S_v[length(Coord_S_v)][2],tiempo_hist[i],0))
      }else if(tipo_sig==6){
      Coord_S_v<-Coord_S_v[-(sig-(2*k[1]+2*k[2]+n_e)),]#Aqui marca error.
      A<-A[,-(sig-(2*k[1]+2*k[2]+n_e))]
      P[,i]<-P[,i-1]+c(0,0,0,0,-1,1)
      Hist_Sist_V[sig-(2*k[1]+2*k[2]+n_e),]<-Hist_Sist_V[sig-(2*k[1]+2*k[2]+n_e),]+c(0,0,0,0,tiempo_hist[i])
    }
    
    k[1]<-P[1,i] 
    k[2]<-P[2,i]
    
    v<-P[5,i]
  }
  
  return(tiempo_hist,Hist_Sist_Cel,Hist_Sist_V) #Regresa el estado del sistema a cada tiempo
} 


#CONTEXTO:
# Aqui las poblaciones son: Poblaciones=(S,E,I,R,S_v,R_v) 
# Estamos suponiendo que cuando comenzamos el proceso el sistema tiene un numero de celulas expuestas e infectadas E(0),I(0)=0 
# El vector de reacciones esta dado con el siguiente orden:
# Reac=(Nacimiento S, Muerte natural de S, Muerte natural de E, Muerte natural de I, Exposicion, Infeccion E->I, Liberacion de S_v, Muerte de S_v) 
# tiempo es el vector que contiene (empezando por 0) los tiempos de espera entre reacciones, los cuales tienen distribucion exponencial.

#ENTRADA: 
# tm es el numero de observaciones que se haran al sistema, por lo que habra tm-1 reacciones en la historia 
# Coord_S y Coord_S_v son matrices de dim |S|x2 y |S_v|x2 que contienen las coordenadas de S y S_v, respectivamente.
# Si un virus esta a una distancia menor que r de una celula entonces la celula puede ser infectada por ese virus. 
# Las nuevas celulas y virus nacen a una distancia menor que r_p de la celula infectada

#SALIDA:  ##########  REVISAR  ###########
#tiempo_hist es el vector que contiene (empezando por 0) los tiempos totales a los que ocurren las reacciones. 
# La diferencia entre cada par de entradas consecutivas de tiempo_hist tiene distribucion exponencial.

# La historia del proceso en dos matrices Hist_Sist_Cel que por renglones describe la historia de cada celula (que siempre empezara siendo sana) y Hist_Sist_V que por renglones describe la historia de cada virus (que siempre empezara vivo).
# Hist_Sist_Cel tiene dim numero de celulas totales en el proceso por 7 columnas: numero que identifica a la celula, su coord en x, su coord en y,S,E,I,R (S se usara para nacimiento de celulas
# en donde se indicaran los tiempos (que son entradas de tiempo_hist) en los cuales la celula se volvio parte de cada una de estas poblaciones, resp. Se usara un -1 si la celula nunca llega a una poblacion (E,I o R) y 0 si estaba en el sistema desde el principio.)
# Hist_Sist_V tiene dim numero de virus totales en el proceso por 5 columnas: numero que identifica al virus, su coord en x, su coord en y,S_v,R_v (S_v se usara para nacimiento de virus
# en donde se indicaran los tiempos (que son entradas de tiempo_hist) en los cuales el virus se volvio parte de cada una de estas poblaciones, resp. Se usara un -1 si el virus nunca llega a R_v y 0 si estaba en el sistema desde el principio.)


Gillespie_Redes_Brote<-function(tm,Coord_S,Coord_S_v,r,r_p,alp,mu_c,lamb,del,gam,mu_v){
  tiempo<-rep(0,tm)
  tiempo_hist<-rep(0,tm)
  P<-matrix(0,6,tm) # Guardara el tamaño de las poblaciones S, E, R_I, R (celulas), S_v, R_v (virus)
  
  Hist_Sist_Cel<-cbind(c(1:length(Coord_S[,1])),Coord_S[,1],Coord_S[,2],rep(0,length(Coord_S[,1])),
                       rep(0,length(Coord_S[,1])),rep(0,length(Coord_S[,1])),rep(0,length(Coord_S[,1]))) #Decimos que comenzamos en tiempo 0?
  Hist_Sist_V<-cbind(c(1:length(Coord_S_v[,1])),Coord_S_v[,1],Coord_S_v[,2],rep(0,length(Coord_S_v[,1])),
                                                                           rep(0,length(Coord_S_v[,1]))) #Decimos que comenzamos en tiempo 0?
  
  Reac<-c(1:8)# Las diferentes reacciones: Nacimiento S
  #                            Muerte natural de S
  #                            Muerte natural de E
  #                            Muerte natural de I
  #                            Exposiciones de celulas
  #                            Infeccion de celulas   
  #                            Liberacion de virus
  #                            Muerte de S_v
  Coord_E<-Coord_I<-matrix(0,1,2) #Inicializa la matriz de coordenadas de E. 
  A<-Mat_Ady(Coord_S,Coord_S_v,r)
  P[,1]<-c(length(Coord_S[,1]),0,0,0,length(Coord_S_v[,1]),0) #Inicializa los tamaños de las poblaciones dim(Coord_S_v)[1]
  
  k<-c(0,3)
  k[1]<-P[1,1]
  k[2]<-P[2,1]
  k[3]<-P[3,1]
  v<-P[5,1]
  
  
  for(i in 2:tm){
    WA<-which(A==1,arr.ind = TRUE)#Regresa un matriz de dos reglones donde cada uno son las coord donde A[i,j]=1.
    n_e<-sum(A)
    if(k[2]==0 && k[3]==0 && n_e==0){
      Propensidades<-c(rep(alp,k[1]),rep(mu_c,k[1]),rep(mu_v,v))
      Tipos_Reacciones<-c(rep(1,k[1]),rep(2,k[1]),rep(8,v))
    }
    else if(k[2]==0 && k[3]==0){
      Propensidades<-c(rep(alp,k[1]),rep(mu_c,k[1]),rep(lamb,n_e),rep(mu_v,v))
     Tipos_Reacciones<-c(rep(1,k[1]),rep(2,k[1]),rep(5,n_e),rep(8,v))
    }
    else if(k[3]==0 && n_e==0){
      Propensidades<-c(rep(alp,k[1]),rep(mu_c,k[1]),rep(mu_c,k[2]),rep(del,k[2]),rep(mu_v,v))
      Tipos_Reacciones<-c(rep(1,k[1]),rep(2,k[1]),rep(3,k[2]),rep(6,k[2]),rep(8,v))
    }
    else if(k[2]==0 && n_e==0){
      Propensidades<-c(rep(alp,k[1]),rep(mu_c,k[1]),rep(mu_c,k[3]),rep(gam,k[3]),rep(mu_v,v))
      Tipos_Reacciones<-c(rep(1,k[1]),rep(2,k[1]),rep(4,k[3]),rep(7,k[3]),rep(8,v))
    }
    else if(k[2]==0){
      Propensidades<-c(rep(alp,k[1]),rep(mu_c,k[1]),rep(mu_c,k[3]),rep(lamb,n_e),rep(gam,k[3]),rep(mu_v,v))
      Tipos_Reacciones<-c(rep(1,k[1]),rep(2,k[1]),rep(4,k[3]),rep(5,n_e),rep(7,k[3]),rep(8,v))
    }
    else if(k[3]==0){
      Propensidades<-c(rep(alp,k[1]),rep(mu_c,k[1]),rep(mu_c,k[2]),rep(lamb,n_e),rep(del,k[2]),rep(mu_v,v))
      Tipos_Reacciones<-c(rep(1,k[1]),rep(2,k[1]),rep(3,k[2]),rep(5,n_e),rep(6,k[2]),rep(8,v))
    }
    else if(n_e==0){
      Propensidades<-c(rep(alp,k[1]),rep(mu_c,k[1]),rep(mu_c,k[2]),rep(mu_c,k[3]),rep(del,k[2]),rep(gam,k[3]),rep(mu_v,v))
      Tipos_Reacciones<-c(rep(1,k[1]),rep(2,k[1]),rep(3,k[2]),rep(4,k[3]),rep(6,k[2]),rep(7,k[3]),rep(8,v))
    }
    else if(k[2]>0 && k[3]>0 && n_e>0){
    Propensidades<-c(rep(alp,k[1]),rep(mu_c,k[1]),rep(mu_c,k[2]),rep(mu_c,k[3]),rep(lamb,n_e),rep(del,k[2]),
                                                                                  rep(gam,k[3]),rep(mu_v,v))
    Tipos_Reacciones<-c(rep(1,k[1]),rep(2,k[1]),rep(3,k[2]),rep(4,k[3]),rep(5,n_e),rep(6,k[2]),rep(7,k[3]),rep(8,v))
    }
    
    tiempo[i]<-rexp(1,sum(Propensidades)) #Tiempo exponencial de parametro alp_0=sum(Propensidades)
    tiempo_hist[i]<-tiempo_hist[i-1]+tiempo[i]
    sig<-Sig_Reac(Propensidades) #Regresa el indice j de la proxima reaccion
    tipo_sig<-Tipos_Reacciones[sig] #Indica el tipo de reaccion de Propensidades[j]
    if(tipo_sig==1){ #Nacimiento de S
      nueva_cel<-c(runif(1,Coord_S[sig,1]-r_p,Coord_S[sig,1]+r_p),runif(1,Coord_S[sig,2]-r_p,Coord_S[sig,2]+r_p))
      Coord_S<-rbind(Coord_S,nueva_cel)
      P[,i]<-P[,i-1]+c(1,0,0,0,0,0) 
      Hist_Sist_Cel<-rbind(Hist_Sist_Cel,c(length(Hist_Sist_Cel[,1])+1,nueva_cel[1],nueva_cel[2],tiempo_hist[i],0,0,0))

    }else if(tipo_sig==2){ #Muerte de S
      Coord_S<-Coord_S[-(sig-k[1]),]
      A[-(sig-k[1]),]
      P[,i]<-P[,i-1]+c(-1,0,0,0,0,0)
      Hist_Sist_Cel[sig-k[1],]<-Hist_Sist_Cel[sig-k[1],]+c(0,0,0,0,0,0,tiempo_hist[i])
    }else if(tipo_sig==3){ #Muerte de E
      Coord_S<-Coord_S[-(sig-2*k[1]),]
      P[,i]<-P[,i-1]+c(0,-1,0,0,0,0)
      Hist_Sist_Cel[sig-2*k[1],]<-Hist_Sist_Cel[sig-2*k[1],]+c(0,0,0,0,0,0,tiempo_hist[i])
      
      }else if(tipo_sig==4){ #Muerte de I
      Coord_S<-Coord_S[-(sig-(2*k[1]+k[2])),]
      P[,i]<-P[,i-1]+c(0,0,-1,0,0,0)
      Hist_Sist_Cel[sig-(2*k[1]+k[2]),]<-Hist_Sist_Cel[sig-(2*k[1]+k[2]),]+c(0,0,0,0,0,0,tiempo_hist[i])
    }
    
      else if(tipo_sig==5){#Exposicion S+S_v->E
      Coord_E<-rbind(Coord_E,Coord_S[WA[sig,1],])#Aqui
      Coord_S<-Coord_S[-WA[sig,1],]
      Coord_S_v<-Coord_S_v[-WA[sig,2],]
      A[-WA[sig,1],-WA[sig,2]]
      P[,i]<-P[,i-1]+c(-1,1,0,0,-1,1)
      Hist_Sist_Cel[sig-(2*k[1]+k[2]+k[3]),]+c(0,0,0,0,tiempo_hist[i],0,0)
      }else if(tipo_sig==6){ #Infeccion E->I 
      
      Coord_E<-Coord_E[-(sig-(2*k[1]+k[2]+k[3]+n_e)),]
      Coord_I<-rbind(Coord_I,Coord_E[(sig-(2*k[1]+k[2]+k[3]+n_e)),])#Matriz coordenadas de I 
      P[,i]<-P[,i-1]+c(0,-1,1,0,0,0) 
      Hist_Sist_Cel[sig-(2*k[1]+k[2]+k[3]+n_e),]+c(0,0,0,0,0,tiempo_hist[i],0)
      
      }else if(tipo_sig==7){ #Liberacion de S_v 
        coord_cel_daluz<-Coord_S[-(sig-(2*k[1]+2*k[2]+k[3]+n_e)),] 
        Coord_S_v<-rbind(Coord_S_v,cbind(runif(1,coord_cel_daluz[1]-r_p,coord_cel_daluz[1]+r_p),
                                         runif(1,coord_cel_daluz[2]-r_p,coord_cel_daluz[2]+r_p)))
        P[,i]<-P[,i-1]+c(0,0,0,0,1,0) 
        
        Hist_Sist_V<-rbind(Hist_Sist_V,c(length(Hist_Sist_V)+1,Coord_S_v[length(Coord_S_v),1],
                                         Coord_S_v[length(Coord_S_v),2],tiempo_hist[i],0))
      }else if(tipo_sig==8){ #Muerte de S_v
      Coord_S_v<-Coord_S_v[-(sig-(2*k[1]+2*k[2]+2*k[3]+n_e)),] #Aqui marca error.
      A[,-(sig-(2*k[1]+2*k[2]+2*k[3]+n_e))] # Es posible que el error este aqui? porque A puede tener menos columnas.
      P[,i]<-P[,i-1]+c(0,0,0,0,-1,1)
      Hist_Sist_V[sig-(2*k[1]+2*k[2]+2*k[3]+n_e),]<-Hist_Sist_V[sig-(2*k[1]+2*k[2]+2*k[3]+n_e),]+c(0,0,0,0,tiempo_hist[i])
    }
    k[1]<-P[1,i]
    k[2]<-P[2,i]
    k[3]<-P[3,i]
    
    v<-P[5,i]
    
  }
  
  #return(tiempo_hist,Hist_Sist_Cel,Hist_Sist_V) #Regresa el estado del sistema a cada tiempo  
  return(list(tiempo_hist=tiempo_hist,Hist_Sist_Cel=Hist_Sist_Cel,Hist_Sist_V=Hist_Sist_V)) #Regresa el estado del sistema a cada tiempo  
  
  } 


# Revisar la sintaxis del ultimo renglon de una matriz
# Revisar si es necesario usar length o si se puede usar c[i,k] para crear el indice que sigue en Hist_Sist_Cel y Hist_Sist_V
# Falta guardar las coordenadas de R_I durante todo el proceso para Lisis, 
# porque eso se agregara a la matriz Sist_Hist_Cel

# Una vez implementado lo anterior empezar a hacer pruebas para asegurarme que todas las funciones hacen lo que se quiere paso a paso.



##### PARA LAS VISUALIZACIONES DE LISIS#####


d<-70 #La medida del lado de E

Coord_S<-Coloca(50,d)
Coord_S_v<-Coloca(6,d)

tm<-120
m<-100
r<-8 #Alcance de virus
r_p<-3 #Radio donde nacen los virus a partir de la celula infectada
alp<-0.3
mu_c<-0.055
lamb<-0.002
del<-0.24
gam<-0.5
mu_v<-3
  
y<-Gillespie_Redes_Brote(tm,Coord_S,Coord_S_v,r,r_p,alp,mu_c,lamb,del,gam,mu_v)

z<-Gillespie_Redes_Lisis(tm,m,Coord_S,Coord_S_v,r,r_p,alp,mu_c,lamb,del,mu_v)

#w<-z$tiempo_hist

