library("rmgarch")
library("vars")

install.packages("rmgarch")
install.packages("vars")




library(rmgarch)
datos.mexico<-resumen[,c(24,25,26,27)]
dim(datos.mexico)


ndata<-dim(datos.mexico)
nwindow_2<-70


#Estimar un modelo de vector autoregresivo 

vfit   <- varxfit(X=datos.mexico[1:((dim(datos.mexico)[1])-nwindow_2),], p=1, exogen = NULL, robust = FALSE,
               gamma = 0.25, delta = 0.01, nc = 10, ns = 500, postpad = "constant")

uspec  <- ugarchspec(mean.model = list(armaOrder = c(0,0), 
                                     include.mean = FALSE), variance.model = list(garchOrder = c(1,1), model = "sGARCH"),
                   distribution.model = "norm")

spec   <- dccspec(uspec = multispec(replicate(4, uspec)), VAR = TRUE,
               lag = 1, dccOrder = c(1,1), distribution = "mvnorm")#, asymmetric = FALSE

fit1    <- dccfit(spec, data = datos.mexico, fit.control = list(eval.se=TRUE), 
             VAR.fit = vfit,out.sample = nwindow_2) 

tail(datos.mexico)
fitted(fit1)


#NOTAS :
#out.sample tiene efecto directo sobre el tamaño de la muestra que va a vfit
#con la funcion fitted puedes recobrar las predicciones 


dcc.focast     <- list()
i              <- 1
dcc.focast[[i]]<- dccforecast(fit1, n.ahead = 1, n.roll = nwindow_2)

plot(dcc.focast[[1]], which=1) #grafico forecast vs dato real , forecast in sample 



#alpha,beta,omega --> parametros del modelo garch
#dcca1 y dccb1    --> parametros desconocidos por mi 




##################################################################################################
#
# 2da aproximacion
#
##################################################################################################


# 

#install.packages("tseries")
library("tseries")

#install.packages(file.choose(), lib="C://Program Files//R//R-3.6.1//library",repos = NULL) #install.packages("ccgarch")
library("ccgarch")
library(readxl)
library(xts)

setwd("C://Users//JESÚS//Desktop//universidad//2020-1//tesis 2")
datos_series_consolidado    <-read_excel("deamed.xlsx",sheet = "final")
datos_series_consolidado    <-read_excel("BD_price_LN_no_negs_brasil.xlsx",sheet = "ln") #BD dolarizada stocks brasileños

resumen                     <-datos_series_consolidado
datos.mexico                <-resumen[,c(24,25,26,27)]
datos.brasil                <-resumen[,c(2,3,4,5)] # el objeto tibble contabiliza las fechas como una columna del  dataframe mas
dim(datos.mexico)



Ti = nrow(datos.mexico)
Ti
t = 200 #Conjunto de datos para el testing 
n = ncol(datos.mexico)
n
alpha95 = 1.65
alpha99 = 2.33

VaR95 = matrix(, nrow = Ti - t, ncol = 1) #AL observar var95 me di cuenta que Ti-t recien es la ventana, entendido como el espacio donde se va a llevar acabo el movimiento 
dim(VaR95)
VaR99 = matrix(, nrow = Ti - t, ncol = 1)
dim(VaR99)


#simulacion del valor del portafolio (PROVISIONAL)
for(i in 1:(dim(datos.mexico)[1])){
datos.mexico[i,5]=1000+rnorm(1,4,1)
} 

check<-as.vector(as.data.frame(datos.mexico[,5]))
#VER1 value = matrix(check, nrow = Ti) #via este comando no generaba una matriz adecuada 
value = data.matrix(check)

W = matrix(unlist(datos.mexico[1,c(1,2,3,4)]/value[1, ]), nrow = n, ncol = 1)

#ver2: W = matrix(c(0.25,0.25,0.25,0.25), nrow = n, ncol = 1) #para no perder tiempo arbitrariamente fije los pesos

for (j in 1:(Ti - t)) {# longitud de BD original - longitud de ventana =BD de entrenamiento
  
      ts = datos.mexico[j:(j + t -1), ] # se actualiza ts con una muestra del tamaño de la ventana 1:200 2:201 3:202 y asi. de la BASE DE DATOS DE RETORNOS
      
      residual = matrix(, nrow = t, ncol = n) # t tamaño de ventana  y n=numero de columnas de la BD original sin retornos
      
      for (i in 1:n) {#Para cada n columna de la BD ORIGINAL
        residual[, i] = matrix(residuals(arma(ts[, i], order = c(1, 0))))
       #para cada columna de residual( el num de col = num de stock) = residuos de la estimacion de un arma de cada columna de ts osea de una bd con 200 observaciones 
       #En cada iteracion del bucle mas grande se resetea y el segundo bloquees para guardar los residuos de cada i-esimo stock  con longitud = t ( ejm 2000) 
        }
  
      residual = residual[-1, ]           # Elimina la primera fila 
      
      coef = matrix(, nrow = n, ncol = 3) # Matriz de coeficientes con n  filas =  numero de stock , son 3 columnas por el numero de parametros del arma garch 
                                          # Se resetean los coeficientes 
      # initial GARCH model estimation
      for (i in 1:n) { #Aca estaba 10 ahora lo reemplazo por n=numero de stocks
        coef[i, ] = matrix(coef(garch(residual[, i], order = c(1, 1), series = NULL)))
        
      }
      
      # DCC-GARCH model estimation
      a = coef[, 1]       #intercepto
      A = diag(coef[, 2]) #parte glusten
      B = diag(coef[, 3]) #parte arch
      dcc.para = c(0.01, 0.97)
      results = dcc.estimation(inia = a, iniA = A, iniB = B, ini.dcc = dcc.para, dvar = residual, 
                               model = "diagonal")
      h = results$h
      dcc = results$DCC
      v = sqrt(diag(h[t - 1, ]))
      R = matrix(data = dcc[1, ], nrow = n, ncol = n)
      H = v %*% R %*% v
      
      VaR95[j, ] = sqrt(t(W) %*% H %*% W) * alpha95 * value[j]  #aca se almacena defrente los resultos de los value at risk 
      VaR99[j, ] = sqrt(t(W) %*% H %*% W) * alpha99 * value[j]  #se calcula el VaR con lo de 200pts y ya 
      
}

VaR = cbind(VaR95, VaR99)

# backtest profit&loss
valuesd = matrix(value, nrow = Ti)
PL = valuesd[(t + 1):Ti, ] - valuesd[t:(Ti - 1), ]

# install.packages('parallel') install.packages('rugarch')
# library("parallel")
library("rugarch")

ind_test = function(V) {
  J = matrix(ncol = 4, nrow = length(V))
  for (i in 2:length(V)) {
    J[i, 1] = V[i - 1] == 0 & V[i] == 0
    J[i, 2] = V[i - 1] == 0 & V[i] == 1
    J[i, 3] = V[i - 1] == 1 & V[i] == 0
    J[i, 4] = V[i - 1] == 1 & V[i] == 1
  }
  V_00 = sum(J[, 1], na.rm = TRUE)
  V_01 = sum(J[, 2], na.rm = TRUE)
  V_10 = sum(J[, 3], na.rm = TRUE)
  V_11 = sum(J[, 4], na.rm = TRUE)
  p_00 = V_00/(V_00 + V_01)
  p_01 = V_01/(V_00 + V_01)
  p_10 = V_10/(V_10 + V_11)
  p_11 = V_11/(V_10 + V_11)
  hat_p = (V_01 + V_11)/(V_00 + V_01 + V_10 + V_11)
  a = (1 - hat_p)^(V_00 + V_10) * (hat_p)^(V_01 + V_11)
  b = (p_00)^(V_00) * (p_01)^(V_01) * (p_10)^(V_10) * p_11^(V_11)
  #Note: if p_11 and V_11 both are 0, then p_11^(V_11)=1
  return(-2 * log(a/b))
}



V95 = matrix(, nrow = (Ti - t))
for (i in 1:(Ti - t)) {
  if (PL[i] > (-VaR95[i, 1])) {
    V95[i] = 0
  } else {
    V95[i] = 1
  }
}





V99 = matrix(, nrow = (Ti - t))
for (i in 1:(Ti - t)) {
  if (PL[i] > (-VaR99[i, 1])) {
    V99[i] = 0
  } else {
    V99[i] = 1
  }
}

inde95 = ind_test(V95)
inde99 = ind_test(V99)



library("rugarch")

vt95 = VaRTest(alpha = 0.05, as.numeric(PL), as.numeric(-VaR95), conf.level = 0.95)
vt99 = VaRTest(alpha = 0.01, as.numeric(PL), as.numeric(-VaR99), conf.level = 0.95) #VaRTest ojo hay 3 vartest con distintas mayus y minus
vt95
vt99

cc.stat95 = vt95$uc.LRstat + inde95
cc.stat99 = vt99$uc.critical + inde99

if (cc.stat95 > vt95$cc.critical) {
  print("Reject H0")
} else {
  print("Fail to Reject H0")
}

if (cc.stat99 > vt99$cc.critical) {
  print("Reject H0")
} else {
  print("Fail to Reject H0")
}

bt = cbind(PL, -VaR)
colnames(bt) = c("PL", "95%VaR", "99%VaR")
matplot(c(1:(Ti - t)), bt[, 1:3], type = "l", xlab = "time", ylab = "P&L")
legend("topright", colnames(bt)[-1], lwd = 1, col = 2:3)
title("Portfolio P&L and estimated VaR")




### OJO : POR LOS VALORES CALCULADOS PUEDE DEJAR DE FUNCIONAR 
vt95 = VaRTest(alpha = 0.05,c( 0.034107089, -0.434369347 , 0.326791026) , c(-1.66630 , 6.13079 , -5.51437), conf.level = 0.95)
vt95

#intentar con logaritmos 

