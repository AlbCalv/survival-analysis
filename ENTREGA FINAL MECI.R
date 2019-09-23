library(survival)
data(cancer)
help(cancer)
#Obtenemos la submuestra de menores de 55 años
datos<-cancer[rowSums(is.na(cancer))==0,]
#Tenemos que ordenar por tiempo para realizar el análisis correctamente
datos<-datos[order(datos$time),]
funcion_supervivencia<-survfit(Surv(datos$time,datos$status)~1)
plot(funcion_supervivencia,xlab = "t",ylab="S(t)",main="Función de Supervivencia")

#Hacemos el cálculo manual para observar que sale igual
tiempo<-datos$time[datos$status==2]
S<-funcion_supervivencia$surv[datos$status==2]
F<-1-S
plot(tiempo,F,type="s",xlab = "t",ylab="F(t)",main="Función de Distribución")#type="s" la hace escalonada

f<-(c(F,0)-c(0,F))[-(length(F)+1)]
plot(tiempo,f,cex=0.7,xlab = "t",ylab="f(t)",main="Función de Densidad")

h<-f/S
plot(tiempo,h,cex=0.7,xlab = "t",ylab="h(t)",main="Función de Riesgo")
round(h,2)


#Modelo de cox
cox<-coxph(formula = Surv(datos$time, datos$status==2) ~., data = datos)
summary(cox)

#Modelo de cox para las variables significativas únicamente
cox2<-coxph(formula = Surv(datos$time, datos$status==2) ~inst+sex+ph.ecog+ph.karno+wt.loss, data = datos)
summary(cox2)

#Test de Riesgos proporcionales
test.cox<-cox.zph(cox)
test.cox
test.cox2<-cox.zph(cox2)
test.cox2

#Gráficos con los residuos Schoenfeld para validar la hipótesis de RP
x11()
par(mfrow=c(3,3))
plot(test.cox)
x11()
par(mfrow=c(3,3))
plot(test.cox2)

#------Análisis de Residuos con el modelo completo-----------

#Residuos Martingala
scatter.smooth(datos$time, residuals(cox,type="martingale"),xlab="Tiempo",ylab="Residuos",main="Residuos Martingala")
#Utiles para valorar la forma funcional de las variables explicativas
sum(residuals(cox,type="martingale"))

#Residuos Deviance
#Utiles para determinar la existencia de Outliers
#Como se acerca a una aproximación Gaussiana, consideraremos outliers valores 
#fuera de [-3,3] incluso de [-2.5,2.5] que en este caso sucede para 
#la primera observacion
resids.dev<-residuals(cox,type="deviance")
lp <- predict(cox, type="lp" )
plot(lp, resids.dev, xlab="Predictor Lineal", ylab="Residuos",main="Residuos Deviance")
abline(h=2.45,col="red")

#Podemos observar diferencias
summary(datos[,c("time","ph.ecog", "ph.karno", "pat.karno", "meal.cal", "wt.loss" )])
cbind(resids.dev,datos[,c("time","ph.ecog", "ph.karno", "pat.karno", "meal.cal", "wt.loss" )] )[ abs(resids.dev) >= 2.4, ]
#Valores de ph.ecog y ph.karno buenísimos y sin embargo muerte prematura

#Residuos Score
#Utiles para validar la hipótesis de RP
resids.score<-residuals(cox,type="score")
#Plots Score junto con las etiquetas para los puntos influyentes
plot(datos[,1],resids.score[,1],xlab ="Inst",ylab="Residuos Score",main="ResScore*Inst")
text(datos[163,1], resids.score[163,1], labels="163", cex= 1,pos = 3,col="red")
plot(datos[,4],resids.score[,2],xlab ="Edad",ylab="Residuos Score",main="ResScore*Edad")
text(datos[162,4], resids.score[162,2], labels="162", cex= 1,pos = 4,col="red")
plot(datos[,5],resids.score[,3],xlab ="Sexo",ylab="Residuos Score",main="ResScore*Sexo")
plot(datos[,6],resids.score[,4],xlab ="Ph.ecog",ylab="Residuos Score")
plot(datos[,7],resids.score[,5],xlab ="Ph.karno",ylab="Residuos Score")
plot(datos[,8],resids.score[,6],xlab ="Pat.karno",ylab="Residuos Score",main="ResScore*Pat.karno")
text(datos[163,8], resids.score[163,6], labels="163", cex= 1,pos = 4,col="red")
plot(datos[,9],resids.score[,7],xlab ="Meal.cal",ylab="Residuos Score",,main="ResScore*Meal.cal")
text(c(datos[163,9],datos[97,9]), c(resids.score[163,7],resids.score[97,7]), labels=c("163","97"), cex= 1,pos = 4,col="red")
plot(datos[,10],resids.score[,8],xlab ="Wt.loss",ylab="Residuos Score",main="ResScore*Wt.loss")
text(c(datos[163,10],datos[68,10]), c(resids.score[163,8],resids.score[68,8]), labels=c("163","68"), cex= 1,pos = c(4,1),col="red")
#Para ver cuales eran los puntos influyentes
which.min(resids.score[,1]) #163
which.max(resids.score[,2]) #162
which.max(resids.score[,6]) #163
which.min(resids.score[,7]) #97 
which.max(resids.score[,7]) #163
which.min(resids.score[,8]) #163
which.max(resids.score[,8]) #68

#Residuos Schoenfeld
#Utiles para validar la hipótesis de RP
resids.scho<-residuals(cox,type="schoenfeld")
