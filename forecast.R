library(knitr)
library(readxl)
library(MASS)
library(forecast)
library(timeSeries)
library(timeDate)
library(sarima)

setwd("/Users/sunxiaotan/Desktop/University of Toronto/Sta457/STA457 A1")
GFluTrends<-read_excel("case_study_1_fluwatch.xlsx",sheet="Google Flu Trends", skip = 1) 
fluWatch<-read_excel("case_study_1_fluwatch.xlsx", 
sheet="FluWatch-BC", skip = 2)  
tim<-timeSequence(from = "2003-09-28", to = "2015-08-09", by = 
"week") 
tim1<-timeSequence(from = "2003-09-07", to = "2015-08-23", 
by="week") 
GFT<- timeSeries(GFluTrends[,"British Columbia"], charvec = tim)
fluTest<- timeSeries(fluWatch[,"FluTest"], charvec = tim1)

# Split data set into training set and test set
GFT_training<-window(GFT,start = '2003-09-28', end = '2014-08-10')
Flu_training<-window(fluTest,start = '2003-09-28', end = '2014-08-10')
GFT_test<-window(GFT,start = '2014-08-17', end = '2015-08-09')
Flu_test<-window(fluTest,start = '2014-08-17', end = '2015-08-09')

par(mfrow=c(2,1),cex=0.5)
plot(GFT_training,type='b',pch=18,main="Training GFT")
grid()
plot(Flu_training,type='b',pch=18,main="Training Flutest")
grid()

par(mfrow=c(2,1),cex=0.5)
plot(GFT_test,type='b',pch=18,main="Test set GFT")
grid()
plot(Flu_test,type='b',pch=18,main="Test set Flutest")
grid()

# Testing prewhiten functions
PreWhiten.ar<- function(x , ar = NULL){
        if(is.null(ar)) print(" autoregressive coefficients are empty!")
        pwData = numeric(0)
        pwData = filter(x, c(1, -ar),method=c("convo"),sides=1) 
        pwData[!is.na(pwData)]
}

PreWhiten.arma<- function(x , ar = NULL, ma = 0){
        if(is.null(ar) && is.null(ma)) print("both ar and ma coefficients are empty!")
        pwData = numeric(0)
        m = as(modelCoef(new("ArmaModel", ar = ar, ma = ma)), "list")
        eps = numeric(length(x))
        pwData = xarmaFilter(m, x =x, eps = eps, whiten = TRUE) 
        pwData[!is.na(pwData)]
}

LBTest<- function(res, nPQ = 0, m = 24, ifPlot = FALSE){
        stopifnot(nPQ >= 0, m >= 1, m > nPQ)
        n <- length(res)
        lags <- 1:m
        df <- (nPQ+1):m 
        ra <- (acf(res, lag.max = m, plot = FALSE)$acf)[-1]
        QQ <- n * (n + 2) * cumsum((ra^2)/(n - (1:m)))[df]
        
        pv <- 1 - pchisq(QQ, df)
        QQ <- round(QQ, 2)
        a <- matrix(c(df, QQ, pv), ncol = 3)
        dimnames(a) <- list(rep("", length(QQ)), c("m", "Qm", "pvalue"))
        if(ifPlot){
                plot(x = a[,1],y = a[,3],
                     ylim = c(0,1), pch = 15, col =4,
                     ylab = "p-value", xlab = "m",
                     main = "Ljung-Box portmanteau test")
                abline(h =0, col =2)
                grid()
        } else {a}
}

mod.ar <- auto.arima(GFT_training, max.p = 52, max.q = 0, stationary = TRUE)
temp = PreWhiten.arma(GFT_training, ar = mod.ar$coef[1:3], ma = 0) 
temp1 = PreWhiten.ar(GFT_training, ar = mod.ar$coef[1:3])
par(mfrow = c(2,1), cex =0.5)
plot(mod.ar)
plot(c(temp1)-temp[-(1:3)], col = 1, 
     ylab = "Difference", ylim = c(-0.01,0.01), pch = 15,
     main = "Difference of prewhitened time series btw 2 filters")
grid()

#All higher than red line so model okay
mod.arma<-auto.arima(GFT_training, max.p = 52, max.q = 52, stationary = TRUE) 
p = mod.arma$arma[1]; q = mod.arma$arma[2]
coef(mod.arma)
plot(mod.arma)
npq = sum(mod.arma$arma[c(1,2)])
LBTest(mod.arma$residuals, nPQ = npq, m = 52, ifPlot = TRUE)

mod = mod.arma;nAR = mod$arma[1]; nMA = mod$arma[2]

if(nMA!=0){
  xf = PreWhiten.arma(GFT_training, ar = mod$coef[1:nAR], 
                      ma = mod$coef[(1:nMA)+nAR])[-(1:nAR)]
  yf = PreWhiten.arma(Flu_training, ar = mod$coef[1:nAR], 
                      ma=mod$coef[(1:nMA)+nAR])[-(1:nAR)]  
}else{
  xf = PreWhiten.arma(GFT_training, ar = mod$coef[1:nAR], 
                      ma = 0)[-(1:nAR)]
  yf = PreWhiten.arma(Flu_training, ar = mod$coef[1:nAR], 
                      ma=0)[-(1:nAR)] 
}

par(cex=0.75,bg="gray95")
ccf(c(xf), c(yf), lwd=1, ylab="Cross-correlation functions",
    main="CCF of prewhitened GFT and flu test")
abline(v=0, col="gold", lwd=2, lty="dashed")
text(-1, 0.2, "-1", col=2)
text(-2, 0.2, "-2", col=2)

y<-fluTest
x<-GFT
dat<- cbind(y,x, lag(x))[-c(1:4),]
colnames(dat)<-c("fluTest", "GFT", "GFT1")
tim2<-timeSequence(from = '2003-10-05', to = '2015-08-23', by = 'week' )
data<- timeSeries(dat, charvec = tim2)

data.training = window(data, start = "2003-10-05", end = "2014-08-10")
data.test = window(data, start = "2014-08-17", end = "2015-08-09")

mod.tfn = auto.arima(data.training[,1], xreg = data.training[,-1], stationary = TRUE)
coef(mod.tfn)

m = 26
lags = 1:m
df <- (1+2+1):m
n = length(mod.tfn$res)
rccf = ccf(mod$residuals,mod.tfn$residuals, plot = FALSE, lag.max = m)$acf[-(1:m)]

Qm = n* (n + 2) * cumsum((rccf^2)/(n - (0:m)))[df]
pv <- 1 - pchisq(Qm, df)
a = cbind(df, Qm,pv)

par(mfrow = c(1,2))
LBTest(mod.tfn$res, nPQ = 6, ifPlot = TRUE)
plot(x = a[,1],y = a[,3],
     ylim = c(0,1), pch = 15, col =4,
     ylab = "p-value", xlab = "m",
     main = "Cross-correlation check")
abline(h =0.05, col =2)
grid()

#The fit of TFN model
par(mfrow = c(1,1), cex = 0.75)
ts.plot(mod.tfn$fitted, ylab = "", main ="TFN model")
lines(c(Flu_training), pch = 10, col = "green", type ="p")
grid()

# ARIMA Model
mod.arima<-auto.arima(data.training[,1])
coef(mod.arima)
plot(mod.arima)
#Fit of ARIMA
par(mfrow = c(1,1), cex = 0.75)
ts.plot(mod.arima$fitted, ylab = "", main ="ARIMA model")
lines(c(Flu_training), pch = 10, col = "green", type ="p")
grid()
#LB Test for ARIMA
npq<-sum(mod.arima$arma[c(1,2)])
LBTest(mod.arima$res, nPQ = npq,m = 52, ifPlot = TRUE)

# NN Model
mod.nn = forecast::nnetar(Flu_training);mod.nn
mod.nnx = forecast::nnetar(data.training[,1], xreg = 
data.training[,-1]);mod.nnx
#Fit of NN Model
par(mfrow = c(1,1), cex = 0.75)
ts.plot(mod.nn$fitted, ylab = "", main ="Neural network")
lines(c(mod.nnx$fitted), col = "red", lwd = 2, type ="p", pch = 15)
lines(c(Flu_training), pch = 10, col = "green", type ="p")
grid()

# 50 steps ahead forcast
forecast(mod.arima,50)
forecast(mod.nn,50)

Table <- matrix(c(accuracy(mod.arima), accuracy(mod.tfn), accuracy(mod.nn), accuracy(mod.nnx)), nrow = 4, byrow=T)
rownames(Table) <- c("ARIMA", "TFN", "NN", "NNX")
colnames(Table) <- c("ME","RMSE","MAE","MPE","MAPE","MASE","ACF1")
Table

set.seed(100)
#arima
arima_forcast1<-mod.arima %>% forecast(h=1)
arima_forcast4<-mod.arima %>% forecast(h=4)
arima_forcast8<-mod.arima %>% forecast(h=8)
arima_forcast50<-mod.arima %>% forecast(h=50)
#tfn
tfn_forcast1<-mod.tfn %>% forecast(h=1, xreg = data.test[,-1])
tfn_forcast4<-mod.tfn %>% forecast(h=4, xreg = data.test[,-1])
tfn_forcast8<-mod.tfn %>% forecast(h=8, xreg = data.test[,-1])
tfn_forcast50<-mod.tfn %>% forecast(h=50, xreg = data.test[,-1])
#nn
nn_forcast1<-mod.nn %>% forecast(h=1)
nn_forcast4<-mod.nn %>% forecast(h=4)
nn_forcast8<-mod.nn %>% forecast(h=8)
nn_forcast50<-mod.nn %>% forecast(h=50)
#nnx
nnx_forecast1<-mod.nnx %>% forecast(h=1, xreg=data.test[,-1])
nnx_forecast4<-mod.nnx %>% forecast(h=4, xreg=data.test[,-1])
nnx_forecast8<-mod.nnx %>% forecast(h=8, xreg=data.test[,-1])
nnx_forecast50<-mod.nnx %>% forecast(h=50, xreg=data.test[,-1])

# lag=1
set.seed(1)
Table1 <- matrix(c(accuracy(arima_forcast1, x=data.test[,1][1:1]), accuracy(tfn_forcast1, x=data.test[,1][1:1]), accuracy(nn_forcast1, x=data.test[,1][1:1]), accuracy(nnx_forecast1, x=data.test[,1][1:1])), nrow = 4, byrow=T)
rownames(Table1) <- c("ARIMA", "TFN", "NN", "NNX")
Table1<-Table1[,-c(1,3,5,7,9,11,13)]
colnames(Table1) <- c("ME","RMSE","MAE","MPE","MAPE","MASE","ACF1")
Table1

# lag=4
set.seed(4)
Table4 <- matrix(c(accuracy(arima_forcast4, x=data.test[,1][1:4]), accuracy(tfn_forcast4, x=data.test[,1][1:4]), accuracy(nn_forcast4, x=data.test[,1][1:4]), accuracy(nnx_forecast4, x=data.test[,1][1:4])), nrow = 4, byrow=T)
rownames(Table4) <- c("ARIMA", "TFN", "NN", "NNX")
Table4<-Table4[,-c(1,3,5,7,9,11,13)]
colnames(Table4) <- c("ME","RMSE","MAE","MPE","MAPE","MASE","ACF1")
Table4

# lag=8
set.seed(8)
Table8 <- matrix(c(accuracy(arima_forcast8, x=data.test[,1][1:8]), accuracy(tfn_forcast8, x=data.test[,1][1:8]), accuracy(nn_forcast8, x=data.test[,1][1:8]), accuracy(nnx_forecast8, x=data.test[,1][1:8])), nrow = 4, byrow=T)
rownames(Table8) <- c("ARIMA", "TFN", "NN", "NNX")
Table8<-Table8[,-c(1,3,5,7,9,11,13)]
colnames(Table8) <- c("ME","RMSE","MAE","MPE","MAPE","MASE","ACF1")
Table8

# lag=50
set.seed(50)
Table50 <- matrix(c(accuracy(arima_forcast50, x=data.test[,1][1:50]), accuracy(tfn_forcast50, x=data.test[,1][1:50]), accuracy(nn_forcast50, x=data.test[,1][1:50]), accuracy(nnx_forecast50, x=data.test[,1][1:50])), nrow = 4, byrow=T)
rownames(Table50) <- c("ARIMA", "TFN", "NN", "NNX")
Table50<-Table50[,-c(1,3,5,7,9,11,13)]
colnames(Table50) <- c("ME","RMSE","MAE","MPE","MAPE","MASE","ACF1")
Table50
