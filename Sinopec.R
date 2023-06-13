library("forecast")
library("lmtest")
library("tseries")
library("vars")
library("urca")
library("TSA")
library("Matrix")
library("matlib")
library("portes")
library("xts")
library("PerformanceAnalytics")
library("quantmod")
library("patchwork")
library("ggplot2")
library("rugarch")
library("tsDyn")

library ("broom")
library("dplyr")    
library("tidyr") 
library("astsa")

library("corrr")
library("ggcorrplot")
library("FactoMineR")
library("devtools")
library("factoextra")
library("dynamac")
library("Metrics")
library("stargazer")

library("xtable")

data.initial <-read.csv("datamock.csv")

data = data.initial[order(nrow(data.initial):1),]

#LLadj <- dataLaShSi$LLadj
#LShadj <- dataLaShSi$adjLSh
#LSiadj <- dataLaShSi$adjLSi

#data <- cbind(data, LLadj,LShadj,LSiadj)
plot.ts(data$LL)
plot.ts(data$LSh)
plot.ts(data$LSi)
summary(data)


stargazer(data, 
          type = "html",
          out = "description_China.htm",
          out.header = TRUE)

write.csv(data, "data_for_python.csv")
#----------------------------------stationarity--------------------------------


Pacf(data$LSi)
ur.df(data$LSi, type="drift", lags = 1, 
      selectlags = "Fixed")

d1L<-diff(data$LSi, differences=1)

plot.ts(d1L)

Pacf(d1L)
ura <- ur.df(d1L, type="drift", lags = 14, 
             selectlags = "Fixed") 
#-8.1852
summary (ura)

d1OG<-diff(data$OGlog, differences=1)
d1E<-diff(data$EnergyL, differences=1)
d1I<-diff(data$IndustL, differences=1)
d1Mine<-diff(data$MineL, differences=1)
d1T<-diff(data$TEchL, differences=1)
d1U<-diff(data$UtilL, differences=1)
d1Cong<-diff(data$ConglL, differences=1)
d110Y<-diff(data$L10Y, differences=1)
d1COil<-diff(data$COilL, differences=1)
d1S300<-diff(data$S300L, differences=1)
d1MSCI<-diff(data$MSCIL, differences=1)
d1S100<-diff(data$S100L, differences=1)
d1Comp<-diff(data$CompL, differences=1)
d1CNY <- diff(data$CNYL, differences = 1)
d1IPI <- diff(data$dataIC.IPI, differences = 1)
d1CPI <- diff(data$dataIC.CPI, differences = 1)
d1Gold<- diff(data$GoldL, differences = 1)
d1DJIA <- diff(data$DJIAL, differences = 1)

Pacf(d1OG)
ur.df(d1OG, type="drift", lags = 23, 
      selectlags = "Fixed")
Pacf(d1E)
ur.df(d1E, type="drift", lags = 14, 
      selectlags = "Fixed") 
Pacf(d1I)
ur.df(d1I, type="drift", lags = 24, 
      selectlags = "Fixed") 
Pacf(d1Mine)
ur.df(d1Mine, type="drift", lags = 24, 
      selectlags = "Fixed") 
Pacf(d1T)
ur.df(d1T, type="drift", lags = 24, 
      selectlags = "Fixed") 
Pacf(d1U)
ur.df(d1U, type="drift", lags = 27, 
      selectlags = "Fixed") 
Pacf(d1Cong)
ur.df(d1Cong, type="drift", lags = 24, 
      selectlags = "Fixed") 
Pacf(d110Y)
ur.df(d110Y, type="drift", lags = 8, 
      selectlags = "Fixed") 
Pacf(d1COil)
ur.df(d1COil, type="drift", lags = 18, 
      selectlags = "Fixed") 
Pacf(d1S300)
ur.df(d1S300, type="drift", lags = 27, 
      selectlags = "Fixed") 
Pacf(d1MSCI)
ur.df(d1MSCI, type="drift", lags = 24, 
      selectlags = "Fixed") 
Pacf(d1S100)
ur.df(d1S100, type="drift", lags = 26, 
      selectlags = "Fixed") 
Pacf(d1Comp)
ur.df(d1Comp, type="drift", lags = 24, 
      selectlags = "Fixed") 
Pacf(d1CNY)
ur.df(d1CNY, type="drift", lags = 4, 
      selectlags = "Fixed") 

Pacf(d1IPI)
ur.df(d1IPI, type="drift", lags = 4, 
      selectlags = "Fixed") 
Pacf(d1CPI)
ur.df(d1CPI, type="drift", lags = 23, 
      selectlags = "Fixed") 
Pacf(d1Gold)
ur.df(d1Gold, type="drift", lags = 12, 
      selectlags = "Fixed") 
Pacf(d1DJIA)
ur.df(d1DJIA, type="drift", lags = 24, 
      selectlags = "Fixed") 




data2 <- data[-c(0:1), ]
df = data.frame(data2$data7.Date, d1L, d1OG, d1E, d1I, d1Mine, d1T, d1U, 
                d1Cong, d110Y, d1COil, d1S300, d1MSCI, d1S100, d1Comp, d1CNY,
                d1CPI, d1Gold, d1DJIA)

#--------------------------------------Arima&Garch------------------------------
eac1 <- eacf(d1L)

ARMAmodel <- Arima(data$LSi, c(2,1,1), include.constant =TRUE, method = c("CSS-ML"))  
coeftest(ARMAmodel)
summary(ARMAmodel)
SinLBT <- Box.test(residuals(ARMAmodel), lag = 7, type = c("Ljung-Box"), fitdf = 3)
forecast_ARMA<-forecast(ARMAmodel, h=70, level = c(0.90,0.95))
plot(forecast_ARMA)
forecast_ARMA
write.csv(forecast_ARMA, "ARMA_Sinopec_forecast.csv")

stargazer(ARMAmodel$coef, 
          type = "html",
          out = "ARMA_Sinopec.html",
          out.header = TRUE)

stargazer(eac1,
          type = "html",
          out = "eacSin.html",
          out.header = TRUE)

#ARCH эффект
Acf(residuals(ARMAmodel))
Pacf(residuals(ARMAmodel))
Box.test(residuals(ARMAmodel)^2, lag = 7, type = c("Ljung-Box"), fitdf = 1)
shapiro.test(residuals(ARMAmodel))



specse<- ugarchspec(mean.model = list(armaOrder = c(2,1), include.mean = TRUE),
                    variance.model = list(model = 'sGARCH', garchOrder = c(3,1)),
                    distribution.model = 'std')

modele <- ugarchfit(data = data$LSi, spec = specse)
modele@fit$LLH


Acf(residuals(modele))
Acf(residuals(modele)^2)
Acf(residuals(modele, standardize="TRUE"))
Acf(residuals(modele, standardize="TRUE")^2)

write.csv(residuals(modele, standardize="TRUE"), "garch_modele1.csv")
strese <-read.csv("garch_modele1.csv", sep=',', dec='.')

Box.test(strese$V1, lag = 7, type = c("Ljung-Box"), fitdf = 3) #p-value = 0.1182
Box.test(strese$V1^2, lag = 7, type = c("Ljung-Box"), fitdf = 3) #p-value = 0.7837

#----------------------------------graphs for comparison------------------------
data$data7.Date <- anytime(data$data7.Date)
data$LL <- as.numeric(data$LL)
xtsD <- as.xts(data$LL, order.by = data$data7.Date, frequency=NULL, RECLASS=FALSE)
xtsD
chartSeries(xtsD)
returnD <- CalculateReturns(xtsD)
returnD <- returnD[-1]
returnD
hist(returnD)
plot(returnD)
chart.Histogram(returnD,
                methods = c('add.density', 'add.normal'),
                colorset = c('blue', 'green','red'))

strese$V1
hist(strese$V1)
plot(strese$V1)
chart.Histogram(strese$V1,
                methods = c('add.density', 'add.normal'),
                colorset = c('blue', 'green','red'))

#--------------------------------------Granger test-----------------------------


#H0 - Time series X does not cause time series Y to Granger-cause itself.
#Ha - Time series X  cause time series Y to Granger-cause itself
# if *** we refuse H0, meaning x defines y
#in ccf negative for x leads to y (ccf(x,y,lag.max, type))
#grangertest (x,y)


#2
ccf(d1L, d1OG, lag.max = 23, type = c("correlation"),
    plot = TRUE)

grangertest(d1L, d1OG, order = 5) # 0.205

grangertest(d1OG, d1L, order = 8) #  6.765e-09 OG granger causes L


#10
ccf(d1L, d1COil, lag.max = 27, type = c("correlation"),
    plot = TRUE)

grangertest(d1L, d1COil, order = 10) # 0.0001647

grangertest(d1COil, d1L, order = 19) # 0.2774

#14
ccf(d1L, d1CNY, lag.max = 14, type = c("correlation"),
    plot = TRUE)

grangertest(d1L, d1CNY, order = 1) # 0.01803*

grangertest(d1Comp, d1L, order = 3) # 0.6242

#15
ccf(d1L, d1CPI, lag.max = 23, type = c("correlation"),
    plot = TRUE)

grangertest(d1L, d1CPI, order = 6) # 0.294-

grangertest(d1CPI, d1L, order = 6) # -

#16
ccf(d1L, d1Gold, lag.max = 14, type = c("correlation"),
    plot = TRUE)

grangertest(d1L, d1Gold, order = 8) # 0.4622-

grangertest(d1Gold, d1L, order = 12) # 0.4219-

#17
ccf(d1L, d1DJIA, lag.max = 24, type = c("correlation"),
    plot = TRUE)

grangertest(d1L, d1DJIA, order = 22) # -

grangertest(d1DJIA, d1L, order = 13) # 0.7842-

#18
ccf(d1L, d1IPI, lag.max = 24, type = c("correlation"),
    plot = TRUE)

grangertest(d1L, d1IPI, order = 16) # -

grangertest(d1IPI, d1L, order = 3) # 0.3812-

#19
ccf(d1L, d1CPI, lag.max = 24, type = c("correlation"),
    plot = TRUE)

grangertest(d1L, d1CPI, order = 6) # -

grangertest(d1CPI, d1L, order = ) # -

#20
ccf(d1L, df_ex$ISF, lag.max = 24, type = c("correlation"),
    plot = TRUE)

grangertest(d1L, df_ex$ISF, order = 9) # 0.0002221-

grangertest(df_ex$ISF, d1L, order = 8) # 2.326e-09 -


#--------------------------------Principal component method---------------------
df_ex = data.frame(d1E, d1I, d1Mine, d1T, d1U, d1S300, d1S100, d1Comp, d1DJIA,
                   d1Cong, d110Y, d1MSCI, d1Gold, d1CPI, d1IPI)

#industrial sector factors
dfforpac = data.frame(d1E, d1I, d1Mine, d1T, d1U)
data_normalized <- scale(dfforpac)
corr_matrix <- cor(data_normalized)
ggcorrplot(corr_matrix)
data.pca <- princomp(corr_matrix)
summary(data.pca)
dataPCAISF <- data.pca$loadings[, 1:2]
fviz_eig(data.pca, addlabels = TRUE)
fviz_cos2(data.pca, choice = "var", axes = 1:2)

df_ex <- df_ex %>% mutate (ISF = (0.4075852*d1E - 0.3293501*d1I + 0.5712718*d1Mine - 0.6205376*d1T -
                                    d1U*0.1182556) + 0.1558730*d1E + 0.2055331*d1I - 0.1428964*d1Mine -
                             0.3104718*d1T +  0.9036856*d1U)
stargazer(dataPCAISF, 
          type = "html",
          out = "PCA_ISF.html",
          out.header = TRUE)


#Blue-chips & major companies
dfforpac1 = data.frame(d1S300, d1S100, d1Comp, d1DJIA, d1Cong)
data_normalized1 <- scale(dfforpac1)
corr_matrix1 <- cor(data_normalized1)
ggcorrplot(corr_matrix1)
data.pca1 <- princomp(corr_matrix1)
summary(data.pca1)
dataPCABCC <- data.pca1$loadings[, 1:1]
fviz_eig(data.pca1, addlabels = TRUE)
fviz_cos2(data.pca1, choice = "var", axes = 1:2)

df_ex <- df_ex %>% mutate (BCC = 0.4229833*d1S300 + 0.4640462*d1S100  +
                             0.4629385*d1Comp -0.4936296*d1DJIA + 0.3844010*d1Cong )

stargazer(dataPCABCC, 
          type = "html",
          out = "PCA_BCC.html",
          out.header = TRUE)




#Financial security
dfforpac3 = data.frame(d110Y, d1MSCI, d1Gold)
data_normalized3 <- scale(dfforpac3)
corr_matrix3 <- cor(data_normalized3)
ggcorrplot(corr_matrix3)
data.pca3 <- princomp(corr_matrix3)
summary(data.pca3)
dataPCAFC <- data.pca3$loadings[, 1:2]
fviz_eig(data.pca3, addlabels = TRUE)
fviz_cos2(data.pca3, choice = "var", axes = 1:2)

df_ex <- df_ex %>% mutate (FC = 0.4912615*d110Y + 0.3513476*d1MSCI -0.7970050*d1Gold + 
                             0.65821637*d1Gold - 0.74903184*d1MSCI +0.07551502*d110Y)


stargazer(dataPCAFC, 
          type = "html",
          out = "PCA_FC.html",
          out.header = TRUE)

exog <- data.frame(df_ex$ISF, df_ex$BCC, df_ex$FC)

#---------------------------xts format, dataframes and exogeneous dataframe-------------------
df$data2.data7.Date <- anytime(df$data2.data7.Date)

df$d1L <- as.numeric(df$d1L)
ydL <- as.xts(df$d1L, order.by = df$data2.data7.Date, frequency=NULL, RECLASS=FALSE)
xOG <- as.xts(df$d1OG, order.by = df$data2.data7.Date, frequency=NULL, RECLASS=FALSE)
xE <- as.xts(df$d1E, order.by = df$data2.data7.Date, frequency=NULL, RECLASS=FALSE)
xI <- as.xts(df$d1I, order.by = df$data2.data7.Date, frequency=NULL, RECLASS=FALSE)
xMine <- as.xts(df$d1Mine, order.by = df$data2.data7.Date, frequency=NULL, RECLASS=FALSE)
xT <- as.xts(df$d1T, order.by = df$data2.data7.Date, frequency=NULL, RECLASS=FALSE)
xU <- as.xts(df$d1U, order.by = df$data2.data7.Date, frequency=NULL, RECLASS=FALSE)
xCong <- as.xts(df$d1Cong, order.by = df$data2.data7.Date, frequency=NULL, RECLASS=FALSE)
x10Y <- as.xts(df$d110Y, order.by = df$data2.data7.Date, frequency=NULL, RECLASS=FALSE)
xOil <- as.xts(df$d1COil, order.by = df$data2.data7.Date, frequency=NULL, RECLASS=FALSE)
x300 <- as.xts(df$d1S300, order.by = df$data2.data7.Date, frequency=NULL, RECLASS=FALSE)
xMSCI <- as.xts(df$d1MSCI, order.by = df$data2.data7.Date, frequency=NULL, RECLASS=FALSE)
x100 <- as.xts(df$d1S100, order.by = df$data2.data7.Date, frequency=NULL, RECLASS=FALSE)
xComp <- as.xts(df$d1Comp, order.by = df$data2.data7.Date, frequency=NULL, RECLASS=FALSE)
xCNY <- as.xts(df$d1CNY, order.by = df$data2.data7.Date, frequency=NULL, RECLASS=FALSE)
xGold <- as.xts(df$d1Gold, order.by = df$data2.data7.Date, frequency=NULL, RECLASS=FALSE)
xDJIA <- as.xts(df$d1DJIA, order.by = df$data2.data7.Date, frequency=NULL, RECLASS=FALSE)
xCPI <- as.xts(df$d1CPI, order.by = df$data2.data7.Date, frequency=NULL, RECLASS=FALSE)
xISF<- as.xts(df_ex$ISF, order.by = df$data2.data7.Date, frequency=NULL, RECLASS=FALSE)
xBCC<- as.xts(df_ex$BCC, order.by = df$data2.data7.Date, frequency=NULL, RECLASS=FALSE)
xFC<- as.xts(df_ex$FC, order.by = df$data2.data7.Date, frequency=NULL, RECLASS=FALSE)


df1 = data.frame(ydL, xOG)

df2 = data.frame(ydL, xOG, xFC)
#---------------------------------------VAR model-------------------------------


ARMAmodel <- Arima(data$LSi, c(0,1,4), include.constant =TRUE, method = c("CSS-ML"))  

specse<- ugarchspec(mean.model = list(armaOrder = c(0,4), include.mean = TRUE),
                    variance.model = list(model = 'sGARCH', garchOrder = c(2,2)),
                    distribution.model = 'std')

modele <- ugarchfit(data = d1L, spec = specse)
modele 

#--------------------------------------VAR select-------------------------------------

VARselect(df1, lag.max = 30, type="const")
aa<-VARselect(df1, lag.max = 30, type="const")
aa$selection

var<-VAR(df1, p = 6, type = "const")
summary(var)
coef(var)

Hosking(var, lags=1.5*var$p) # 0.4414947
LiMcLeod(var, lags=1.5*var$p) # 0.4393276

resvar1<-restrict(var, method = c("ser"), thresh = 1)
summary(resvar1)

#-------------------------ARIMA(rss) vs VAR(2v, thresh = 0, ess)----------------
n<-length(coef(ARMAmodel))-1

x<-var$p-n+1
rss<-sum(ARMAmodel$residuals[x:length(residuals(ARMAmodel))]^2)

n1<-resvar1$varresult$ydL$rank-1
ess<-sum(resvar1$varresult$ydL$residuals^2)

Fstat<-((rss-ess)/(n1-n))/(ess/(length(ydL)-2*n1-1))
Fstat #12.13944

#'F крит 1%'
qf(0.99, df1=n1-n, df2=length(ydL)-2*n1-1) # 4.6303
#'p-val'
pf(Fstat, n1-n, length(ydL)-2*n1-1, lower.tail=F) # 6.339501e-06


#-----------------------------VAR(2v, rss) vs VAR(3v, ess)--------------------------------

#модель с 2-мя переменными

VARselect(df1, lag.max = 30, type="const")
aa<-VARselect(df1, lag.max = 30, type="const")
aa$selection

var<-VAR(df1, p = 2, type = "const")
summary(var)
coef(var)

Hosking(var, lags=1.5*var$p) # 0.2619134
LiMcLeod(var, lags=1.5*var$p) # 0.254382

resvar1<-restrict(var, method = c("ser"), thresh = 1.5)
summary(resvar1)

n01<-resvar1$varresult$ydL$rank-1
x2<-30-2+1
rss01<-sum(resvar1$varresult$ydL$residuals[x2:length(resvar1$varresult$ydL$residuals)]^2)

#'модель с 3-мя переменными'

a2 <- VARselect(df2, lag.max = 30, type="const")
a2$selection
var2<-VAR(df2, p = 2, type = "const")
summary(var2)
Hosking(var2, lags=1.124*var2$p) #0.07808312
LiMcLeod(var2, lags=1.124*var2$p) #0.06928767


resvar2<-restrict(var2, method = c("ser"), thresh = 1)
summary(resvar2)
Hosking(resvar2, lags=1.35*var2$p) # 0.9324441
LiMcLeod(resvar2, lags=1.45*var2$p) # 0.06928767


n11<-resvar2$varresult$ydL$rank-1
ess11<-sum(resvar2$varresult$ydL$residuals^2)

#Fstat
Fstat3<-((rss01-ess11)/(n11-n01))/(ess11/(length(ydL)-2*n11-1))
Fstat3 #-18.74093

#'F крит 1%'
qf(0.99, df1=n11-n01, df2 = length(ydL)-2*n01-1) # 3.34115
#'p-val'
pf(Fstat3, n11-n01, length(ydL)-2*n01-1, lower.tail=F)

#--------------------------------ARDL same as VAR with exogeneous variables-------------------------------------------

VARselect(df1, lag.max = 30, type="const", exogen = exog)
aa<-VARselect(df1, lag.max = 30, type="const", exogen = exog)
aa$selection

var<-VAR(df1, p = 6, type = "const", exogen = exog)
summary(var$varresult$ydL)
coef(var)

resvar1<-restrict(var, method = c("ser"), thresh = 1.1)
summary(resvar1)
coef(resvar1)

Hosking(resvar1, lags= 1.5*var$p) # 0.1162797
LiMcLeod(resvar1, lags=1.5*var$p) # 0.1105711



s<-irf(resvar1, impulse = "xOG", response = c("ydL"), cumulative=TRUE, boot =
         TRUE)
plot(s)


modelVARe <- VAR(df1, p = 6, exogen = exog, type = "const", ic = c("AIC", "HQ", "SC", "FPE"))
tmp <- summary(modelVARe)

xtable(tmp$varresult$ydL)
#-------------------------------------test Johansen----------------------------------

dfJ <- data.frame(df1, exog)
sjd.vecm <- ca.jo(dfJ, ecdet = "trend", type="eigen", K= 25, spec="transitory",
                  season=NULL)
summary(sjd.vecm)

sjd.vecm2 <- ca.jo(dfJ, ecdet = "const", type="eigen", K= 25, spec="transitory",
                   season=NULL)
summary(sjd.vecm2)

#-------------------------------------VECM---------------------------------------

varvecm <- VARselect(dfJ, lag.max = 30, type = "const")
varvecm$selection

vecm.eg<-VECM(dfJ, lag = 2, include = c("none"), LRinclude = c("none"))
summary(vecm.eg)

vecm.con<-VECM(dfJ, lag = 2, include = c("const"), estim="ML", LRinclude = c("none"))
summary(vecm.con)

vecm.tr<-VECM(dfJ, lag=2, include = c("trend"), estim="ML", LRinclude = c("none"))
summary(vecm.tr)

vecm.both <- VECM(dfJ, lag = 2, include = c("both"), estim="ML", LRinclude = c("none"))
summary(vecm.both)


#-------------------------------forecast for VAR with exogeneous variables--------------------------
datamock_after <- read.csv("datamock_after.csv")
dm_after <- datamock_after[order(nrow(datamock_after):1),]
dm_after$data7.Date <- anytime(dm_after$data7.Date)
dm_after <- dm_after[1:70, ]
#--------------------------------forecast VAR with ARIMA ex vars----------------
eacf(d1OG)

ARMAmodelOG <- Arima(d1OG, c(1,0,1), include.constant =TRUE, method = c("CSS-ML"))  
Box.test(residuals(ARMAmodelOG), lag = 7, type = c("Ljung-Box"), fitdf = 2) #p-value = 0.8552

forecast_ARMAOG<-forecast(ARMAmodelOG, h=70)

eacf(df_ex$ISF)

ARMAmodelISF <- Arima(df_ex$ISF, c(0,0,4), include.constant =TRUE, method = c("CSS-ML"))  
Box.test(residuals(ARMAmodelISF), lag = 7, type = c("Ljung-Box"), fitdf = 4) #p-value = 0.2402

forecast_ARMAISF<-forecast(ARMAmodelISF, h=70)

eacf(df_ex$BCC)

ARMAmodelBCC <- Arima(df_ex$BCC, c(0,0,0), include.constant =TRUE, method = c("CSS-ML"))  
Box.test(residuals(ARMAmodelBCC), lag = 7, type = c("Ljung-Box"), fitdf = 0) #p-value =  0.2364

forecast_ARMABCC<-forecast(ARMAmodelBCC, h=70)

eacf(df_ex$FC)

ARMAmodelFC <- Arima(df_ex$FC, c(0,0,1), include.constant =TRUE, method = c("CSS-ML"))  
Box.test(residuals(ARMAmodelFC), lag = 7, type = c("Ljung-Box"), fitdf = 1) #p-value = 0.1566

forecast_ARMAFC<-forecast(ARMAmodelFC, h=70)




df_ex.ISF <- forecast_ARMAISF$mean
df_ex.BCC <- forecast_ARMABCC$mean
df_ex.FC <- forecast_ARMAFC$mean

ef <- data.frame(dm_after$data7.Date, df_ex.ISF)
ef$df_ex.ISF <- as.numeric(ef$df_ex.ISF)

df_ex.ISF <- as.xts(ef$df_ex.ISF, order.by = ef$dm_after.data7.Date, frequency=NULL, RECLASS=FALSE)

exogf = data.frame(df_ex.ISF, df_ex.BCC, df_ex.FC)
exogf$df_ex.ISF <- as.numeric(exogf$df_ex.ISF)
exogf$df_ex.BCC <- as.numeric(exogf$df_ex.BCC)
exogf$df_ex.FC <- as.numeric(exogf$df_ex.FC)


forecastVARef <- predict(resvar1, dumvar = exogf, n.ahead = 70, ci = 0.95)
plot.ts(forecastVARef$fcst$ydL)
write.csv(forecastVARef$fcst$ydL, "forecast_VARef_Sinopec.csv")
VARef <- data.frame(forecastVARef$fcst$ydL)

#-----------------------forecast VECM-----------------------------------------

modelVV <- vec2var(sjd.vecm)

modelVV2 <- vec2var(sjd.vecm2)

forecast_trend <- predict(modelVV, n.ahead = 70, ci = 0.95)
plot.ts(forecast_trend$fcst$ydL)
write.csv(forecast_trend$fcst$ydL, "forecast_VECM_trend.csv")
VECM_tr <- data.frame(forecast_trend$fcst$ydL)

forecast_const <- predict(modelVV2, n.ahead = 70, ci = 0.95)
plot.ts(forecast_const$fcst$ydL)
write.csv(forecast_const$fcst$ydL, "forecast_VECM_const.csv")
VECM_const <- data.frame(forecast_const$fcst$ydL)
#----------------------------------------------------------------
dataf <- read.csv("datamock_after.csv")
dataf = dataf[order(nrow(dataf):1),]
dataf$data7.Date <- anytime(dataf$data7.Date)
dataf$LL <- as.numeric(dataf$LL)
dataf <- dataf[1:70, ]
dL = diff(dataf$LL, differences=1)
dataERROR = data.frame(dataf$data7.Date, VECM_tr$fcst, VECM_const$fcst,
                       VARef$fcst)
dataER = data.frame(dataERROR[2:70, ], dL)
write.csv(dataER, "data_Lanzhou_forecast.csv")

#--------------------------RMSE MAE MAPE--------------------------------------


RMSE_VAR_ef <- sqrt(mean((dataER$dL - dataER$VARef.fcst)^2))
RMSE_VECMc <- sqrt(mean((dataER$dL - dataER$VECM_const.fcst)^2))
RMSE_VECMt <- sqrt(mean((dataER$dL - dataER$VECM_tr.fcst)^2))

RMSE <- data_frame(RMSE_VAR_ef, RMSE_VECMc, RMSE_VECMt)
RMSE

MAEef <- mae(actual = dataER$dL, predicted = dataER$VARef.fcst)
MAEvecmc <- mae(actual = dataER$dL, predicted = dataER$VECM_const.fcst)
MAEvecmt <- mae(actual = dataER$dL, predicted = dataER$VECM_tr.fcst)

MAE <- data.frame(MAEef, MAEvecmc, MAEvecmt)
MAE

MAPEef <- MAPE(dataER$VARef.fcst, dataER$dL)
MAPEvecmc <- MAPE(dataER$dL, dataER$VECM_const.fcst)
MAPEvecmt <- MAPE(dataER$dL, dataER$VECM_tr.fcst)

MAPE <- data.frame(MAPEef, MAPEvecmc, MAPEvecmt)
MAPE

#---------------------prophet act-MA(27)--------------------
dataMA215 <- read.csv("MA(27)Sin.csv")

dataMA215$dif <- gsub(",", ".", dataMA215$dif) 
dataMA215$dif <- as.numeric(dataMA215$dif)
dataMA215 <- dataMA215 %>% mutate (difr = dif*100)
dataMA215 <- dataMA215 %>% mutate (difro = signif(difr, digits = 2))

datapr <- data.frame(dataMA215$dif, dataMA215$Sinopec)
tail(datapr)
datapr <- datapr %>% 
  rename("ds" = "dataMA215.Sinopec",
         "y" = "dataMA215.dif")
datapr$ds <- as.Date(datapr$ds, "%d.%m.%Y")


model1 <- prophet(datapr)
future1 <- make_future_dataframe(model1, periods = 90)
tail(future1)

forecast1 <- predict(model1, future1)
tail(forecast1[c('ds','yhat','yhat_lower','yhat_upper')])
write.csv(forecast1, "prophet_d1_Sinopec.csv")

dyplot.prophet(model1, forecast1)
prophet_plot_components(model1, forecast1)















