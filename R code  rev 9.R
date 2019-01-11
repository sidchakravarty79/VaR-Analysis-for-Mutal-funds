# load packages
library(xts)
library(timeSeries)
library(tseries)
library(TTR)
library(quantmod)
library(PerformanceAnalytics)
library(extRemes)
library(fExtremes)
library(evd)
library(ismev)
library(evir)
library(POT)
library(ismev)

#VaR analysis
ret.date<-read.csv("OPGIX.csv",header=T,sep=",") #YYYY-mm-dd

data<-read.zoo("OPGIX.csv",header=T,sep=",",stringsAsFactors = T) #YYYY-mm-dd
fund_returns <- na.trim(CalculateReturns(data$Adj.Close),method="log")

chartSeries(data$Adj.Close, theme=chartTheme("white"), type="candlesticks",
            bar.type = "ohlc",up.col = "black",name="OPGIX",ylab="Price")

chartSeries(fund_returns,theme="white")
charts.PerformanceSummary(fund_returns,colorset=dark6equal,main ="Combined Performance summary") 
chart.Histogram(fund_returns,breaks=40,methods = c("add.density","add.risk"),
                ylab = "",main ="OPGIX fund return distribution")
chart.Boxplot(fund_returns)

kurtosis(fund_returns)
skewness(fund_returns)
stdev(fund_returns)
mean(fund_returns)
Return.annualized(fund_returns,geometric = T)

mean(fund_returns)
#Monthly returns
data.month<-(read.zoo("OPGIX.month.csv",header=T,sep=",",stringsAsFactors = T))
month.returns.data<- na.trim(CalculateReturns(data.month$Adj.Close,method = "log"))
data.SP<-read.zoo("SP.csv",header=T,sep=",",stringsAsFactors = T)
SP.returns.data <- na.trim(CalculateReturns(data.SP$Adj.Close,method = "log"))

# Standard risk metrics
SharpeRatio(month.returns.data,Rf=0.0123,FUN = "StdDev")*12
CAPM.alpha(month.returns.data,SP.returns.data,Rf=0.0123)*12*100
CAPM.beta(month.returns.data,SP.returns.data*1.3,Rf=0.0123)


#Tail distribution
chart.Drawdown(fund_returns)
chart.Histogram(fund_returns,breaks=100,
                methods = c("add.normal"),ylab = "",xlim = c(min(fund_returns), 
                 VaR(fund_returns,alpha = 0.01,tail = "lower")),
                ylim = c(0, 4),main = "OPGIX tail dist of 1-Day @ 99% VaR Portfolio Returns")
# Adjust the ylim values if the histrograms doesn't show up

chart.Histogram(fund_returns,breaks=100,
                methods = c("add.normal"),ylab = "",xlim = c(min(fund_returns), 
                                                             VaR(fund_returns,alpha = 0.05,tail = "lower")),
                ylim = c(0, 10),main = "OPGIX tail dist 1-Day @ 95% VaR Portfolio Returns")

#Var Estimates
chart.VaRSensitivity(fund_returns,methods = c("GaussianVaR","ModifiedVaR","HistoricalVaR"))
chart.VaRSensitivity(fund_returns,methods = c("GaussianES", "ModifiedES", "HistoricalES"))

#Convert data from Zoo to Dataframe
fund.time<-as.Date(ret.date[-1,1])
fund.ret<-as.data.frame(diff(log(data$Adj.Close)))
fund.return<-(cbind(fund.time,-fund.ret))
colnames(fund.return)<-c("time","obs")

mean<-mean(fund.return[,2])
sd<-stdev(fund.return[,2])

#Gaussin VaR
x <- seq(0.9,0.999,length=1000)
var.norm<-0
yES.norm<-0
yVaR.norm<-0
for (i in 1:length(x) )
{
  var.norm[i]<- qnorm(x[i],mean=mean,sd=sd)
  yVaR.norm[i] <- (-mean + sd * qnorm(x[i])) * 10000
  yES.norm[i] <- ES(fund_returns,x[i],method="gaussian")* 10000
}

plot.new()
plot(x, -yVaR.norm, type="l",
     xlab=expression(alpha), ylab="$ Amount",main="Gaussian VaR OPGIX")
lines(x, yES.norm, lty=2, col=2)
legend("bottomleft", legend=c("VaR","ES"),col=1:2, lty=1:2)

#Historical VaR
var.hist<-0
yVaR.hist<-0
yES.hist<-0
for (i in 1:length(x) )
{
  var.hist[i]<- quantile(fund_returns,probs=1-x[i])
  yVaR.hist[i] <- -quantile(fund_returns,probs=1-x[i]) * 10000
  yES.hist[i] <- -ES(fund_returns,x[i],method="historical")* 10000
}
plot.new()
plot(x,-var.norm,type="l", xlab=expression(alpha), 
     ylab="%",main="Gaussian and Historical VaR OPGIX",ylim=c(-0.03,-0.01))
lines(x,var.hist,lty=2, col=2)
legend("bottomleft", legend=c("Gaussian","Historical"),col=1:2, lty=1:2)

plot.new()
plot(x,- yVaR.hist, type="l",
     xlab=expression(alpha), ylab="$ Amount",main="Historical VaR OPGIX")
lines(x, -yES.hist, lty=2, col=2)
legend("bottomleft", legend=c("VaR","ES"),col=1:2, lty=1:2)

################################################End of Normal VaR####
#GPD Data exploration

qplot(-fund_returns,cex=0.8)
title(main = "OPGIX QQplot for Extreme Value Analysis",outer = F,font.main=3)

qplot(-fund_returns,threshold=0.02,cex=0.8)
title(main = "OPGIX QQplot for Extreme Value Analysis (threshold =0.02)",outer = F,font.main=3)

#GEV distrition

fit0 <- fevd(obs,fund.return, type = "GEV", units = "Return%")

# it will give warnings, Ignore it !
plot.new()
par(mfcol=c(1,1))
plot(fit0,type = c("qq"),main="QQ plot for GEV OPGIX")
plot(fit0,type = c("density"),main="Density OPGIX")


#Block approach with 42 blocks
blocks=42
var.block.42<-gev(-fund.return[,2],blocks)

plot.gev(var.block.42,main="QQ residual plot for 42 blocks: OPGIX")
var.boxvalue.42<-0
for (i in 1:length(x) )
{
  b2<-(var.block.42$par.ests[2]/var.block.42$par.ests[1])
  b3<-(-blocks*log(x[i]))
  b3a<-b3^(-var.block.42$par.ests[1]) 
  var.boxvalue.42[i]<-var.block.42$par.ests[3]-b2*(1-b3a)
}

###Block approach with 21 blocks
blocks=21
var.block.21<-gev(-fund.return[,2],blocks)
plot.gev(var.block.21,main="QQ residual plot for 21 blocks: OPGIX")
var.boxvalue.21<-0
for (i in 1:length(x) )
{
  b2<-(var.block.21$par.ests[2]/var.block.21$par.ests[1])
  b3<-(-blocks*log(x[i]))
  b3a<-b3^(-var.block.21$par.ests[1]) 
  var.boxvalue.21[i]<-var.block.21$par.ests[3]-b2*(1-b3a)
}

plot.new()
plot(x,-var.norm,type="l", xlab=expression(alpha), 
     ylab="%",main="Gaussian,Historical, and Block VaR OPGIX",ylim=c(-0.03,-0.01),lty=1,col=1)
lines(x,var.hist,lty=2, col=2)
lines(x,-var.boxvalue.21,lty=3, col=3)
lines(x,-var.boxvalue.42,lty=4, col=4)
legend("bottomleft", legend=c("Gaussian","Historical","Block 21","Block 42"),col=1:4, lty=1:4)


#GPD 
fund.time<-as.Date(ret.date[-1,1])
fund.ret<-as.data.frame(-diff(log(data$Adj.Close)))
fund.return<-(cbind(fund.time,-fund.ret))
colnames(fund.return)<-c("time","obs")

fit1 <- fevd(obs,fund.return, threshold = 0.01, type = "GP", units = "Return%")
plot(fit1)
#POT
pot.01<-gpd(-fund_returns,threshold=0.01) # Check your threshold values
plot(pot.01,main="POT with 0.01 threhold")
rm.01<-riskmeasures(pot.01,x)

pot.02<-gpd(-fund_returns,threshold=0.02) # Check your threshold values
plot(pot.02,main="POT with 0.02 threhold")
rm.02<-riskmeasures(pot.02,x)

pot.03<-gpd(-fund_returns,threshold=0.04) # Check your threshold values
plot(pot.03,main="POT with 0.04 threhold")
rm.03<-riskmeasures(pot.03,x)

pot.04<-gpd(-fund_returns,threshold=0.008) # Check your threshold values
plot(pot.04,main="POT with 0.008 threhold")
rm.04<-riskmeasures(pot.04,x)

plot.new()
plot(x,-rm.01[,2],type="l", xlab=expression(alpha), 
     ylab="%",main="POT VaR for different threhold: OPGIX",ylim=c(-0.03,-0.01),lty=1,col=1)
lines(x,-rm.02[,2],lty=100,col=4)
lines(x,-rm.03[,2],lty=100,col=10)
lines(x,-rm.04[,2],lty=100,col=11)
legend("bottomleft", legend=c("Th 0.01","Th 0.02", "Th 0.04","Th 0.008"),
       col=c(1,4,10,11), lty=c(1,100,100,100))

plot.new()
plot(x,-var.norm,type="l", xlab=expression(alpha), 
     ylab="%",main="Gaussian,Historical,Block and POT VaR OPGIX",ylim=c(-0.03,-0.01),lty=1,col=1)
lines(x,var.hist,lty=1, col=2)
lines(x,-var.boxvalue.21,lty=2, col=3)
lines(x,-var.boxvalue.42,lty=2, col=4)
lines(x,-rm.01[,2],lty=100,col=6)
lines(x,-rm.02[,2],lty=100,col=81)
lines(x,-rm.03[,2],lty=100,col=10)
lines(x,-rm.04[,2],lty=100,col=11)
legend("bottomleft", legend=c("Gaussian","Historical","Block 21","Block 42",
                              "Th 0.01","Th 0.02", "Th 0.04","Th 0.008"),
       col=c(1:5,6,81,10,11), lty=c(1,2,2,100,100,100,100))

return.level(pot3)
## Expected Shortfall
plot.new()
plot(x,-rm.01[,2]*10000,type="l", xlab=expression(alpha), 
     ylab="$",main="Expected short fall Gaussian,Historical,and POT VaR OPGIX",
     lty=1,col=1)
lines(x,-rm.02[,2]*10000,lty=1,col=2)
lines(x,-rm.03[,2]*10000,lty=2,col=3)
lines(x,-rm.04[,2]*10000,lty=2,col=4)
lines(x,yES.norm,lty=1, col=5)
lines(x,-yES.hist, lty=2, col=6)
legend("bottomleft", legend=c("Th 0.01","Th 0.02", "Th 0.04","Th 0.008",
                              "Gaussian","Historical"),
       col=c(1:6), lty=c(1,1,2,2,5,6))


