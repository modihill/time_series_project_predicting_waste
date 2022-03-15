# Clearing the environment
rm(list=ls())

# libraries
library(readr)
library(dplyr)
library(tidyr)
library(TSA)
library(tseries)
library(forecast)
library(fUnitRoots)
library(lmtest)
library(fGarch)
library(CombMSC)
library(stats)
library(FitAR)

# data preparation
waste <- read_csv("Waste_collected_per_month.csv")

waste <- waste %>% separate(date, into = c("month", "day","year"), sep = "/")  %>% 
  separate(year, into = c("year", "time"), sep = " ") %>% 
  arrange(year, month, day) %>% select(month, year, residential)

View(head(waste))
View(tail(waste))


# converting dataframe into time series
wasteTS = ts(waste$residential, start = c(2009,4), end = c(2020,3), frequency = 12)
wasteTS

#----------------------------------------------------------------------

# Functions for Time Series
plot_ts <- function(ts, transformation)
{
  win.graph(width = 45,
            height = 30,
            pointsize = 15)
  plot(ts,
       ylab = "Waste",
       main = c(paste0(toString(transformation),
                       " plot of Waste")),
       type="o")
  points(y=ts,x=time(ts), pch=as.vector(season(ts)))
}


# Function for ACF and PACF
autoCorrelation <- function(ts, time_series)
{
  win.graph(width = 20, 
            height = 15,
            pointsize = 15)
  
  par(mfrow=c(2,1))
  acf(ts,
      lag.max = 48,
      main = c(paste0("ACF plot of ",toString(time_series))))
  
  pacf(ts, 
       lag.max = 48,
       main = c(paste0("PACF plot of ",toString(time_series))))
  par(mfrow=c(1,1))
}


# Functions for its Residual Analysis
residual_model <- function(model, modelname)
{
  res.model = rstandard(model)
  win.graph(width = 20,
            height = 15,
            pointsize=15)
  
  par(mfrow=c(3,2))
  
  plot(y = res.model, 
       x = as.vector(time(wasteTS)),
       main = c(paste0("Time Series plot of Standardised residuals of 
                       fitted ", toString(modelname))),
       xlab = 'Time', 
       ylab = 'Standardized Residuals',
       type = 'o')
  
  hist(res.model,
       main = c(paste0("Histogram of Standardised residuals of 
                       fitted ", toString(modelname))),
       xlab='Standardized Residuals') 
  
  qqnorm(y=res.model, 
         main = c(paste0("QQ plot of standardised residuals of 
                         fitted ", toString(modelname))))
  
  qqline(y=res.model, 
         col = 2, 
         lwd = 1, 
         lty = 2)
  
  print(shapiro.test(res.model))
  
  acf(res.model,
      lag.max = 48,
      main = c(paste0("ACF plot of standardised residuals of 
                      fitted ", toString(modelname))))
  
  pacf(res.model,
       lag.max = 48,
       main = c(paste0("PACF plot of standardised residuals of 
                       fitted ", toString(modelname))))    
  
  
  k=0
  LBQPlot(res.model,lag.max = 30, StartLag = k + 1, k = 0, SquaredQ = FALSE)
  par(mfrow=c(1,1))
}


# Function for SARIMA model
SARIMA <- function(p,d,q, modelname)
{
  ml = arima(wasteTSBC,order=c(p,d,q),
             seasonal = list(order = c(2,1,1), period = 12),
             method = "ML")
  
  print(coeftest(ml))
  residual_model(ml,toString(modelname))
  
  Model <- modelname
  
  AIC <- AIC(ml)
  aic <- cbind.data.frame(Model,AIC)
  AIC_SCORE <- rbind(AIC_SCORE, aic)
  AIC_SCORE <<- AIC_SCORE
  
  BIC <- AIC(ml,k=(log(length(wasteTS))))
  bic <- cbind.data.frame(Model,BIC)
  BIC_SCORE <- rbind(BIC_SCORE, bic)
  BIC_SCORE <<- BIC_SCORE
}

#----------------------------------------------------------------------
# Descriptive Analysis
#----------------------------------------------------------------------

# Original Time Series

par(mar=c(1,1,1,1))
plot_ts(wasteTS, "Time Series")


# ACF and PACF of original Time Series
autoCorrelation(wasteTS, "Time Series")


win.graph(width = 20, 
          height = 15, 
          pointsize=15)
# Scatter Plot of Consecutive Years
plot(y=wasteTS,
     x=zlag(wasteTS),
     ylab='Residential Waste (Tonnes)', 
     xlab='Previous Month Residential Waste (Tonnes)', 
     main = "Scatter plot of Residential Waste (Tonnes) in consequtive months")


# Correlation of Consecutive Years
y = wasteTS             
x = zlag(wasteTS)       
index = 2:length(x)    
cor(y[index],x[index]) %>% round(3)

# -----------------------------------------------------------------------
# Box-cox Transformation
# -----------------------------------------------------------------------

win.graph(width = 20, 
          height = 15, 
          pointsize=15)
BC <- BoxCox.ar(wasteTS)

# box-cox confidence Interval
BC$ci

# maximum likelihood
lambda <- BC$lambda[which(max(BC$loglike) == BC$loglike)]
lambda

# box-cox transformation
wasteTSBC = ((wasteTS^lambda)-1)/lambda

plot_ts(wasteTSBC, "Box-Cox Transformation Time Series")

autoCorrelation(wasteTSBC, "Box-Cox Transformation Time Series")

# -----------------------------------------------------------------------
# Residual Approach SARIMA Models
# -----------------------------------------------------------------------

m1.temp = arima(wasteTSBC,order=c(0,0,0),
                seasonal = list(order = c(0,1,0), period = 12))
res.m1 = residuals(m1.temp)  
plot_ts(res.m1,"Residuals SARIMA(0,0,0)x(0,1,0)")
autoCorrelation(res.m1,"Residuals SARIMA(0,0,0)x(0,1,0)")


m2.temp = arima(wasteTSBC,order=c(0,0,0),
                seasonal = list(order = c(2,1,1), period = 12))
res.m2 = residuals(m2.temp)
plot_ts(res.m2,"Residuals SARIMA(0,0,0)x(2,1,1)")
autoCorrelation(res.m2,"Residuals SARIMA(0,0,0)x(2,1,1)")

# P = 2, D = 1, Q = 1

# p = (1,2)
# q = (1,2,3)

# SARIMA(1,0,1)x(2,1,1), 
# SARIMA(1,0,2)x(2,1,1), 
# SARIMA(1,0,3)x(2,1,1),
# SARIMA(2,0,1)x(2,1,1),
# SARIMA(2,0,2)x(2,1,1), 
# SARIMA(2,0,3)x(2,1,1),

m3.temp = arima(wasteTSBC,order=c(2,0,2),
                seasonal = list(order = c(2,1,1), period = 12))
res.m3 = residuals(m3.temp)  
plot_ts(res.m3,"Residuals SARIMA(2,0,2)x(2,1,1)")
autoCorrelation(res.m3,"Residuals SARIMA(2,0,2)x(2,1,1)")

# -----------------------------------------------------------------------
# EACF
# -----------------------------------------------------------------------
eacf(res.m2)
# SARIMA(0,0,0)x(2,1,1),
# SARIMA(0,0,1)x(2,1,1),

# -----------------------------------------------------------------------
# BIC
# -----------------------------------------------------------------------
win.graph(width = 10, 
          height = 10,
          pointsize = 15)
res = armasubsets(y=res.m2,
                  nar=8,
                  nma=8,
                  y.name='test',
                  ar.method='ols')
plot(res)
# SARIMA(1,0,3)x(2,1,1),
# SARIMA(3,0,3)x(2,1,1),

# -----------------------------------------------------------------------
# set of possible SARIMA models;
# -----------------------------------------------------------------------

# SARIMA(0,0,0)x(2,1,1),
# SARIMA(0,0,1)x(2,1,1),
# SARIMA(1,0,1)x(2,1,1), 
# SARIMA(1,0,2)x(2,1,1), 
# SARIMA(1,0,3)x(2,1,1),
# SARIMA(2,0,1)x(2,1,1),
# SARIMA(2,0,2)x(2,1,1), 
# SARIMA(2,0,3)x(2,1,1),
# SARIMA(3,0,3)x(2,1,1)

# -----------------------------------------------------------------------
# creating null data frame for AIC and BIC
# -----------------------------------------------------------------------

AIC_SCORE <- data.frame(Model = c(),
                        AIC = c())

BIC_SCORE <- data.frame(Model = c(),
                        BIC = c())

# -----------------------------------------------------------------------
# Model Fitting
# -----------------------------------------------------------------------

# SARIMA(0,0,0)x(2,1,1),
SARIMA(0,0,0,"SARIMA(0,0,0) X (2,1,1)")

# SARIMA(0,0,1)x(2,1,1)
SARIMA(0,0,1,"SARIMA(0,0,1) X (2,1,1)")

# SARIMA(1,0,1)x(2,1,1)
SARIMA(1,0,1,"SARIMA(1,0,1) X (2,1,1)")

# SARIMA(1,0,2)x(2,1,1)
SARIMA(1,0,2,"SARIMA(1,0,2) X (2,1,1)")

# SARIMA(1,0,3)x(2,1,1)
SARIMA(1,0,3,"SARIMA(1,0,3) X (2,1,1)")

# SARIMA(2,0,1)x(2,1,1)
SARIMA(2,0,1,"SARIMA(2,0,1) X (2,1,1)")

# SARIMA(2,0,2)x(2,1,1)
SARIMA(2,0,2,"SARIMA(2,0,2) X (2,1,1)")

# SARIMA(2,0,3)x(2,1,1)
SARIMA(2,0,3,"SARIMA(2,0,3) X (2,1,1)")

# SARIMA(3,0,3)x(2,1,1)
SARIMA(3,0,3,"SARIMA(3,0,3) X (2,1,1)")


# -----------------------------------------------------------------------
# Finding Best Model
# -----------------------------------------------------------------------
AIC_SCORE %>% unique() %>%  arrange(AIC)
BIC_SCORE %>% unique() %>%  arrange(BIC)

# -----------------------------------------------------------------------
# SARIMA(3,0,3)x(2,1,1) is best model, now, let's try overfit
# -----------------------------------------------------------------------

# SARIMA(4,0,3)x(2,1,1)
SARIMA(4,0,3,"SARIMA(4,0,3) X (2,1,1)")

# SARIMA(3,0,4)x(2,1,1)
SARIMA(3,0,4,"SARIMA(3,0,4) X (2,1,1)")

# -----------------------------------------------------------------------
# SARIMA(3,0,4)x(2,1,1) is best model, now, let's try overfit
# -----------------------------------------------------------------------

# SARIMA(3,0,5)x(2,1,1)
SARIMA(3,0,5,"SARIMA(3,0,5) X (2,1,1)")

# SARIMA(4,0,4)x(2,1,1)
SARIMA(4,0,4,"SARIMA(4,0,4) X (2,1,1)")


# -----------------------------------------------------------------------
# SARIMA(3,0,5)x(2,1,1) is best model, now, let's try overfit
# -----------------------------------------------------------------------

# SARIMA(3,0,6)x(2,1,1)
SARIMA(3,0,6,"SARIMA(3,0,6) X (2,1,1)")

# SARIMA(4,0,5)x(2,1,1)
SARIMA(4,0,5,"SARIMA(4,0,5) X (2,1,1)")


# -----------------------------------------------------------------------
# SARIMA(4,0,5)x(2,1,1) is best model, now, let's try overfit
# -----------------------------------------------------------------------

# SARIMA(4,0,6)x(2,1,1)
SARIMA(4,0,6,"SARIMA(4,0,6) X (2,1,1)")

# SARIMA(5,0,5)x(2,1,1)
SARIMA(5,0,5,"SARIMA(5,0,5) X (2,1,1)")

# -----------------------------------------------------------------------
# Forecasting
# -----------------------------------------------------------------------
forecasting = Arima(wasteTS,order=c(4,0,5),
                    seasonal=list(order=c(2,1,1), period=12), 
                    lambda = -0.1,
                    method = "ML")
prediction = forecast(forecasting, h = 10)

win.graph(width = 10, 
          height = 10,
          pointsize = 15)
plot(prediction)

prediction
