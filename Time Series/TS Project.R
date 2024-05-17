# MAT482 Data Analysis: Times Series Project
# This is an analysis of Mississippi average temperatures from 1950 through 2020.
# A prediction will be made for 2021 temperatures. 
# Then, the prediction will be compared to real values for 2021.

rm(list = ls())

# Get the data [Stage 1]
setwd("/Users/razor/Desktop/MAT 482/TS Project")

# msat = MS Average Temperatures (recorded monthly)
data <- ts(read.csv("MS_avg_temp_1950-2022.csv")[,4], st=1950, fr=12)
msat.ts <- window(data, st=1950, end=2020+(11/12))                      # saving the last period for later
#plot(decompose(msat.ts))


# Labeled time series plot of the data [Stage 2]
plot(msat.ts, main='Avgerage Monthly Temperature of MS', ylab='Temperature (ºF)', xlab='Years')
plot(decompose(msat.ts)$trend, main='Trend from Decomposition', ylab='Temperature (ºF)', xlab='Years')


# Linear regression [Stage 3]
trend <- decompose(msat.ts)$trend
ms.t <- time(msat.ts) # time variable for msat.ts [1950,2021)
linfit <- lm(trend ~ ms.t)
I <- linfit$coef[1]   # intercept
S <- linfit$coef[2]   # slope
linmodel <- S*ms.t + I
summary(linfit)

# Checking some confidence intervals
#confint(linfit, 'ms.t')
#confint(linfit, '(Intercept)') 

# Plot of the data with linear regression
ts.plot(cbind(msat.ts,trend), lty=c(4,1),col=c('gray','black'), main='Linear Fit on the Trend', ylab='Temperature (ºF)', xlab='Years')
lines(linmodel, col='blue')
legend('bottomleft', inset=0.05, c('Data', 'Trend', 'Linear fit'), col=c('gray','black', 'blue'), lty=c(4,1,1), title='Legend')

# Nonlinear regression  
a_s <- 0.02   # a start
b_s <- -18    # b start
c_s <- 6.28   # c start
d_s <- -5     # d start
e_s <- 20     # e start
nlfit <- nls(msat.ts ~ a*ms.t + b*sin(c*ms.t + d) + e, start=list(a=a_s, b=b_s, c=c_s, d=d_s, e=e_s))
A <- coef(nlfit)[1]
B <- coef(nlfit)[2]
C <- coef(nlfit)[3]
D <- coef(nlfit)[4]
E <- coef(nlfit)[5]
fit <- A*ms.t + B*sin(C*ms.t + D) + E
summary(nlfit)

# Checking some confidence intervals
#confint(nlfit, 'a')
#confint(nlfit, 'b')
#confint(nlfit, 'c')
#confint(nlfit, 'd')
#confint(nlfit, 'e')

# Messing with time
future.t <- time(ts(fit, st=2021, end=2022+(1/6), fr=12))  # time variable for the future fit [2021,2022)
past <- A*ms.t + B*sin(C*ms.t + D) + E
future <- A*future.t + B*sin(C*future.t + D) + E

# Combined plot of msat.ts, fit.ts, and future.ts (with legend)
ts.plot(cbind(msat.ts, past, future), col=c('black','blue', 'red'), xlim=c(1950,2022), main='Nonlinear Model', ylab='Temperature (ºF)', xlab='Years')
legend('bottomleft', inset=0.05, c('Data', 'Fit', 'Prediction'), col=c('black','blue', 'red'), lty=1, title='Legend')


# Residuals (nlfit v. data) [Stage 4]
res <- residuals(nlfit)
#plot(past + res.ts - msat.ts)
res.ts <- ts(res, st=1950, fr=12)
plot(res.ts, main='Nonlinear Fit v. Input Data', xlab='Years', ylab='Residuals')
library(forecast)

plot(decompose(res.ts), xlab='Years')


# Autocorrelation function [Stage 5]
acf(res.ts, main='Nonlinear Model')
# acf(decompose(res.ts)$seas, main='Seasonal of NL Model')
# comparing these correlelograms confirms that there is significant seasonal correlation in the residuals


# AR model [Stage 6]                           
res.ar <- ar(res.ts)                              # ar of the residuals time series
ar_fit <- past + res.ar$resid                     # nlfit of the input data + ar residuals
ar_fut_res <- predict(res.ar, n.ahead=15)$pred    # prediction for future residuals
ar_future <- future + ar_fut_res                  # future fit + ar residuals

# Generate a plot of all three parts
ts.plot(cbind(msat.ts,ar_fit,ar_future), col=c('black','blue','red'), main='AR Model', ylab='Temperature (ºF)', xlab='Years')
legend('bottomleft', inset=0.05, c('Data', 'Fit', 'Prediction'), col=c('black','blue', 'red'), lty=1, title='Legend')
plot(res.ar$resid,main='AR Fit v. Input Data', xlab='Years', ylab='Residuals')
summary(res.ar$resid)
summary(res.ts)

# AR residual [Stage 7]
acf(na.omit((res.ar$resid)), main='AR Model')

# The get.best.arima function
get.best.arima <- function(x.ts, maxord = c(1,1,1,1,1,1))
{
  best.aic <- 1e8
  n <- length(x.ts)
  
  for(p in 0:maxord[1])
    for (d in 0:maxord[2])
      for(q in 0:maxord[3])
        for(P in 0:maxord[4])
          for(D in 0:maxord[5])
            for(Q in 0:maxord[6])
            {
              try( 
                {
                  fit <- arima(x.ts, order=c(p,d,q), seas = list(order=c(P,D,Q), frequency(x.ts) ), method = "CSS") 
                  fit.aic <- -2*fit$loglik + (log(n)+1) * length(fit$coef) #suppressing code after the + gives better fits (using loglik as the only criterion that way)
                  if (fit.aic < best.aic)
                  {
                    best.aic <- fit.aic
                    best.fit <- fit
                    best.model <- c(p,d,q,P,D,Q)
                  } #end if
                } # end first argument of try
                , FALSE) #end try
              
              print(c(p,d,q,P,D,Q))
              flush.console()
            } # end for	
  
  
  
  dev.new()
  
  
  over <- paste("Process Fit with ARIMA(", toString(best.model[1]), ",", toString(best.model[2]), ",", toString(best.model[3]), ") Process. \n Coefficients:", toString(round(best.fit$coef, digits=3))) 
  under <- paste("Periodic Coefficients:", toString(best.model[4]), ",", toString(best.model[5]), ",", toString(best.model[6]))
  
  acf(na.omit(best.fit$resid), lag.max=100, main=over, xlab=under)
  
  
  
  list(akaike=best.aic, data=best.fit, ordersbest_arima=best.model)
} # end get.best.arima	

#get.best.arima(res.ts, c(24,1,1,1,1,2))
best_arima <- arima(res.ts, order=c(1,0,1), seas=list(order=c(1,1,2), 12), method = "CSS")
summary(best_arima)

  
# ARIMA model [Stage 8]
arima_res <- best_arima$res
arima_fut_res <- predict(best_arima, n.ahead=15)$pred
arima_fit <- past + arima_res
arima_future <- future + arima_fut_res
#plot(best_arima$data)

# Generate a plot of all three parts
ts.plot(cbind(msat.ts,arima_fit,arima_future), col=c('black','blue','red'), main='ARIMA Model', ylab='Temperature (ºF)', xlab='Years')
legend('bottomleft', inset=0.05, c('Data', 'Fit', 'Prediction'), col=c('black','blue', 'red'), lty=1, title='Legend')
plot(arima_res, main='ARIMA Fit v. Input Data', xlab='Years', ylab='Residuals')
summary(arima_res)
summary(res.ar$resid)


# ARIMA residual [Stage 9]
acf(arima_res, main='ARIMA Model')


# Reload the original data set [Stage 10]
plot(data)


# Extract the last period [Stage 11]
vault <- window(data, st=2021, end=2022+(1/6))
plot(vault)


# Write a function to compute differences [Stage 12]
difference <- function(x, y){
  length.x <- length(x)
  length.y <- length(y)
  
  if(length.x == length.y){
    numerator <- sum(abs(x-y)[1:length.x])^2
    denominator <- sum(abs(y)[1:length.y])^2
    
    return((1-(numerator/denominator))*100)
  } 
  else{
    stop('The length of vectors x and y must be equal')
  }
}


# Compute differences for every prediction [Stage 13]
difference(vault, future)
difference(vault, ar_future)
difference(vault, arima_future)


# Plot showing all the predictions v. the input data
ts.plot(cbind(vault, future, ar_future, arima_future), col=c('black','red','blue','purple'), main='These Three Models v. the Data', xlab='Years', ylab='Temperature (ºF)')
legend('bottom', inset=0.05, c('Data','NL Model', 'AR Model', 'ARIMA Model'), col=c('black','blue', 'red','purple'), lty=1, title='Legend')

# Calculate explained variance for all three models
exp_var <- function(data, estimator){
  length.data <- length(data)
  length.est <- length(estimator)
  data.avg <- mean(data)
  
  if(length.data == length.est){
    numerator <- sum(((data - estimator)^2)[1:length.est])
    #print('numerator:'); print(numerator)
    denominator <- sum(((data - data.avg)^2)[1:length.data])
    #print('denominator:'); print(denominator)
    
    return((1-(numerator/denominator))*100)
  }
  else{
    stop('The length of the two vectors must be equal')
  }
}

exp_var(msat.ts, past)
exp_var(window(msat.ts, st=1951+(11/12)), window(ar_fit, st=1951+(11/12)))
exp_var(msat.ts, arima_fit)
