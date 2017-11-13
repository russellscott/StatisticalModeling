# Statistical Modelling Final Project Part B
library(MuMIn)
library(MASS)
library(xtable)
library(glmnet)
library(mgcv)
library(SemiPar)
setwd("C:/Users/Russell/Documents/KUL/Year2_Winter/Stat Modeling/final")

full = read.table("car_accidents.txt",header=T)
set.seed(603723) # student number r0603723
rownumbers = sample(1:1286,size=300)
mydata.B = full[rownumbers,]
rownames(mydata.B) = seq(1, nrow(mydata.B))

X = mydata.B[,-8]
pairs(X)
hist(mydata.B$Severity, main = "Distribution of the Outcome, Severity", 
     xlab = "Severity", breaks = )
par(mfrow=c(1,2))
scatter.smooth(mydata.B$Year, mydata.B$Severity, family = "gaussian", 
               xlab = "Year", ylab = "Severity",
               main = "Scatter Plot")
scatter.smooth(mydata.B$Safety, mydata.B$Severity, family = "gaussian",
               xlab = "Safety", ylab = "Severity",
               main = "Scatter Plot")
par(mfrow=c(1,1))
scatter.smooth(mydata.B$Speed, mydata.B$Severity, 
               family = "gaussian", xlab = "Speed")
scatter.smooth(mydata.B$Power, mydata.B$Severity, 
               family = "gaussian", xlab = "Power")
scatter.smooth(mydata.B$Airbag, mydata.B$Severity, 
               family = "gaussian", xlab = "Airbag")
scatter.smooth(mydata.B$Temp, mydata.B$Severity, 
               family = "gaussian", xlab = "Temp")
scatter.smooth(mydata.B$Wind, mydata.B$Severity, 
               family = "gaussian", xlab = "Wind")
# non linear relationships with Severity, expect higher order terms

fit1 = lm(Severity ~ Year + Safety + Speed + Power + Airbag + Temp + Wind, 
          data = mydata.B)
summary(fit1)

# center variables to make sense of intercept
mydata.B.center = mydata.B
mydata.B.center$Year = mydata.B$Year - mean(mydata.B$Year)
mydata.B.center$Safety = mydata.B$Safety - mean(mydata.B$Safety)
mydata.B.center$Speed = mydata.B$Speed - mean(mydata.B$Speed)
mydata.B.center$Power = mydata.B$Power - mean(mydata.B$Power)
mydata.B.center$Airbag = mydata.B$Airbag - mean(mydata.B$Airbag)
mydata.B.center$Temp = mydata.B$Temp - mean(mydata.B$Temp)
mydata.B.center$Wind = mydata.B$Wind - mean(mydata.B$Wind)

fit1.centered = lm(Severity ~ Year + Safety + Speed + Power + Airbag + Temp + Wind, 
                   data = mydata.B.center)
summary(fit1.centered)

fit1.AIC = stepAIC(fit1.centered, k = 2, scope = list(upper=~., lower=~1))
summary(fit1.AIC) # candidate 1
plot(fit1.AIC)
boxcox(fit1.AIC) # even though the residuals vs fitted produces a curved band, 
                  # no transformation of Y is recommended to correct. 
                  # therefore the model probably need more terms. => interactions

# model 2 ####################
# try squared terms
fit2 = lm(Severity ~ Year + I(Year^2) + Safety + I(Safety^2)+ Speed + I(Speed^2) + 
            Power + I(Power^2)+ Airbag + I(Airbag^2)+ Temp + I(Temp^2)
          + Wind + I(Wind^2) , data = mydata.B.center)
fit2.AIC = stepAIC(fit2, k = 2, 
                   scope = list(upper=~.,
                                lower=~Year + Safety + Speed + Power + Airbag))
summary(fit2.AIC) # takes out first order terms but keeps some second order
plot(fit2.AIC)# still no correction for res vs fitted

# model 3 ###########
# try just first order plus interactions
fit4 = lm(Severity ~ .*., data = mydata.B.center)
fit4.AIC = stepAIC(fit4, direction = "both", k = 2, 
                   scope = list(upper=~., lower=~1))

summary(fit4.AIC) 
# candidate 3, but can be reduced further by removing insignificant interaction terms
plot(fit4.AIC)

fit4a = lm(formula = Severity ~ Year + Safety + Speed + Power + Airbag + 
     Temp +  Year:Speed + Year:Power + Year:Airbag + Safety:Speed + 
     Safety:Temp + Speed:Power , data = mydata.B.center)
summary(fit4a)
plot(fit4a)

par(mfrow=c(1,2))
scatter.smooth(x=fit4a$fitted.values, y=fit4a$residuals,
               main = "Model 3 Residual Plot",
               xlab = "Fitted Values", ylab = "Residuals")
qqnorm(fit4a$residuals, main = "Model 3 Q-Q Plot")
qqline(fit4a$residuals)
par(mfrow=c(1,1))

tablefit1.AIC = data.frame("AIC" = AIC(fit1.AIC), 
                           "Adj.R2" = round(summary(fit1.AIC)$adj.r.squared,4),
                           "Number of Covariates" = length(fit1.AIC$assign), 
                           "degrees of freedom" = fit1.AIC$df.residual,
                           "Residual Std Error" = summary(fit1.AIC)$sigma)

tablefit2.AIC = data.frame("AIC" = AIC(fit2.AIC), 
                           "Adj.R2" = round(summary(fit2.AIC)$adj.r.squared,4),
                           "Number of Covariates" = length(fit2.AIC$assign),
                           "degrees of freedom" = fit2.AIC$df.residual,
                           "Residual Std Error" = summary(fit2.AIC)$sigma)

tablefit4a = data.frame("AIC" = AIC(fit4a), 
                        "Adj.R2" = round(summary(fit4a)$adj.r.squared,4),
                        "Number of Covariates" = length(fit4a$assign),
                        "degrees of freedom" = fit4a$df.residual,
                        "Residual Std Error" = summary(fit4)$sigma)
alltable = rbind(tablefit1.AIC, tablefit2.AIC, tablefit4a)
rownames(alltable) = c("Model 1", "Model 2", "Model 3")

#interpretation
newdata0 = data.frame(Year = min(mydata.B.center$Year),
                      Safety = 0, 
                      Speed=min(mydata.B.center$Speed), 
                      Power=min(mydata.B.center$Power), 
                      Airbag=min(mydata.B.center$Airbag),
                      Temp=0)

newdata1 = data.frame(Year = min(mydata.B.center$Year)+1,
                      Safety = 0, 
                      Speed=min(mydata.B.center$Speed), 
                      Power=min(mydata.B.center$Power), 
                      Airbag=min(mydata.B.center$Airbag),
                      Temp=0)
predict(fit4a, newdata = newdata1) - predict(fit4a, newdata = newdata0)

fit4a$coefficients[which(names(fit4a$coefficients)=="Year")]+
fit4a$coefficients[which(names(fit4a$coefficients)=="Year:Speed")]*min(mydata.B.center$Speed)+
fit4a$coefficients[which(names(fit4a$coefficients)=="Year:Power")]*min(mydata.B.center$Power)+
fit4a$coefficients[which(names(fit4a$coefficients)=="Year:Airbag")]*min(mydata.B.center$Airbag)

######
# B2
#####
design = mydata.B.center
attach(design)
design$Year2 = Year*Year
design$Year3 = Year*Year*Year
design$Safety2 = Safety*Safety
design$Safety3 = Safety*Safety*Safety
design$Speed2 = Speed*Speed
design$Speed3 = Speed*Speed*Speed
design$Power2 = Power*Power
design$Power3 = Power*Power*Power
design$Airbag2 = Airbag*Airbag
design$Airbag3 = Airbag*Airbag*Airbag
design$Temp2 = Temp*Temp
design$Temp3 = Temp*Temp*Temp
design$Wind2 = Wind*Wind
design$Wind3 = Wind*Wind*Wind
design$YearSafety = Year*Safety
design$YearSpeed = Year*Speed
design$YearPower = Year*Power
design$YearAirbag = Year*Airbag
design$YearTemp = Year*Temp
design$YearWind = Year*Wind
design$SafetySpeed = Safety*Speed
design$SafetyPower = Safety*Power
design$SafetyAirbag = Safety*Airbag
design$SafetyTemp = Safety*Temp
design$SafetyWind = Safety*Wind
design$SpeedPower = Speed*Power
design$SpeedAirbag = Speed*Airbag
design$SpeedTemp = Speed*Temp
design$SpeedWind = Speed*Wind
design$PowerAirbag = Power*Airbag
design$PowerTemp = Power*Temp
design$PowerWind = Power*Wind
design$AirbagTemp = Airbag*Temp
design$AirbagWind = Airbag*Wind
design$TempWind = Temp*Wind
detach(design)

#setup
set.seed(603723)
permut=sample(nrow(design))
myMSPE = matrix(200*10, nrow = 200, ncol = 10) # column is a fold
myCVerror = c()
lambda = seq(0.01, 2, 0.01)

#loop
for (k in 1:10){
#k=1  
  rownumbers = permut[(30*(k-1)+1):(30*k)] #length 30
  mytest = as.matrix(design[rownumbers,]) # size of 30
  mytraining = as.matrix(design[-rownumbers,]) # size 300 - 30
  
  # (lasso) alpha = 1
  myfit = glmnet(x = mytraining[,-8], y = mytraining[,8], alpha = 1, lambda = lambda)
  # every column corresponds to a lambda value: 30 x 200
  # predicting values of Y for every value of lambda
  mypredictions = predict.glmnet(myfit, newx = mytest[,-8])
   
  for(i in 1:length(lambda)){
    myMSPE[i,k] = 1/30*sum((mytest[,8]-mypredictions[,i])^2)
  }
}

myCVerror = data.frame(cbind(lambda = lambda, "CVerror" = rowMeans(myMSPE)))
# plot lambda v CV error
myCVerror[which(myCVerror$CVerror==min(myCVerror$CVerror)),] # 90.81642 lambda = 1.67
par(mfrow=c(1,1))
plot(x=myCVerror$lambda, y=myCVerror$CVerror, type = 'l', ylab = 'CV error', 
     xlab = expression(lambda))
abline(v=1.67, col = "blue")
text(1.4, 110, expression(lambda == 1.67), cex = 1.2)

# fit
glmnet.fit = glmnet(x=as.matrix(design[,-8]), y=as.matrix(design[,8]), alpha = 1)
coeffs = as.data.frame(as.matrix(coef(glmnet.fit, s = 1.67)))
coeffs$name = row.names(coeffs)
colnames(coeffs) = c("Estimates", "Parameters")
coeffs = coeffs[,c(2,1)]

nonzeroCoeffs = coeffs[which(coeffs$Estimates!=0),]
rownames(nonzeroCoeffs) = seq(1,nrow(nonzeroCoeffs),1)
nonzeroCoeffs$Estimates = format(nonzeroCoeffs$Estimates, digits = 4, scientific=T)
xtable(cbind(nonzeroCoeffs[1:6,], nonzeroCoeffs[7:12,]))
plot(glmnet.fit, label = T)

######
# B3
#####
fit.gam <- gam(Severity ~ s(Year) +  s(Safety) + s(Speed) + 
                 s(Power) + s(Airbag) + s(Temp) + s(Wind), 
               family = gaussian, data = mydata.B.center, method = "ML")

pobject = c(format(summary(fit.gam)$p.table[,1:3], digits=3), 
                format(summary(fit.gam)$p.table[,4], scientific = T, digits=5))

sobject = cbind(format(summary(fit.gam)$s.table[,1:3],digits=3), 
                format(summary(fit.gam)$s.table[,4], scientific = T, digits=5))
psobject = rbind("Intercept" = pobject, sobject)
colnames(psobject) = colnames(summary(fit.gam)$p.table)

fit.gam.2 <- gam(Severity ~ s(Year) +  s(Safety) + s(Speed) + 
                 s(Power) + s(Airbag) + s(Temp), 
               family = gaussian, data = mydata.B.center, method = "ML")

fit.gam.3 <- gam(Severity ~ s(Year) +  s(Safety) + s(Speed) + 
                   s(Power) + s(Airbag), 
                 family = gaussian, data = mydata.B.center, method = "ML")

fit.gam.4 <- gam(Severity ~ s(Year) +  s(Safety) + s(Speed) + 
                   s(Power), 
                 family = gaussian, data = mydata.B.center, method = "ML")

fit.gam.5 <- gam(Severity ~ s(Year) + Safety + Speed + 
                   Power  , 
                 family = gaussian, data = mydata.B.center, method = "ML")

fit.gam.6 <- gam(Severity ~ s(Year) + Safety + Speed + 
                   Power + Airbag  , 
                 family = gaussian, data = mydata.B.center, method = "ML")

# variable selection with AIC (gam)
#fit.gam
fit.gam$aic; fit.gam.2$aic; fit.gam.3$aic; fit.gam.4$aic; fit.gam.5$aic;fit.gam.6$aic

which.min(c(fit.gam$aic, fit.gam.2$aic, fit.gam.3$aic, fit.gam.4$aic, 
            fit.gam.5$aic, fit.gam.6$aic))

pobject = cbind(format(summary(fit.gam.6)$p.table[,1:3], digits=3), 
            format(summary(fit.gam.6)$p.table[,4], scientific = T, digits=5))
colnames(pobject) = colnames(summary(fit.gam)$p.table)

sobject = c(format(summary(fit.gam.6)$s.table[,1:3],digits=3), 
                format(summary(fit.gam.6)$s.table[,4], scientific = T, digits=5))
psobject = rbind(pobject, "s(Year)" = sobject)
colnames(psobject) = colnames(summary(fit.gam)$p.table)

# selected model:
# Severity ~ s(Year) + Safety + Speed + Power + Airbag
plot(fit.gam.6)
scatter.smooth(fit.gam.6$fitted.values, fit.gam.6$residuals)
qqnorm(fit.gam.6$residuals);qqline(fit.gam.6$residuals)

# polynomial model:
# Severity ~ Year + I(Year^2) + Safety + Speed + Power + Airbag
poly.fit = lm(Severity ~ Year + I(Year^2) + Safety + Speed + 
                Power + Airbag, data = mydata.B.center)
object = summary(poly.fit)$coefficients
str(summary(poly.fit))
table = cbind(format(object[,1:3], digits=2), format(object[,4],scientific=T, digits=5))
colnames(table) = colnames(summary(poly.fit)$coefficients)

plot(poly.fit)
AICpoly = AIC(poly.fit, k=2)
polyrsq = summary(poly.fit)$adj.r.squared
compare = data.frame("AIC" = c(AICgam5, AICpoly), "AdjR2" = c(gamrsq, polyrsq))
rownames(compare) = c("GAM Model", "Polynomial Model")
