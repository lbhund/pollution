library("dplyr")
library("MASS")
library(pscl)
library(xtable)

# LOAD DATA --------------

load("Data/dat.Rda")
dat$bedperpop <- dat$beds/dat$population
dat$logpop <- log(dat$population)
dat$pmq <- cut(dat$mean_pm25, breaks=c(quantile(dat$mean_pm25, probs = seq(0, 1, by = 0.20))))
subdat <-  dat[dat$state != "NY",]
subdat2 <-  dat[dat$Deaths > 0,]
subdat3 <-  dat[dat$Deaths > 10 & dat$state != "NY",]
subdat4 <-  dat[dat$Deaths > 0 & dat$state != "NY",]




# LINEAR PREDICTOR -----------------

linpred_nooffset <- formula(Deaths ~ mean_pm25 + factor(q_popdensity)
+ scale(poverty)  + scale(log(medianhousevalue))
+ scale(log(medhouseholdincome)) + scale(pct_owner_occ) 
+ scale(education) + scale(pct_blk) + scale(hispanic)
+ scale(older_pecent) + scale(prime_pecent) + scale(mid_pecent) 
+ scale(date_since_social) + scale(date_since)
+ scale(beds/population) 
+ scale(obese) + scale(smoke)
+ scale(mean_summer_temp) + scale(mean_winter_temp) 
+ scale(mean_summer_rm) + scale(mean_winter_rm))

linpred <- update(linpred_nooffset, . ~ . + offset(log(population)))

# summary stats -------------
dat$deathrate <- 100000*dat$Deaths/dat$population
dat$bedpopscale <- dat$bedperpop*100000
dat$medhouseholdincomes <- dat$medhouseholdincome/1000
dat$medianhousevalues <- dat$medianhousevalue/1000
dat$smokes <- dat$smoke*100
dat$obeses <- dat$obese*100
dat$povertys <- dat$poverty*100
dat$hispanics <- dat$hispanic*100
dat$pct_blks <- dat$pct_blk*100
dat$older_pecents <- dat$older_pecent*100
dat$prime_pecents <- dat$prime_pecent*100
dat$mid_pecents <- dat$mid_pecent*100
dat$educations <- dat$education*100
dat$pct_owner_occs <- dat$pct_owner_occ*100
varnames <- c("deathrate",
              "mean_pm25", 
              "bedpopscale", 
              "date_since",
              "date_since_social",
              "smokes",
              "obeses",
              "povertys",
              "educations",
              "pct_owner_occs",
              "hispanics",
              "pct_blks",
              "older_pecents",
              "mid_pecents",
              "prime_pecents",
              "popdensity",
              "medhouseholdincomes",
              "medianhousevalues",
              "mean_summer_temp_F",
              "mean_winter_temp_F",
              "mean_summer_rm",
              "mean_winter_rm")

#convert temp to fahrenheit
dat$mean_summer_temp_F <- dat$mean_summer_temp*9/5 - 459.67
dat$mean_winter_temp_F <- dat$mean_winter_temp*9/5 - 459.67

demographicstable <- read.csv("Data/demographicstable.txt", header=FALSE)

demographicstable[,c("V4", "V5")] <- round(t(sapply(varnames, function(aa) 
                                      c(mean(dat[,aa]), sd(dat[,aa])))), digits=1)
colnames(demographicstable) <- c("", "Mean", "SD", "Mean", "SD")
print(xtable(demographicstable, label = "tab_demo", caption="Mean and standard deviation of variables in analysis
             for original data (left) and my attempt to recreate the data (right)."), file="LH/Table3.tex", include.rownames=F)






# GLM MODELS -----------------

makemod <- function(dat, name){
  out <- list()
  out$name <- name
  
  mod1 = glm(linpred, data = dat, family='quasipoisson')
  out$c1 <- exp(coef(mod1)[2])
  out$fit1 <-   sd(exp(predict(mod1)) - dat$Deaths)
  
  mod2 <- glm.nb(linpred,  data = dat)
  out$c2 <- exp(coef(mod2)[2])
  out$fit2 <- sd(exp(predict(mod2)) - dat$Deaths)
  return(out)
}

subdat <-  dat[dat$state != "NY",]
subdat2 <-  dat[dat$Deaths > 0,]
subdat3 <-  dat[dat$Deaths > 0 & dat$state != "NY",]
subdat4 <-  dat[dat$Deaths > 10 & dat$state != "NY",]

m1 <- makemod(dat, "All")
m2 <- makemod(subdat, "No NY")
m3 <- makemod(subdat2, "> 0")
m4 <- makemod(subdat3, ">0 & no NY")
m5 <- makemod(subdat4, ">10, no NY")

mod5 <- zeroinfl(linpred_nooffset, data = dat, offset =logpop, dist="poisson")
c5 <- exp(coef(mod5))[2]
n5 <- "Zero P: all"
mod6 <- zeroinfl(linpred_nooffset, data = dat, offset =logpop, dist="negbin")
c6 <- exp(coef(mod6))[2]
n6 <- "Zero NB: all"

mod7 <- zeroinfl(linpred_nooffset, data = subdat, offset =logpop, dist="poisson")
c7 <- exp(coef(mod7))[2]
n7 <- "Zero P: no NY"
mod8 <- zeroinfl(linpred_nooffset, data = subdat, offset =logpop, dist="negbin")
c8 <- exp(coef(mod8))[2]
n8 <- "Zero NB: no NY"

#table of model MRRs
names <- c("Zero: all", "Zero: no NY", m1$name, m2$name, m3$name, m4$name, m5$name)
estsp  <- c(c5, c7, m1$c1, m2$c1, m3$c1, m4$c1, m5$c1)
estsnb  <- c(c6, c8, m1$c2, m2$c2, m3$c2, m4$c2, m5$c2)
print(xtable(data.frame(Analysis=names, Poisson=estsp, NegBin=estsnb), label = "tab_sensitivity",
             caption="Sensitivity analysis comparing models on different subsets of data."), include.rownames=F, file="LH/Table2.tex")





# MAKE TABLE OF RESULTS ---------------------
rownames.vec <- c("PM 2.5 (ug/m3)",
                  "Population density (Q2)",
                  "Population density (Q3)",
                  "Population density (Q4)",
                  "Population density (Q5)",
                  "Pct Poverty",
                  "log(Median house value)",
                  "log(Median household income)",
                  "Pct Owner-occupied housing",
                  "Pct Less than high school education",
                  "Pct Black",
                  "Pct Hispanic",
                  "Pct ≥65 years of age",
                  "Pct 15-44 years of age",
                  "Pct 45-64 years of age",
                  "Days since stay-at-home order",
                  "Days since first case",
                  "Rate of hospital beds",
                  "Pct Obese",
                  "Pct Smokers",
                  "Average summer temperature (F)",
                  "Average winter temperature (F)",
                  "Average summer relative humidity (Pct)",
                  "Average winter relative humidity (Pct)")
original.ests <- c(1.08,
0.86,
0.58,
0.47,
0.52,
1.02,
1.17,
1.28,
1.12,
1.36,
1.45,
1.00,
1.15,
0.93,
0.96,
1.28,
2.96,
1.12,
0.94,
1.08,
0.96,
1.18,
0.84,
1.00)
nparam <- length(exp(coef(mod8)))/2
table.params <- cbind(exp(coef(mod7)[2:nparam]), exp(coef(mod8)[2:nparam]), original.ests)
rownames(table.params) <- rownames.vec
colnames(table.params) <- c("Poisson", "Negative Bin", "Original")
print(xtable(table.params, label = "tab_coefs", caption="Mortality rate ratios for all variables in the main analysis; my zero-inflated models
      using the Poisson and negative binomial are in the first two columns (fit to data excluding NY state); the original 
      results from April 22 are in the far right column."), include.rownames=T, include.colnames=T, file="LH/Table1.tex")






# RESIDUAL PLOT --------------

mod1 = glm(linpred, data = subdat, family='quasipoisson')
mod2 <- glm.nb(linpred,  data = subdat)

#residual plot
  mus <- exp(predict(mod1))
  resids <- (mus - subdat$Deaths)
  mus2 <- exp(predict(mod2))
  resids2 <- (mus2 - subdat$Deaths)
  col1 <- rgb(.1,.5,.5,.8)
  col2 <- rgb(.6,.2,.4,.8)
  
  png(filename="LH/Figure_Residuals.png", width=480, height=480)
  plot(subdat$Deaths, resids, cex=.7, log="", col=col1, pch=16, 
       xlab="Observed", ylab="(Observed - fitted)", cex.lab=1.5, cex.axis=1.5,
       ylim=range(c(resids, resids2)), xlim=range(c(mus, mus2)))
  points(subdat$Deaths, resids2, cex=.7, col=col2, pch=16)
    legend("bottomleft", pch=c(16, 16), c("Quasi-poisson", "Neg. Bin."),
           col=c(col1, col2), cex=1.3)
  dev.off()
  
# predicted variance versus residual variance  
  model.pred <- mus*summary(mod1)$dispersion
  model.pred2 <- mus2 + mus2^2/mod2$theta
  
  png(filename="LH/Figure_Variance.png", width=480, height=480)
  plot(mus, model.pred, pch=16, ylab="Predicted model variance", 
       xlab="Predicted mean",cex.lab=1.5, cex.axis=1.5,
       ylim=c(0, 300000), xlim=range(c(mus, mus2)), col=col1)
  points(mus2, model.pred2, pch=16, col=col2)
  
  breakpoints <- c(0, 5, 30, 75, 150, 350, 1000)
  subdat$cats <- cut(subdat$Deaths, breaks=c(-1, 0, 10, 50, 100, 200, 500, 1500))
  varpts1 <- tapply(resids, subdat$cats, var)
  varpts2 <- tapply(resids2, subdat$cats, var)
  points(breakpoints, varpts1, col=col1, pch=1, cex=4)
  points(breakpoints, varpts2, col=col2, pch=1, cex=4)
  dev.off()
  
  #weights
  png(filename="LH/Figure_Weights.png", width=480, height=480)
  plot(mus, mus/summary(mod1)$dispersion, pch=15, col=col1,cex.lab=1.5, cex.axis=1.5,
       xlab="mean", ylab="weight", cex.axis=2, cex.lab=2)
  points(mus2, mus2/(1+ mus2/mod2$theta), pch=15, col=col2)
  dev.off()
  
  


# CHECK LINEAR MODEL -----------------
subdat$deathrate <- subdat$Deaths/subdat$population
lm(deathrate ~ mean_pm25 + factor(q_popdensity) + scale(poverty) + 
     scale(log(medianhousevalue)) + scale(log(medhouseholdincome)) + 
     scale(pct_owner_occ) + scale(education) + scale(pct_blk) + 
     scale(hispanic) + scale(older_pecent) + scale(prime_pecent) + 
     scale(mid_pecent) + scale(date_since_social) + scale(date_since) + 
     scale(beds/population) + scale(obese) + scale(smoke) + scale(mean_summer_temp) + 
     scale(mean_winter_temp) + scale(mean_summer_rm) + scale(mean_winter_rm) + 
     population, data=subdat, weights=mus/summary(mod1)$dispersion)

lm(deathrate ~ mean_pm25 + factor(q_popdensity) + scale(poverty) + 
     scale(log(medianhousevalue)) + scale(log(medhouseholdincome)) + 
     scale(pct_owner_occ) + scale(education) + scale(pct_blk) + 
     scale(hispanic) + scale(older_pecent) + scale(prime_pecent) + 
     scale(mid_pecent) + scale(date_since_social) + scale(date_since) + 
     scale(beds/population) + scale(obese) + scale(smoke) + scale(mean_summer_temp) + 
     scale(mean_winter_temp) + scale(mean_summer_rm) + scale(mean_winter_rm) + 
     population, data=subdat, weights=mus2/(1+ mus2/mod2$theta))




# CHECK CATEGORICAL predictor ----------------
linpred2 <- update(linpred, .~. - mean_pm25 + pmq)
moda = glm(linpred2, data = dat, family='quasipoisson')
modb = glm(linpred2, data = subdat, family='quasipoisson')
modc = glm(linpred2, data = subdat2, family='quasipoisson')
modd = glm(linpred2, data = subdat3, family='quasipoisson')
mode = glm(linpred2, data = subdat4, family='quasipoisson')
modd = glm.nb(linpred2, data = subdat)
summary(modd)
summary(modb)