rm(list=ls(all=TRUE))

#LOAD PACKAGES
library(crayon)
library(backports)
library(vctrs)
library(pillar)
library(dplyr)
library(withr)
library(ggplot2)
library(bayesAB)
library(labeling)
library(matrixStats)
library(magrittr)
library(ggpubr)
library(ggsci)
library(digest)
library(lubridate)
library(survival)
library(plyr)


#SET PARAMETERS
Number.Runs = 200 #Number of runs to perform of stochastic model. Model parameters are constant across runs.
Run.Num.Output = 20 #Output progress every X runs
Flu.Season.Length = 166 #Last day to include cases. 166 = March 31st, 0 = October 16th (earliest flu cases across all years reported in Ray et al. CID 2019)
Non.Flu.N = 43000 #Total number of non-flu cases to create data for
Flu.N = 7000 #Total number of flu cases to create data for
Flu.Vacc.Start = -76 #First day of flu vaccination. -76 = August 1st
Flu.Vacc.End = 166 #Last day of flu vaccination. 166 = March 31st
Flu.Vacc.Delay = 14 #Exclude people who become cases within this many days after receiving the flu vaccine
Flu.Vacc.Cutoff = 45 #Exclude people vaccinated on or after this date. 45 = December 1st
Prop.Early.Vacc = 0.2 #Proportion of flu cases that were vaccinated earlier in the season than average
Early.Vacc.Days = 20 #Number days earlier that the 'Prop.Early.Vacc' participants were vaccinated. NOTE: The average vaccination time of flu cases will typically be brought forward by less than Prop.Early.Vacc*Early.Vacc.Days due to application of the exclusion criteria used by Ray et al. CID 2019
Analysis.Start = 0 #Exclude participants who became cases before this time

Vacc.Mean.Log.Norm = 4.36 #Mean parameter of the log-normal distribution for the dates of influenza vaccination
Vacc.SD.Log.Norm = 0.453 #Standard deviation parameter of the log-normal distribution for the dates of influenza vaccination
Non.Flu.Mean.Log.Norm = 4.7 #Mean parameter of the log-normal distribution for the dates of non-influenza medically-attended illness
Non.Flu.SD.Log.Norm = 0.35 #Standard deviation parameter of the log-normal distribution for the dates of non-influenza medically-attended illness
Flu.Mean.Log.Norm = 4.75 #Mean parameter of the log-normal distribution for the dates of influenza medically-attended illness
Flu.SD.Log.Norm = 0.32 #Standard deviation parameter of the log-normal distribution for the dates of influenza medically-attended illness

dec.places = 3 #Number of decimal places for rounding results


#INITIALIZE VARIABLES
Full.Max.Weekly.Pct.Flu = rep(NaN, Number.Runs)
Number.Pcpts = matrix(nrow = Number.Runs, ncol = 2); colnames(Number.Pcpts) = c('Flu', 'Non.Flu')

Full.Case.Time.Means = matrix(nrow = Number.Runs, ncol = 2); colnames(Full.Case.Time.Means) = c('Flu', 'Non.Flu')
Full.Case.Time.Medians = matrix(nrow = Number.Runs, ncol = 2); colnames(Full.Case.Time.Medians) = c('Flu', 'Non.Flu')
Full.Case.Time.Modes = matrix(nrow = Number.Runs, ncol = 2); colnames(Full.Case.Time.Modes) = c('Flu', 'Non.Flu')

Full.Vacc.Time.Means = matrix(nrow = Number.Runs, ncol = 2); colnames(Full.Vacc.Time.Means) = c('Flu', 'Non.Flu')
Full.Vacc.Time.Medians = matrix(nrow = Number.Runs, ncol = 2); colnames(Full.Vacc.Time.Medians) = c('Flu', 'Non.Flu')
Full.Vacc.Time.Modes = matrix(nrow = Number.Runs, ncol = 2); colnames(Full.Vacc.Time.Modes) = c('Flu', 'Non.Flu')

Full.Vacc.Case.Lag.Means = matrix(nrow = Number.Runs, ncol = 2); colnames(Full.Vacc.Case.Lag.Means) = c('Flu', 'Non.Flu')
Full.Vacc.Case.Lag.Medians = matrix(nrow = Number.Runs, ncol = 2); colnames(Full.Vacc.Case.Lag.Medians) = c('Flu', 'Non.Flu')
Full.Vacc.Case.Lag.Modes = matrix(nrow = Number.Runs, ncol = 2); colnames(Full.Vacc.Case.Lag.Modes) = c('Flu', 'Non.Flu')


Analysis.Case.Time.Means = matrix(nrow = Number.Runs, ncol = 2); colnames(Analysis.Case.Time.Means) = c('Flu', 'Non.Flu')
Analysis.Case.Time.Medians = matrix(nrow = Number.Runs, ncol = 2); colnames(Analysis.Case.Time.Medians) = c('Flu', 'Non.Flu')
Analysis.Case.Time.Modes = matrix(nrow = Number.Runs, ncol = 2); colnames(Analysis.Case.Time.Modes) = c('Flu', 'Non.Flu')

Analysis.Vacc.Time.Means = matrix(nrow = Number.Runs, ncol = 2); colnames(Analysis.Vacc.Time.Means) = c('Flu', 'Non.Flu')
Analysis.Vacc.Time.Medians = matrix(nrow = Number.Runs, ncol = 2); colnames(Analysis.Vacc.Time.Medians) = c('Flu', 'Non.Flu')
Analysis.Vacc.Time.Modes = matrix(nrow = Number.Runs, ncol = 2); colnames(Analysis.Vacc.Time.Modes) = c('Flu', 'Non.Flu')

Analysis.Vacc.Case.Lag.Means = matrix(nrow = Number.Runs, ncol = 2); colnames(Analysis.Vacc.Case.Lag.Means) = c('Flu', 'Non.Flu')
Analysis.Vacc.Case.Lag.Medians = matrix(nrow = Number.Runs, ncol = 2); colnames(Analysis.Vacc.Case.Lag.Medians) = c('Flu', 'Non.Flu')
Analysis.Vacc.Case.Lag.Modes = matrix(nrow = Number.Runs, ncol = 2); colnames(Analysis.Vacc.Case.Lag.Modes) = c('Flu', 'Non.Flu')

Month.VE.Crude = rep(NaN, Number.Runs)
Month.VE.Month.Adj = rep(NaN, Number.Runs)
Month.VE.Week.Adj = rep(NaN, Number.Runs)
Month.VE.Cond.Linear = rep(NaN, Number.Runs)
Month.VE.Cond.Cat = matrix(nrow = Number.Runs, ncol = 5)


#DEFINE FUNCTION FOR CALCULATING A MODE
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


t0 <- proc.time() #Start timer to track model running time
cat('Number.Runs =', Number.Runs, '\n')
for (i in 1:Number.Runs) {
  
  #Set flu vaccination dates
  Flu.Vacc.Times = rlnorm(Flu.N + Non.Flu.N, meanlog = Vacc.Mean.Log.Norm, sdlog = Vacc.SD.Log.Norm) + Flu.Vacc.Start
  
  #Set case dates for non-influenza medically-attended infection
  Non.Flu.Case.Times = as.data.frame(floor(rlnorm(n = Non.Flu.N, meanlog = Non.Flu.Mean.Log.Norm, sdlog = Non.Flu.SD.Log.Norm)))
  Non.Flu.Case.Times = cbind(Non.Flu.Case.Times, 'Non.Flu')
  colnames(Non.Flu.Case.Times) = c('Case.Time', 'Flu.Case')
  
  #Set case dates for influenza medically-attended infection
  Flu.Case.Times = as.data.frame(floor(rlnorm(n = Flu.N, meanlog = Flu.Mean.Log.Norm, sdlog = Flu.SD.Log.Norm)))
  Flu.Case.Times = cbind(Flu.Case.Times, 'Flu')
  colnames(Flu.Case.Times) = c('Case.Time', 'Flu.Case')
  
  #Combine simulated data on dates of vaccination, non-flu infection, and flu infection
  Combined.Case.Times = rbind(Flu.Case.Times, Non.Flu.Case.Times)
  Full.Data = cbind(Flu.Vacc.Times, Combined.Case.Times)
  colnames(Full.Data) = c('Vacc.Time', 'Case.Time', 'Flu.Case')
  
  
  #Shift vaccination time earlier in season for a proportion of influenza cases
  Flu.Cases = which(Full.Data[,'Flu.Case'] == 'Flu')
  Full.Data[Flu.Cases,'Vacc.Time'] = Full.Data[Flu.Cases,'Vacc.Time'] - sample(c(0,Early.Vacc.Days), size = length(Flu.Cases), replace = TRUE, prob = c(1-Prop.Early.Vacc, Prop.Early.Vacc)) #Include unmeasured confounder of earlier vaccination for flu cases
  
  #Create variable that stores the time between vaccination and becoming a case
  Full.Data = cbind(Full.Data, Full.Data[,'Case.Time'] - Full.Data[,'Vacc.Time'])
  colnames(Full.Data) = c('Vacc.Time', 'Case.Time', 'Flu.Case', 'Vacc.Case.Lag')
  
  #Calculate mean, median, and mode describing vaccination and case dates and lag between them for all individuals
  Full.Vacc.Time.Means[i,] = tapply(Full.Data[,'Vacc.Time'], Full.Data[,'Flu.Case'], mean)
  Full.Vacc.Time.Medians[i,] = tapply(Full.Data[,'Vacc.Time'], Full.Data[,'Flu.Case'], median)
  Full.Vacc.Time.Modes[i,] = tapply(Full.Data[,'Vacc.Time'], Full.Data[,'Flu.Case'], getmode)
  
  Full.Case.Time.Means[i,] = tapply(Full.Data[,'Case.Time'], Full.Data[,'Flu.Case'], mean)
  Full.Case.Time.Medians[i,] = tapply(Full.Data[,'Case.Time'], Full.Data[,'Flu.Case'], median)
  Full.Case.Time.Modes[i,] = tapply(Full.Data[,'Case.Time'], Full.Data[,'Flu.Case'], getmode)
  
  Full.Vacc.Case.Lag.Means[i,] = tapply(Full.Data[,'Vacc.Case.Lag'], Full.Data[,'Flu.Case'], mean)
  Full.Vacc.Case.Lag.Medians[i,] = tapply(Full.Data[,'Vacc.Case.Lag'], Full.Data[,'Flu.Case'], median)
  Full.Vacc.Case.Lag.Modes[i,] = tapply(Full.Data[,'Vacc.Case.Lag'], Full.Data[,'Flu.Case'], getmode)
  
  
  #Calculate the weekly % of medically-attended caases that were flu cases
  Weekly.Pct.Flu = vector(mode = "numeric", length = ceiling(Flu.Season.Length/7))
  for (j in 1:ceiling(Flu.Season.Length/7)) {
    Week = which(Full.Data[,'Case.Time'] >= j*7-7 & Full.Data[,'Case.Time'] < j*7)
    Weekly.Flu.Cases = length(which(Full.Data[Week,'Flu.Case'] == 'Flu'))
    Weekly.Pcpts = sum(table(Full.Data[Week,'Case.Time']))
    Weekly.Pct.Flu[j] = 100*Weekly.Flu.Cases/Weekly.Pcpts
  }
  Full.Max.Weekly.Pct.Flu[i] = max(Weekly.Pct.Flu, na.rm=T)
  
  
  #Determine which individuals meet the inclusion criteria used by Ray et al. (Clin Infect Dis 2019, 68(10):1623-1630) + create analysis dataset
  Pcpts = which(Full.Data[,'Case.Time'] <= Flu.Season.Length & Full.Data[,'Case.Time'] >= (Full.Data[,'Vacc.Time'] + Flu.Vacc.Delay) & Full.Data[,'Vacc.Time'] <= Flu.Vacc.Cutoff & Full.Data[,'Case.Time'] >= Analysis.Start) # 
  Analysis.Data.Frame = as.data.frame(Full.Data[Pcpts,])
  Number.Pcpts[i,'Flu'] = length(which(Analysis.Data.Frame[,'Flu.Case'] == 'Flu'))
  Number.Pcpts[i,'Non.Flu'] = length(which(Analysis.Data.Frame[,'Flu.Case'] == 'Non.Flu'))
  
  #Create month variable to use to adjust for calendar time
  for (j in 1:ceiling(Flu.Season.Length/30)) {
    Month = which(Analysis.Data.Frame[,'Case.Time'] >= j*30-30 & Analysis.Data.Frame[,'Case.Time'] < j*30)
    Analysis.Data.Frame[Month,5] = j
  }
  Analysis.Data.Frame[,5] = factor(Analysis.Data.Frame[,5])
  
  
  #Create week variable to use to adjust for calendar time
  for (j in 1:ceiling(Flu.Season.Length/7)) {
    Week = which(Analysis.Data.Frame[,'Case.Time'] >= j*7-7 & Analysis.Data.Frame[,'Case.Time'] < j*7)
    Analysis.Data.Frame[Week,6] = j
  }
  Analysis.Data.Frame[,6] = factor(Analysis.Data.Frame[,6])
  
  
  #Create categorical variable for the lag between vaccination and case dates
  Analysis.Data.Frame[,dim(Analysis.Data.Frame)[2]+1] = NA
  Cat.1 = which(Analysis.Data.Frame[,'Vacc.Case.Lag'] >= 14 & Analysis.Data.Frame[,'Vacc.Case.Lag'] <= 41)
  Cat.2 = which(Analysis.Data.Frame[,'Vacc.Case.Lag'] >= 42 & Analysis.Data.Frame[,'Vacc.Case.Lag'] <= 69)
  Cat.3 = which(Analysis.Data.Frame[,'Vacc.Case.Lag'] >= 70 & Analysis.Data.Frame[,'Vacc.Case.Lag'] <= 97)
  Cat.4 = which(Analysis.Data.Frame[,'Vacc.Case.Lag'] >= 98 & Analysis.Data.Frame[,'Vacc.Case.Lag'] <= 125)
  Cat.5 = which(Analysis.Data.Frame[,'Vacc.Case.Lag'] >= 126 & Analysis.Data.Frame[,'Vacc.Case.Lag'] <= 153)
  Cat.6 = which(Analysis.Data.Frame[,'Vacc.Case.Lag'] >= 154)
  
  Analysis.Data.Frame[Cat.1,dim(Analysis.Data.Frame)[2]] = '1'
  Analysis.Data.Frame[Cat.2,dim(Analysis.Data.Frame)[2]] = '2'
  Analysis.Data.Frame[Cat.3,dim(Analysis.Data.Frame)[2]] = '3'
  Analysis.Data.Frame[Cat.4,dim(Analysis.Data.Frame)[2]] = '4'
  Analysis.Data.Frame[Cat.5,dim(Analysis.Data.Frame)[2]] = '5'
  Analysis.Data.Frame[Cat.6,dim(Analysis.Data.Frame)[2]] = '6'
  
  colnames(Analysis.Data.Frame) = c('Vacc.Time', 'Case.Time', 'Flu.Case', 'Vacc.Case.Lag', 'Month.Factor', 'Week.Factor', 'Cat.Factor')
  
  
  #Calculate mean, median, and mode describing vaccination and case dates and lag between them for the analysis dataset
  Analysis.Vacc.Time.Means[i,] = tapply(Analysis.Data.Frame[,'Vacc.Time'], Analysis.Data.Frame[,'Flu.Case'], mean)
  Analysis.Vacc.Time.Medians[i,] = tapply(Analysis.Data.Frame[,'Vacc.Time'], Analysis.Data.Frame[,'Flu.Case'], median)
  Analysis.Vacc.Time.Modes[i,] = tapply(Analysis.Data.Frame[,'Vacc.Time'], Analysis.Data.Frame[,'Flu.Case'], getmode)
  
  Analysis.Case.Time.Means[i,] = tapply(Analysis.Data.Frame[,'Case.Time'], Analysis.Data.Frame[,'Flu.Case'], mean)
  Analysis.Case.Time.Medians[i,] = tapply(Analysis.Data.Frame[,'Case.Time'], Analysis.Data.Frame[,'Flu.Case'], median)
  Analysis.Case.Time.Modes[i,] = tapply(Analysis.Data.Frame[,'Case.Time'], Analysis.Data.Frame[,'Flu.Case'], getmode)

  Analysis.Vacc.Case.Lag.Means[i,] = tapply(Analysis.Data.Frame[,'Vacc.Case.Lag'], Analysis.Data.Frame[,'Flu.Case'], mean)
  Analysis.Vacc.Case.Lag.Medians[i,] = tapply(Analysis.Data.Frame[,'Vacc.Case.Lag'], Analysis.Data.Frame[,'Flu.Case'], median)
  Analysis.Vacc.Case.Lag.Modes[i,] = tapply(Analysis.Data.Frame[,'Vacc.Case.Lag'], Analysis.Data.Frame[,'Flu.Case'], getmode)
  
  
  #Create logical dependent variable for regression models
  Analysis.Data.Frame$Flu.Outcome = Analysis.Data.Frame$Flu.Case == "Flu"
  
  
  #Analyze simulated data using multiple approaches
  Model.Crude <- glm(Flu.Outcome ~ Vacc.Case.Lag, family=binomial(link='logit'), data=Analysis.Data.Frame)
  Month.VE.Crude[i] = exp(28*Model.Crude$coefficients[2])
  
  Model.Month.Adj<- glm(Flu.Outcome ~ Vacc.Case.Lag + Month.Factor, family=binomial(link='logit'), data=Analysis.Data.Frame)
  Month.VE.Month.Adj[i] = exp(28*Model.Month.Adj$coefficients[2])
  
  Model.Week.Adj<- glm(Flu.Outcome ~ Vacc.Case.Lag + Week.Factor, family=binomial(link='logit'), data=Analysis.Data.Frame)
  Month.VE.Week.Adj[i] = exp(28*Model.Week.Adj$coefficients[2])
  
  Model.Cond.Logit.Linear = clogit(Flu.Outcome ~ Vacc.Case.Lag + strata(Case.Time), Analysis.Data.Frame, method = "efron")
  Month.VE.Cond.Linear[i] = exp(28*Model.Cond.Logit.Linear$coefficients[1])
  
  Model.Cond.Logit.Categ = clogit(Flu.Outcome ~ Cat.Factor + strata(Case.Time), Analysis.Data.Frame, method = "efron")
  Month.VE.Cond.Cat[i,] = exp(Model.Cond.Logit.Categ$coefficients[1:5])
  
  
  #Output status of model runs
  if (i%%Run.Num.Output == 0 | i%%Number.Runs == 0) {
    cat(i," out of ",Number.Runs, "runs\ttime =", proc.time() - t0, "\n")
    
    #Create summary plots of simulated data if there is only 1 parameter set
    if (Number.Runs == 1) {
      pgamma(41, shape = Non.Flu.Shape, scale = Non.Flu.Scale, lower.tail = TRUE, log.p = FALSE)
      qgamma(0.975, Non.Flu.Shape, scale = Non.Flu.Scale, lower.tail = TRUE, log.p = FALSE)
      pgamma(41, shape = Flu.Shape, scale = Flu.Scale, lower.tail = TRUE, log.p = FALSE)
      qgamma(0.975, Flu.Shape, scale = Flu.Scale, lower.tail = TRUE, log.p = FALSE)
      
      plot(Weekly.Pct.Flu)
      
      #Plot histogram of case dates stratified by flu vs. non-flu case status
      p.Case.Time<-ggplot(Analysis.Data.Frame, aes(x=Case.Time, fill=Flu.Case, color=Flu.Case)) +
        geom_histogram(position="identity", alpha=0.5, breaks=seq(13,max(Analysis.Data.Frame$Case.Time),7))
      # Add mean lines
      mu.case <- ddply(Analysis.Data.Frame, "Flu.Case", summarise, grp.mean=mean(Case.Time))
      head(mu.case)
      p.Case.Time+geom_vline(data=mu.case, aes(xintercept=grp.mean, color=Flu.Case), linetype="dashed")
      
      #Plot histogram of vaccination dates stratified by flu vs. non-flu case status
      p.Vacc.Time<-ggplot(Analysis.Data.Frame, aes(x=Vacc.Time, fill=Flu.Case, color=Flu.Case)) +
        geom_histogram(position="identity", alpha=0.5, breaks=seq(Flu.Vacc.Start,Flu.Vacc.End,7))
      # Add mean lines
      mu.vacc <- ddply(Analysis.Data.Frame, "Flu.Case", summarise, grp.mean=mean(Vacc.Time))
      head(mu.vacc)
      p.Vacc.Time+geom_vline(data=mu.vacc, aes(xintercept=grp.mean, color=Flu.Case), linetype="dashed")
      
      #Plot histogram of time between vaccination date and case date stratified by flu vs. non-flu case status
      p.Vacc.Case.Lag<-ggplot(Analysis.Data.Frame, aes(x=Vacc.Case.Lag, fill=Flu.Case, color=Flu.Case)) +
        geom_histogram(position="identity", alpha=0.5, breaks=seq(13,max(Analysis.Data.Frame$Vacc.Case.Lag),7))
      # Add mean lines
      mu.lag <- ddply(Analysis.Data.Frame, "Flu.Case", summarise, grp.mean=mean(Vacc.Case.Lag))
      head(mu.lag)
      p.Vacc.Case.Lag+geom_vline(data=mu.lag, aes(xintercept=grp.mean, color=Flu.Case), linetype="dashed")
      
      #Plot histogram of time between vaccination date and case time stratified by period (as per Ray et al. CID 2019)
      p.Vacc.Case.Lag.Cat<-ggplot(Analysis.Data.Frame, aes(x=Vacc.Case.Lag, fill=Cat.Factor, color=Cat.Factor)) +
        geom_histogram(position="identity", alpha=0.5, breaks=seq(13,max(Analysis.Data.Frame$Vacc.Case.Lag),7))
      # Add mean lines
      mu.cat <- ddply(Analysis.Data.Frame, "Cat.Factor", summarise, grp.mean=mean(Vacc.Case.Lag))
      head(mu.cat)
      p.Vacc.Case.Lag.Cat+geom_vline(data=mu.cat, aes(xintercept=grp.mean, color=Cat.Factor), linetype="dashed")
      
      aggregate(Vacc.Case.Lag ~ Flu.Case + Cat.Factor, data = Analysis.Data.Frame, FUN = mean)
      aggregate(Vacc.Case.Lag ~ Flu.Case + Cat.Factor, data = Analysis.Data.Frame, FUN = median)
      aggregate(Vacc.Case.Lag ~ Flu.Case + Cat.Factor, data = Analysis.Data.Frame, FUN = getmode)
    } #End of summary plots
  } #End of output status of model runs
} #End for loop over parameter sets


#OUTPUT SUMMARY OF RUNS
#Output means, medians, and modes for vaccination and case dates and lag between them
median(Full.Vacc.Time.Means[,'Flu']); median(Full.Vacc.Time.Means[,'Non.Flu'])
median(Full.Vacc.Time.Medians[,'Flu']); median(Full.Vacc.Time.Medians[,'Non.Flu'])
median(Full.Vacc.Time.Modes[,'Flu']); median(Full.Vacc.Time.Modes[,'Non.Flu'])

median(Full.Case.Time.Means[,'Flu']); median(Full.Case.Time.Means[,'Non.Flu'])
median(Full.Case.Time.Medians[,'Flu']); median(Full.Case.Time.Medians[,'Non.Flu'])
median(Full.Case.Time.Modes[,'Flu']); median(Full.Case.Time.Modes[,'Non.Flu'])

median(Full.Vacc.Case.Lag.Means[,'Flu']); median(Full.Vacc.Case.Lag.Means[,'Non.Flu'])
median(Full.Vacc.Case.Lag.Medians[,'Flu']); median(Full.Vacc.Case.Lag.Medians[,'Non.Flu'])
median(Full.Vacc.Case.Lag.Modes[,'Flu']); median(Full.Vacc.Case.Lag.Modes[,'Non.Flu'])


median(Analysis.Vacc.Time.Means[,'Flu']); median(Analysis.Vacc.Time.Means[,'Non.Flu'])
median(Analysis.Vacc.Time.Medians[,'Flu']); median(Analysis.Vacc.Time.Medians[,'Non.Flu'])
median(Analysis.Vacc.Time.Modes[,'Flu']); median(Analysis.Vacc.Time.Modes[,'Non.Flu'])

median(Analysis.Case.Time.Means[,'Flu']); median(Analysis.Case.Time.Means[,'Non.Flu'])
median(Analysis.Case.Time.Medians[,'Flu']); median(Analysis.Case.Time.Medians[,'Non.Flu'])
median(Analysis.Case.Time.Modes[,'Flu']); median(Analysis.Case.Time.Modes[,'Non.Flu'])

median(Analysis.Vacc.Case.Lag.Means[,'Flu']); median(Analysis.Vacc.Case.Lag.Means[,'Non.Flu'])
median(Analysis.Vacc.Case.Lag.Medians[,'Flu']); median(Analysis.Vacc.Case.Lag.Medians[,'Non.Flu'])
median(Analysis.Vacc.Case.Lag.Modes[,'Flu']); median(Analysis.Vacc.Case.Lag.Modes[,'Non.Flu'])


#Summarize numbers of flu and non-flu participants across runs
round(quantile(Number.Pcpts[,'Flu'], probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), na.rm=TRUE))
round(quantile(Number.Pcpts[,'Non.Flu'], probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), na.rm=TRUE))
round(quantile(Full.Max.Weekly.Pct.Flu, probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), na.rm=TRUE), digits = dec.places)


#Summarize results of simulated studies quantifying waning of vaccine effectiveness
round(quantile(Month.VE.Crude, probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), na.rm=TRUE), digits = dec.places)
round(quantile(Month.VE.Month.Adj, probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), na.rm=TRUE), digits = dec.places)
round(quantile(Month.VE.Week.Adj, probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), na.rm=TRUE), digits = dec.places)
round(quantile(Month.VE.Cond.Linear, probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), na.rm=TRUE), digits = dec.places)
round(colMedians(Month.VE.Cond.Cat, na.rm=TRUE), digits = dec.places)
round(quantile(Month.VE.Cond.Cat[,1], probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), na.rm=TRUE), digits = dec.places)
round(quantile(Month.VE.Cond.Cat[,2], probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), na.rm=TRUE), digits = dec.places)
round(quantile(Month.VE.Cond.Cat[,3], probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), na.rm=TRUE), digits = dec.places)
round(quantile(Month.VE.Cond.Cat[,4], probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), na.rm=TRUE), digits = dec.places)
round(quantile(Month.VE.Cond.Cat[,5], probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), na.rm=TRUE), digits = dec.places)
