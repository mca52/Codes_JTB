source("C:/Users/Admin/Desktop/Oxford/Oxford BDI/SetupRuns.R")
# Codes for a single Population 
# Set up population characteristics
# In this analysis considering only one population so popData will only have one row
nPops <- 1 #number of towns/settlements (set to 1 if not using movement functionality)
popSize <- 1500 #mean population size
popData <- generatePops(nPops,popSize)
popData$prev <- NA
popData$nation <- 1 #number of countries = 1
dMatrix <- matrix(30,1,1) # only one population, so just a zero matrix # 1500 per population, we consider 30 sites so the population size is that of an EU


# Up population characteristics to generate a new population
newPop <- Population$new(popData,dMatrix)

# Burn in period defaults to 100 years (1200 months) to allow prevalence to equilibriate
newPop$burnin()
Popburnin<- newPop$burnin()
# Retrieve population/prevalence statistics:
VtH <- newPop$VtH #Vector to host ratio
#prevalence <- length(which(newPop$Mf>0))/newPop$nHosts #calculate prevalence
aggk <- newPop$k #aggregation parameter k
N <- newPop$nHost #population size

##################
# Example for running annual MDA until EPHP threshold (1% mf) reached:
HALT=FALSE

callrunsMDA <- function(){
  # run with MDA until mf < 1%, storing monthly prevalence
  
  restart <- 0
  simPop <- newPop$clone(deep=T) # clone the population
  #simPop$AssignLowVectorVars()
  prevalence <-  (length(which((simPop$age/12)>0  & simPop$Mf>0))/length(which((simPop$age/12)>0 )))*100
  print(prevalence)
  ict_prevalence <- 0.97*(length(which((simPop$age/12)>0  & (simPop$WF+simPop$WM)>0)))/length(which((simPop$age/12)>0 ))*100
  prevs_ict<-ict_prevalence
  print(prevs_ict)
  mosql_prevalence <- length(which((simPop$age/12)>0 & exp(-simPop$larvae)>0 ))/length(which(simPop$age>0))
  mosq_larvae_prevs<- mosql_prevalence
  HALT <- FALSE
  month_ref <- 59 # start counting months just before MDA (MDA on month 12 in this example)
  month_total <- 1 # overall month count
  simPop$bnCov <- 0 # fractional proportion
  prevs <- prevalence # first entry is baseline prevalence
  maxTime <- 20# maximum number of years to run for
  MDArounds <- 0 # number of MDA rounds completed
  treatment = "IA" # Ivermection and albendazole
  MDAcov <- 0.65 # 65% coverage, can be any value from 0-1
  compliance = 0.2 # can be any value from 0-1
  threshold_prevalence <- 0.1
  prev_sequela <- (length(which((simPop$WF+simPop$WM)>0 & simPop$susceptibility>11))/length(which((simPop$age/12)>0 )))*100
  prevs_susceptibility<- prev_sequela # only change
  prev_sequelah <- (length(which((simPop$WF+simPop$WM)>0 & simPop$susceptibilityh>11))/length(which((simPop$age/12)>0 )))*100 # only change
  prevs_susceptibilityh<- prev_sequelah # only change
  no_surveys<-0
  EPHP<-FALSE  
  totalworms <- simPop$WM+simPop$WF
  #totalworms <- if(length(which(totalworms==0)!=0)) totalworms[-which(totalworms==0)]
  cumtotworms <- totalworms[1]
  dt <- 1
  prob_elim <- ifelse(prevs==0,1,0)
  morbidityl <- c()
  morbidityh <- c()
  a_death <- c()
  a_disability <- c()
  incident_death <- c()
  incident_cases <- c()
  incident_cases_h <- c()
  
  
  while(HALT == FALSE){
    a_disability=simPop$age
    simPop$runTimestep() # run a time step
    #simPop$simPop$AssignLowVectorVars()
    
    if(month_total==12|month_total==24|month_total==36|month_total==48|month_total==60){ # (Change here)if multiple rounds of MDA then do MDA
      simPop$runMDA(coverage=MDAcov, drug=treatment, compliance=compliance)
      MDArounds <- MDArounds +1 # record MDA round
    }
    if (month_total==66){
      print("Wait for 6 months")
    }
    mfprev <-  (length(which((simPop$age/12)>0  & simPop$Mf>0))/length(which((simPop$age/12)>0 )))*100
    month_total <- month_total +1
    prevs[month_total] <- mfprev
    ict_prevalence <- 0.97*(length(which((simPop$age/12)>=0  & (simPop$WF+simPop$WM)>0)))/length(which((simPop$age/12)>=0 ))*100
    totalworms<-(simPop$WM+simPop$WF)
    totalworms <- if(length(which(totalworms==0)!=0)) totalworms[-which(totalworms==0)]
    dt<- month_total-dt
    prevs_ict[month_total]<-ict_prevalence
    cumtotworms[month_total] <- 1+ ((totalworms)*dt)
    mosql_prevalence <- length(which((simPop$age/12)>0 & exp(-simPop$larvae)>0 ))/length(which(simPop$age>0))
    mosq_larvae_prevs[month_total] <- mosql_prevalence
    mfprevtas <-  (length(which((simPop$age/12)>=6 & (simPop$age)/12 <=7 & simPop$Mf>0))/length(which((simPop$age/12)>=6 & (simPop$age)/12 <=7)))*100
    ictprevtas <-  (length(which((simPop$age/12)>=6 & (simPop$age)/12 <=7 & ict_prevalence>0))/length(which((simPop$age/12)>=6 & (simPop$age)/12 <=7)))*100
    prev_sequela <- (length(which((simPop$WF+simPop$WM)>0 & simPop$susceptibility>11))/length(which((simPop$age/12)>0)))*100
    prevs_susceptibility[month_total]<- prev_sequela # only change
    prev_sequelah <- (length(which((simPop$WF+simPop$WM)>0 & simPop$susceptibilityh>11))/length(which((simPop$age/12)>0)))*100 # only change
    prevs_susceptibilityh[month_total]<- prev_sequelah # only change
    prob_elim[month_total]<-ifelse(prevs==0,1,0)
    incident_cases= c(incident_cases,sum(ifelse(prevs_susceptibility>11,1,0)))
    incident_cases_h= c(incident_cases_h,sum(ifelse(prevs_susceptibilityh>11,1,0)))
    incident_death= c(incident_death,sum(ifelse(prevs_susceptibility>11 & (simPop$age/12)>150,1,0)))
    a_death <- c(a_death,a_disability+simPop$age)
    
    if(mfprevtas <= threshold_prevalence && month_total==66){ # (Change here )# If <1% just before next round then stop
      EPHP <- TRUE # Assumed EPHP achieved
      HALT <- FALSE
    }
    if(month_total ==66){ # If exceed maxTime then stop # change here when you change the threshold prevalence
      
      HALT <- TRUE
    }
    #print(month_total)
    
    morbidityl=c(morbidityl,sum(ifelse(prevs_susceptibility>11,1,0))/simPop$nHost)
    morbidityh=c(morbidityh,sum(ifelse(prevs_susceptibilityh>11,1,0))/simPop$nHost)
  }
  if (EPHP==TRUE){
    no_surveys=no_surveys+1
  }
  HALT=FALSE
  no_surveys=no_surveys+1
  if (EPHP==FALSE){
    print("pre-TAS failed")
    while (HALT==FALSE){
      simPop$runTimestep() # run a time step
      if(month_total==72|month_total==84){ # (Change here)if multiple rounds of MDA then do MDA
        simPop$runMDA(coverage=MDAcov, drug=treatment, compliance=compliance)
        MDArounds <- MDArounds +1 # record MDA round
      }
      if (month_total==90){
        #print("Wait for 6 months")
      }
      mfprev <-  (length(which((simPop$age/12)>0  & simPop$Mf>0))/length(which((simPop$age/12)>0)))*100
      month_total <- month_total +1
      prevs[month_total] <- mfprev
      ict_prevalence <- 0.97*(length(which((simPop$age/12)>=0  & (simPop$WF+simPop$WM)>0)))/length(which((simPop$age/12)>=0))*100
      prevs_ict[month_total]<- ict_prevalence
      totalworms<-(simPop$WM+simPop$WF)
      prev_sequela <- (length(which((simPop$WF+simPop$WM)>0 & simPop$susceptibility>11))/length(which((simPop$age/12)>0)))*100
      prevs_susceptibility[month_total]<- prev_sequela # only change
      prev_sequelah <- (length(which((simPop$WF+simPop$WM)>0 & simPop$susceptibilityh>11))/length(which((simPop$age/12)>0)))*100 # only change
      prevs_susceptibilityh[month_total]<- prev_sequelah # only change
      
      dt<- month_total-dt
      cumtotworms[month_total] <- ((totalworms)*dt)
      mosql_prevalence <- length(which((simPop$age/12)>15 & exp(-simPop$larvae)>0 ))/length(which(simPop$age>0))
      mosq_larvae_prevs[month_total] <- mosql_prevalence
      mfprevtas <-  (length(which((simPop$age/12)>=20 & simPop$Mf>0))/length(which((simPop$age/12)>=20)))*100
      ictprevtas <-  (length(which((simPop$age/12)>=6 & (simPop$age)/12 <=7 & ict_prevalence>0))/length(which((simPop$age/12)>=6 & (simPop$age)/12 <=7)))*100
      morbidityl=c(morbidityl,sum(ifelse(prevs_susceptibility>11,1,0))/simPop$nHost)
      morbidityh=c(morbidityh,sum(ifelse(prevs_susceptibilityh>11,1,0))/simPop$nHost)
      incident_cases= c(incident_cases,sum(ifelse(prevs_susceptibility>11,1,0)))
      incident_death= c(incident_death,sum(ifelse(prevs_susceptibility>11 & (simPop$age/12)>11,1,0)))
      incident_cases_h= c(incident_cases_h,sum(ifelse(prevs_susceptibilityh>11,1,0)))
      a_death <- c(a_death,a_disability+simPop$age)
      
      if(mfprevtas <= threshold_prevalence && month_total==90 ){ # (Change here )# If <1% just before next round then stop
        EPHP <- TRUE # Assumed EPHP achieved
        
        HALT <- TRUE
      }
      if(mfprevtas >= threshold_prevalence && month_total == 90){ # If exceed maxTime then stop # change here when you change the threshold prevalence
        EPHP <- FALSE
        HALT <- TRUE
      }
      
    }
    
    print("TAS-1 failed, do 2 more rounds of MDA and check the TAS 6 months later") 
    #cat("TAS-1 MDA rounds=",MDArounds)
  } else { print("pre-TAS passed")
    while(HALT==FALSE){
      simPop$runTimestep() # run a time step
      if(month_total==102)
      {
        print("wait for 3 years")
      }
      
      mfprev <-  (length(which((simPop$age/12)>0  & simPop$Mf>0))/length(which((simPop$age/12)>0)))*100
      month_total <- month_total +1
      prevs[month_total] <- mfprev
      ict_prevalence <- 0.97*(length(which((simPop$age/12)>=0  & (simPop$WF+simPop$WM)>0)))/length(which((simPop$age/12)>=0))*100
      prevs_ict[month_total]<-ict_prevalence
      totalworms<-(simPop$WM+simPop$WF)
      prev_sequela <- (length(which((simPop$WF+simPop$WM)>0 & simPop$susceptibility>11))/length(which((simPop$age/12)>0)))*100
      prevs_susceptibility[month_total]<- prev_sequela # only change
      prev_sequelah <- (length(which((simPop$WF+simPop$WM)>0 & simPop$susceptibilityh>11))/length(which((simPop$age/12)>0)))*100 # only change
      prevs_susceptibilityh[month_total]<- prev_sequelah # only change
      dt<- month_total-dt
      
      cumtotworms[month_total] <- ((totalworms)*dt)
      l3 <- length(which((simPop$age/12)>5 & simPop$L3>0 ))/length(which(simPop$age>=0))
      mosql_prevalence <- length(which((simPop$age/12)>0 & exp(-simPop$larvae)>0 ))/length(which(simPop$age>0))
      mosq_larvae_prevs[month_total] <- mosql_prevalence
      mfprevtas <-  (length(which((simPop$age/12)>=20 & simPop$Mf>0))/length(which((simPop$age/12)>=20)))*100
      ictprevtas <-  (length(which((simPop$age/12)>=6 & (simPop$age)/12 <=7 & ict_prevalence>0))/length(which((simPop$age/12)>=6 & (simPop$age)/12 <=7)))*100
      morbidityl=c(morbidityl,sum(ifelse(prevs_susceptibility>11,1,0))/simPop$nHost)
      morbidityh=c(morbidityh,sum(ifelse(prevs_susceptibilityh>11,1,0))/simPop$nHost)
      incident_cases= c(incident_cases,sum(ifelse(prevs_susceptibility>11,1,0)))
      incident_death= c(incident_death,sum(ifelse(prevs_susceptibility>11 & (simPop$age/12)>0,1,0)))
      incident_cases_h= c(incident_cases_h,sum(ifelse(prevs_susceptibilityh>11,1,0)))
      a_death <- c(a_death,a_disability+simPop$age)
      
      if(mfprevtas <= threshold_prevalence && month_total==102 ){ # (Change here )# If <1% just before next round then stop 
        EPHP <- TRUE # Assumed EPHP achieved
        
        HALT <- TRUE
      }
      if(mfprevtas >= threshold_prevalence && month_total == 102){ # If exceed maxTime then stop # change here when you change the threshold prevalence
        EPHP <- FALSE
        HALT <- TRUE
      }
    }
    print("TAS-1 passed, wait for 3 years")
    
    if (EPHP==FALSE){
      print("Need to restart the MDA")
      restart=restart+1
    }
  }
  #TAS-2
  no_surveys<-no_surveys+1
  #print(EPHP)
  #print(month_total)
  #print(EPHP)
  #print(mfprevtas)
  #print(threshold_prevalence)
  month_total1<-month_total+12
  month_total2<-month_total1+12
  month_total5<-month_total+36
  HALT=FALSE
  if (EPHP==FALSE) {
    
    while(HALT==FALSE){
      simPop$runTimestep() # run a time step
      if(month_total==month_total1|month_total==month_total2){ # (Change here)if multiple rounds of MDA then do MDA
        simPop$runMDA(coverage=MDAcov, drug=treatment, compliance=compliance)
        MDArounds <- MDArounds +1 # record MDA round
      }
      if (month_total==month_total2+6){
        
        print("Wait for 6 months")
      }
      mfprev <-  (length(which((simPop$age/12)>0  & simPop$Mf>0))/length(which((simPop$age/12)>0 )))*100
      month_total <- month_total +1
      prevs[month_total] <- mfprev
      ict_prevalence <- 0.97*(length(which((simPop$age/12)>0 & (simPop$WF+simPop$WM)>0)))/length(which((simPop$age/12)>0))*100
      prevs_ict[month_total]<-ict_prevalence
      
      totalworms<-(simPop$WM+simPop$WF)
      mfprevtas <-  (length(which((simPop$age/12)>=20 & simPop$Mf>0))/length(which((simPop$age/12)>=20)))*100
      ictprevtas <-  (length(which((simPop$age/12)>=6 & (simPop$age)/12 <=7 & ict_prevalence>0))/length(which((simPop$age/12)>=6 & (simPop$age)/12 <=7)))*100
      #totalworms <- if(length(which(totalworms==0)!=0)) totalworms[-which(totalworms==0)]
      dt<- month_total-dt
      cumtotworms[month_total] <- ((totalworms)*dt)
      mosql_prevalence <- length(which((simPop$age/12)>5 & exp(-simPop$larvae)>0 ))/length(which(simPop$age>0))
      mosq_larvae_prevs[month_total] <- mosql_prevalence
      prev_sequela <- (length(which((simPop$WF+simPop$WM)>0 & simPop$susceptibility>11))/length(which((simPop$age/12)>5 )))*100
      prevs_susceptibility[month_total]<- prev_sequela # only change
      prev_sequelah <- (length(which((simPop$WF+simPop$WM)>0 & simPop$susceptibilityh>11))/length(which((simPop$age/12)>5 )))*100 # only change
      prevs_susceptibilityh[month_total]<- prev_sequelah # only change
      morbidityl=c(morbidityl,sum(ifelse(prevs_susceptibility>11,1,0))/simPop$nHost)
      morbidityh=c(morbidityh,sum(ifelse(prevs_susceptibilityh>11,1,0))/simPop$nHost)
      incident_cases= c(incident_cases,sum(ifelse(prevs_susceptibility>11,1,0)))
      incident_death= c(incident_death,sum(ifelse(prevs_susceptibility>11 & (simPop$age/12)>0,1,0)))
      incident_cases_h= c(incident_cases_h,sum(ifelse(prevs_susceptibilityh>11,1,0)))
      a_death <- c(a_death,a_disability+simPop$age)
      
      if(mfprevtas <= threshold_prevalence && month_total==month_total2+6 ){ # (Change here )# If <1% just before next round then stop
        EPHP <- TRUE # Assumed EPHP achieved
        
        HALT <- TRUE
      }
      if(mfprevtas >= threshold_prevalence && month_total == month_total2+6){ # If exceed maxTime then stop # change here when you change the threshold prevalence
        EPHP <- FALSE
        HALT <- TRUE
      }
    }
    print("TAS-2 failed , do 2 more rounds of MDA and wait 6 months for TAS") 
    #cat("TAS-2 MDA rounds=",MDArounds)
  } else { while(HALT==FALSE){
    simPop$runTimestep() # run a time step
    if(month_total==month_total5)
    {
      print("wait for 3 years")
    }
    
    mfprev <-  (length(which((simPop$age/12)>15  & simPop$Mf>0))/length(which((simPop$age/12)>5 )))*100
    month_total <- month_total +1
    prevs[month_total] <- mfprev
    ict_prevalence <- 0.97*(length(which((simPop$age/12)>5  & (simPop$WF+simPop$WM)>0)))/length(which((simPop$age/12)>5 ))*100
    prevs_ict[month_total]<-ict_prevalence
    totalworms<-(simPop$WM+simPop$WF)
    mfprevtas <-  (length(which((simPop$age/12)>=20  & simPop$Mf>0))/length(which((simPop$age/12)>=20)))*100
    ictprevtas <-  (length(which((simPop$age/12)>=6 & (simPop$age)/12 <=7 & ict_prevalence>0))/length(which((simPop$age/12)>=6 & (simPop$age)/12 <=7)))*100
    prev_sequela <- (length(which((simPop$WF+simPop$WM)>0 & simPop$susceptibility>11))/length(which((simPop$age/12)>5 )))*100
    prevs_susceptibility[month_total]<- prev_sequela # only change
    prev_sequelah <- (length(which((simPop$WF+simPop$WM)>0 & simPop$susceptibilityh>11))/length(which((simPop$age/12)>5 )))*100 # only change
    prevs_susceptibilityh[month_total]<- prev_sequelah # only change
    #morbidity[month_total]=sum(ifelse(prev_sequela>11,1,0))/1200
    #totalworms <- if(length(which(totalworms==0)!=0)) totalworms[-which(totalworms==0)]
    dt<- month_total-dt
    #print(totalworms)
    cumtotworms[month_total] <- ((totalworms)*dt)
    mosql_prevalence <- length(which((simPop$age/12)>5 & exp(-simPop$larvae)>0 ))/length(which(simPop$age>0))
    mosq_larvae_prevs[month_total] <- mosql_prevalence
    prev_sequelah<- (simPop$WF+simPop$WM)*simPop$susceptibilityh # only change
    prevs_susceptibilityh[month_total] <- prev_sequelah # only change
    morbidityl=c(morbidityl,sum(ifelse(prevs_susceptibility>11,1,0))/simPop$nHost)
    morbidityh=c(morbidityh,sum(ifelse(prevs_susceptibilityh>11,1,0))/simPop$nHost)
    incident_cases= c(incident_cases,sum(ifelse(prevs_susceptibility>11,1,0)))
    incident_death= c(incident_death,sum(ifelse(prevs_susceptibility>11 & (simPop$age/12)>50,1,0)))
    incident_cases_h= c(incident_cases_h,sum(ifelse(prevs_susceptibilityh>11,1,0)))
    a_death <- c(a_death,a_disability+simPop$age)
    
    if(mfprevtas <= threshold_prevalence && month_total==month_total5 ){ # (Change here )# If <1% just before next round then stop
      EPHP <- TRUE # Assumed EPHP achieved
      
      HALT <- TRUE
    }
    if(mfprevtas >= threshold_prevalence && month_total ==month_total5){ # If exceed maxTime then stop
      EPHP <- FALSE
      HALT <- TRUE
    }
  }
    print("TAS-2 passed, wait for 3 years")
    if (EPHP==FALSE){
      print("Need to restart the MDA")
      restart=restart+1
    }
  }
  #TAS-3
  no_surveys <- no_surveys+1
  month_total3<-month_total+12
  month_total4<-month_total3+12
  month_total6<-month_total+36
  #print(month_total3)
  #print(month_total4)
  
  HALT=FALSE
  if (EPHP==FALSE){
    while(HALT==FALSE){
      simPop$runTimestep() # run a time step
      if(month_total==month_total3|month_total==month_total4){ # (Change here)if multiple rounds of MDA then do MDA
        simPop$runMDA(coverage=MDAcov, drug=treatment, compliance=compliance)
        MDArounds <- MDArounds +1 # record MDA round
      }
      if (month_total==month_total4+6)
      {
        print("wait for 6 months!!")
      }
      mfprev <-  (length(which((simPop$age/12)>0  & simPop$Mf>0))/length(which((simPop$age/12)>0)))*100
      month_total <- month_total +1
      prevs[month_total] <- mfprev
      ict_prevalence <- 0.97*(length(which((simPop$age/12)>0 & (simPop$WF+simPop$WM)>0)))/length(which((simPop$age/12)>0 ))*100
      prevs_ict[month_total]<-ict_prevalence
      totalworms<-(simPop$WM+simPop$WF)
      #totalworms <- if(length(which(totalworms==0)!=0)) totalworms[-which(totalworms==0)]
      dt<- month_total-dt
      cumtotworms[month_total] <- ((totalworms)*dt)
      prev_sequela <- (length(which((simPop$WF+simPop$WM)>0 & simPop$susceptibility>11))/length(which((simPop$age/12)>0 )))*100
      prevs_susceptibility[month_total]<- prev_sequela # only change
      prev_sequelah <- (length(which((simPop$WF+simPop$WM)>0 & simPop$susceptibilityh>11))/length(which((simPop$age/12)>0 )))*100 # only change
      prevs_susceptibilityh[month_total]<- prev_sequelah # only change
      mosql_prevalence <- length(which((simPop$age/12)>5 & exp(-simPop$larvae)>0 ))/length(which(simPop$age>0))
      mosq_larvae_prevs[month_total] <- mosql_prevalence
      mfprevtas <-  (length(which((simPop$age/12)>=20 & simPop$Mf>0))/length(which((simPop$age/12)>=20)))*100
      ictprevtas <-  (length(which((simPop$age/12)>=6 & (simPop$age)/12 <=7 & ict_prevalence>0))/length(which((simPop$age/12)>=6 & (simPop$age)/12 <=7)))*100
      morbidityl=c(morbidityl,sum(ifelse(prevs_susceptibility>11,1,0))/simPop$nHost)
      morbidityh=c(morbidityh,sum(ifelse(prevs_susceptibilityh>11,1,0))/simPop$nHost)
      incident_cases= c(incident_cases,sum(ifelse(prevs_susceptibility>11,1,0)))
      incident_death= c(incident_death,sum(ifelse(prevs_susceptibility>11 & (simPop$age/12)>0,1,0)))
      incident_cases_h= c(incident_cases_h,sum(ifelse(prevs_susceptibilityh>11,1,0)))
      a_death <- c(a_death,a_disability+simPop$age)
      
      if(mfprevtas <= threshold_prevalence && month_total==maxTime*12 ){ # (Change here )# If <1% just before next round then stop
        EPHP <- TRUE # Assumed EPHP achieved
        
        HALT <- TRUE
      }
      if(mfprevtas >= threshold_prevalence && month_total >= maxTime*12){ # If exceed maxTime then stop
        EPHP <- FALSE
        HALT <- TRUE
      }
      
    }
    print("TAS-3 failed , do 2 more rounds of MDA and wait 6 months for TAS !!")
    #cat("TAS-3 MDA rounds=",MDArounds)
  } else { while(HALT==FALSE){
    simPop$runTimestep() # run a time step
    if(month_total==maxTime*12)
    {
      
    }
    
    mfprev <-  (length(which((simPop$age/12)>0  & simPop$Mf>0))/length(which((simPop$age/12)>0 )))*100
    month_total <- month_total +1
    prevs[month_total] <- mfprev
    ict_prevalence <- 0.97*(length(which((simPop$age/12)>0  & (simPop$WF+simPop$WM)>0)))/length(which((simPop$age/12)>0))*100
    prevs_ict[month_total]<-ict_prevalence
    
    
    totalworms<-(simPop$WM+simPop$WF)
    
    prev_sequela <- (length(which((simPop$WF+simPop$WM)>0 & simPop$susceptibility>11))/length(which((simPop$age/12)>0)))*100
    prevs_susceptibility[month_total]<- prev_sequela # only change
    prev_sequelah <- (length(which((simPop$WF+simPop$WM)>0 & simPop$susceptibilityh>11))/length(which((simPop$age/12)>0)))*100 # only change
    prevs_susceptibilityh[month_total]<- prev_sequelah # only change
    morbidityl=c(morbidityl,sum(ifelse(prevs_susceptibility>11,1,0))/simPop$nHost)
    morbidityh=c(morbidityh,sum(ifelse(prevs_susceptibilityh>11,1,0))/simPop$nHost)
    dt<- month_total-dt
    cumtotworms[month_total] <- ((totalworms)*dt)
    mosql_prevalence <- length(which((simPop$age/12)>0 & exp(-simPop$larvae)>0 ))/length(which(simPop$age>0))
    mosq_larvae_prevs[month_total] <- mosql_prevalence
    mfprevtas <-  (length(which((simPop$age/12)>=20  & simPop$Mf>0))/length(which((simPop$age/12)>=20)))*100
    ictprevtas <-  (length(which((simPop$age/12)>=6 & (simPop$age)/12 <=7 & ict_prevalence>0))/length(which((simPop$age/12)>=6 & (simPop$age)/12 <=7)))*100
    incident_cases= c(incident_cases,sum(ifelse(prevs_susceptibility>11,1,0)))
    incident_death= c(incident_death,sum(ifelse(prevs_susceptibility>11 & (simPop$age/12)>0,1,0)))
    incident_cases_h= c(incident_cases_h,sum(ifelse(prevs_susceptibilityh>11,1,0)))
    a_death <- c(a_death,a_disability+simPop$age)
    
    if(mfprevtas <= threshold_prevalence && month_total==maxTime*12 ){ # (Change here )# If <1% just before next round then stop
      EPHP <- TRUE # Assumed EPHP achieved
      
      HALT <- TRUE
    }
    if(mfprevtas >= threshold_prevalence && month_total ==maxTime*12){ # If exceed maxTime then stop
      EPHP <- FALSE
      HALT <- TRUE
    }
    
    
  }
    print("TAS-3 passed")
    if (EPHP==FALSE){
      print("Need to restart the MDA")
      restart=restart+1
    }
  }
  
  return(list(first=month_total,second=EPHP,third=mfprev,fourth=prevs,fifth=MDArounds,sixth=no_surveys,seventh=cumtotworms,eighth=prevs_ict,ninth=prevs_susceptibility,tenth=restart,eleven=mosq_larvae_prevs,twelve=prevs_susceptibilityh,thirteen=prob_elim,fourteen=morbidityl,fifteen=morbidityh,sixteen=a_death,seventeen=a_disability,eighteen=incident_cases,nineteen=incident_death,twenty=incident_cases_h))
}




sample.R <- function(N){
  p <- ggplot()
  q<-ggplot()
  r<-ggplot()
  s<-ggplot()
  t<-ggplot()
  u<-ggplot()
  v<-ggplot()
  w<-ggplot()
  x<-ggplot()
  y<-ggplot()
  z<-ggplot()
  b<-ggplot()
  dat <- vector("list",N)
  dat1 <- vector("list",N)
  dat2<-vector("list",N)
  dat3<-vector("list",N)
  dat4<-vector("list",N)
  dat5<-vector("list",N)
  dat6<-vector("list",N)
  dat7<-vector("list",N)
  dat8<-vector("list",N)
  dat9<-vector("list",N)
  dat10<-vector("list",N)
  dat11<-vector("list",N)
  dat12<-vector("list",N)
  ans1 <- c()
  ans2 <- c()
  ans3 <- c()
  ans4 <- c()
  ans5 <- c()
  ans6 <- c()
  ans7<- c()
  ans8<-c()
  ans9<-c()
  ans10<-c()
  ans11<-c()
  elim <- c()
  ans11 <-c()
  ans12<-c()
  ans13<-c()
  ans14<-c()
  ans15 <- c()
  ans16 <- c()
  for (j in 1:N){
    ans=callrunsMDA()
    dat[[j]]= tibble(time=(1:ans$first)/12,mf_prevalence=ans$fourth)
    dat1[[j]]= tibble(time=(1:ans$first)/12,cum_worm_yrs=log(cumsum(ans$seventh)))
    dat2[[j]]= tibble(time=(1:ans$first)/12,ict_prevalence=ans$eighth)
    dat3[[j]]= tibble(time=(1:ans$first)/12,mosql_prevalence=ans$eleven)
    dat11[[j]]= tibble(time=(1:ans$first)/12,prevalence_l=ans$ninth)
    dat12[[j]]= tibble(time=(1:ans$first)/12,prevalence_h=ans$twelve)
    
    ans1 <- c(ans1, ans$seventh) # cum_worm_yrs
    ans2 <- c(ans2,ans$fifth) #cum_total_MDA_rounds
    ans3 <- c(ans3,ans$sixth) #cum_no_surveys
    ans4=c(ans4,ans$fourth) # cum_prevs_vector
    ans5=c(ans5,ans$first)# cum_month_total
    ans13=c(ans13,ans$ninth) # morbidity lymphodema
    ans14=c(ans14,ans$tweleve)# morbidity hydrocele
    
    
    cost=ans$fifth*7640.92 +ans$sixth*12494.75 
    inf_averted=cost/mean(cumsum(ans$seventh))
    ans12<-c(ans12,inf_averted) # cost per infection averted
    check1<-ifelse(ans$seventh==0,1,0)
    check2<-ifelse(ans$fourth==0,1,0)
    prob_elim<- sum(ifelse(check1==1 & check2==1,1,0))/(100*12)
    dat4[[j]]=tibble(prob_elim=prob_elim,MDArounds=ans$fifth)
    dat5[[j]]=tibble(prob_elim=prob_elim,surveys=ans$sixth)
    dat10[[j]]=tibble(inc=inf_averted,wormyrs=(ans$seven))
    
    
    line_types <- c("LINE1"=1,"LINE2"=3,"LINE3"=2)
    
    p <- p +
      geom_line(data = dat[[j]], mapping = aes(x = time, y = mf_prevalence, color = "Infection(mf) Prevalence", linetype = "Infection(mf) Prevalence")) +
      geom_step(data = dat1[[j]], mapping = aes(x = time, y = cum_worm_yrs, color = "Cumulative Wormburden", linetype = "Cumulative Wormburden")) +
      geom_line(data = dat11[[j]], mapping = aes(x = time, y = prevalence_l, color = "DALY burden", linetype = "DALY burden")) +
      labs(x = 'Time since start of the treatment (years)', y = "Proportion of simulations (%)") +
      scale_linetype_manual(name = "Legend", values = c("Infection(mf) Prevalence" = "solid", "Cumulative Wormburden" = "dotted", "DALY burden" = "dashed")) +
      scale_colour_manual(name = "Legend", values = c("Infection(mf) Prevalence" = "darkred", "Cumulative Wormburden" = "darkgreen", "DALY burden" = "darkblue"))
    
    q <- q +
      geom_line(data = dat2[[j]], mapping = aes(x = time, y = ict_prevalence, color = "Infection(Ag) Prevalence", linetype = "Infection(Ag) Prevalence")) +
      geom_line(data = dat11[[j]], mapping = aes(x = time, y = prevalence_l, color = "DALY burden", linetype = "DALY burden")) +
      labs(x = 'Time since start of the treatment (years)', y = "Proportion of simulations (%)") +
      scale_linetype_manual(name = "Legend", values = c("Infection(Ag) Prevalence" = "solid", "DALY burden" = "dashed")) +
      scale_colour_manual(name = "Legend", values = c("Infection(Ag) Prevalence" = "darkred", "DALY burden" = "darkblue"))
    
    dat6[[j]]= tibble(wormyrs=mean(cumsum(ans$seventh)),MDArounds=ans$fifth)
    dat7[[j]]= tibble(time=(1:ans$first)/12,mosql_prev=ans$eleven)
    dat8[[j]]=tibble(wormyrs=mean(cumsum(ans$seventh)),surveys=ans$sixth)
    dat9[[j]]=tibble(wormyrs=mean(cumsum(ans$seventh)),costs=cost)
    #s = s + geom_point(data = dat3[[j]], mapping = aes(x = costs, y = prob_elim), col = j)
    t = t + geom_point(data = dat4[[j]], mapping = aes(y = prob_elim, x = MDArounds), col = j)
    u = u + geom_point(data = dat5[[j]], mapping = aes(y = prob_elim, x = surveys), col = j)
    v = v + geom_point(data = dat6[[j]], mapping = aes(y = wormyrs, x =MDArounds), col = j)
    x = x + geom_point(data = dat8[[j]], mapping = aes(y = wormyrs, x =surveys), col = j)
    y = y + geom_point(data = dat9[[j]], mapping = aes(x=costs,y = wormyrs), col = j)
    w = w + geom_line(data = dat7[[j]], mapping = aes(x = time, y = mosql_prev), col = j)
    b= b+ geom_point(data = dat10[[j]], mapping = aes(x = inc, y = wormyrs), col = j)
    ans7=c(ans7,cost)
    ans8=c(ans8,ans$ninth)  # mosq prevalence
    ans9=c(ans9,prob_elim)
    ans10=c(ans10,mean(cumsum(ans$seventh)))
    ans11 =c(ans11,ans$tenth) # restart the MDA
    print((mean(ans$sixteen/12)-mean(ans$seventeen/12))*0.109*sum(ans$eighteen)/100+(mean(ans$sixteen/12)-80)*sum(ans$nineteen)/100)
    ans15=c(ans15,(mean(ans$sixteen/12)-mean(ans$seventeen/12))*0.109*sum(ans$eighteen)/100)
    ans16=c(ans16,(mean(ans$sixteen/12)-mean(ans$seventeen/12))*0.128*sum(ans$twenty)/100)
    
    print(j)
  }
  return(list(first=p,second=ans5,third=ans4,fourth=ans2,fifth=ans3,sixth=ans1,seventh=ans7,eighth=ans13,ninth=ans14,tenth=s,eleven=t,twelve=u,thirteen=v,fourteen=w,fifteen=x,sixteen=y,seventeen=ans9,eighteen=ans10,nineteen=elim,twenty=ans11))
  
  
}
sample<-sample.R(1000)

#-------Outcomes----

check1<-ifelse(sample$sixth==0,1,0)
check2<-ifelse(sample$third==0,1,0)
prob_elim<- (ifelse(check1==1 & check2==1,1,0))
cat("The prob of elimination=",sum(sample$nineteen)/1000)
cat("The prob of restart=",sum(sample$twenty)/1000)
#cat("Probability of person with lymphodema=",sample$twentythree)
#cat("Probability of person with hydrocele=",sample$twentyfour)
# write the output into text files
cat(sample$fourth, file = "MDA_0.02_6_7_5_10.txt",sep=",")
cat(sample$eighteen ,file="wormyrs_0.02_6_7_5_10.txt",sep=",")
cat(sample$seventh,file="cost_0.02_6_7_5_10.txt",sep=",")
cat(sample$twentythree,file="lymphodema_0.02_6_7_5_10.txt",sep=",")
cat(sample$twentyfour,file="hydrocele_0.02_6_7_5_10.txt",sep=",")

sample$first
a<-c(1:50)
b<-24*a
sample$third[b]
prob_elim<- sum(ifelse(sample$third[b]<=1,1,0)/50)

