# k is the array job number supplied to the HPC, refers to scenario number
# These two lines are specific to the HPC, the rest are not
kenv <- Sys.getenv("SGE_TASK_ID")
k <- as.numeric(kenv)

#############################################################
# Scenarios dataframes
#############################################################
Nscenario <- 144

# Scenarios dataframe is created. Note that this is for when I>C; code for C>I can be obtained by swapping round ch_coef and ind_coef
scenarios <- data.frame(scenario_no = 1:Nscenario, no_clusters=numeric(Nscenario), clust_size_m=numeric(Nscenario),  
                        clust_pop_k=numeric(Nscenario), design=numeric(Nscenario), tp=numeric(Nscenario), 
                        md="MCAR", maxfup = 78, icc=0.1, cac=0.5, iac=0.5, sigma=1, ch_coef=-0.001285, ind_coef=-0.003855, beta1=-7.05, 
                        beta2=NA, beta3=NA, nonlin=1)
scenarios$no_clusters <- rep(c(44,32,30), c(72,48,24))
scenarios$clust_size_m <- rep(c(15,50,100), c(72,48,24))
scenarios$clust_pop_k <- rep(c(15,50,100,50,100,100), each=24)
alls <- rep(rep(c("cc","cs","oc"),4), times=2) 
subs <- rep(c(rep(c("cc","oc"), times=4),rep("cs",4)),2)
scenarios$design <- c(alls, subs, subs, alls, subs, alls)
allz <- rep(rep(c(2,3,5,20),2), each=3) 
subz <- rep(c(rep(c(2,3,5,20), each=2),c(2,3,5,20)),2)
scenarios$tp <- c(allz, subz, subz,allz,subz,allz)

#############################################################
# Analysis code in batches
#############################################################

library(dplyr)
library(lme4)
library(lmerTest) 

Nsim <- 650
finalresults <- data.frame(scenario = k, avg_ind_estimate = numeric(1),
                           avg_ind_se = numeric(1), avg_ind_se2 = numeric(1), ind_power = numeric(1),
                           avg_ch_estimate = numeric(1),
                           avg_ch_se = numeric(1), avg_ch_se2 = numeric(1), ch_power = numeric(1),
                           avg_total_estimate = numeric(1), 
                           avg_total_se = numeric(1), avg_total_se2 = numeric(1),
                           ind_emp_se = numeric(1), ch_emp_se = numeric(1), total_emp_se = numeric(1),
                           oc52emp_se = numeric(1), oc26emp_se = numeric(1),  # These are in OC estimand analysis code only
                           oc52_mc_se_empse = numeric(1), oc26_mc_se_empse = numeric(1),  # These are in OC estimand analysis code only
                           oc52_mc_se_bias = numeric(1), oc26_mc_se_bias = numeric(1),  # These are in OC estimand analysis code only
                           ind_mc_se_bias = numeric(1), ch_mc_se_bias = numeric(1), total_mc_se_bias = numeric(1),
                           error_count2 = numeric(1),msg_count2 = numeric(1),
                           ind_mc_se_power = numeric(1), ch_mc_se_power = numeric(1),
                           ind_mc_se_empse = numeric(1), ch_mc_se_empse = numeric(1), total_mc_se_empse = numeric(1),
                           avg_estimate = numeric(1), power = numeric(1),
                           avg_se = numeric(1), avg_se2 = numeric(1), emp_se = numeric(1),
                           mc_se_bias = numeric(1), error_count = numeric(1),msg_count = numeric(1),
                           mc_se_power = numeric(1), mc_se_empse = numeric(1))

# Create final results data frame which will have 144 rows for the 144 scenarios of Chapter 5
flag <- integer()
msgflag <- integer()
flag2 <- integer()
msgflag2 <- integer()

# Manipulation of the results dataframe
# comb52 and comb26 are in OC estimand analysis code only
results <- data.frame(simno = 1:Nsim, ind.estimate = numeric(Nsim), ind.se = numeric(Nsim), ind.sesq = numeric(Nsim), ind.pval = numeric(Nsim),
                      ch.estimate = numeric(Nsim), ch.se = numeric(Nsim), ch.sesq = numeric(Nsim), ch.pval = numeric(Nsim), 
                      total = numeric(Nsim), comb52 = numeric(Nsim), comb26 = numeric(Nsim), covariance = numeric(Nsim), total_se = numeric(Nsim), total_sesq = numeric(Nsim),
                      total_pval= numeric(Nsim), estimate = numeric(Nsim), se = numeric(Nsim), sesq = numeric(Nsim), pval = numeric(Nsim)) 

# Though there are 144 scenarios, the data generation only needs to be run for 18 scenarios
# The full set of 144 scenarios can be manipulated from these
if(k<=12){
  kk <- 1
} else if(k>12 & k<=24){
  kk <- 2
} else if(k %in% c(25:32)){
  kk <- 3
} else if(k %in% c(33:36)){ 
  kk <- 4
} else if(k %in% c(37:44)){
  kk <- 5
} else if(k %in% c(45:48)){  
  kk <- 6
} else if(k %in% c(49:56)){
  kk <- 7
} else if(k %in% c(57:60)){  
  kk <- 8
} else if(k %in% c(61:68)){
  kk <- 9
} else if(k %in% c(69:72)){  
  kk <- 10
} else if(k>72 & k<=84){
  kk <- 11
} else if(k>84 & k<=96){  
  kk <- 12
} else if(k %in% c(97:104)){
  kk <- 13
} else if(k %in% c(105:108)){  
  kk <- 14
} else if(k %in% c(109:116)){
  kk <- 15
} else if(k %in% c(117:120)){
  kk <- 16
} else if(k>120 & k<=132){
  kk <- 17
} else if(k>132 & k<=144){  
  kk <- 18
}

mergedf <- readRDS(paste(paste("rdatafile",kk,sep=""),".Rds",sep=""))

for (i in 1:Nsim){
  
  # cc, 2 tp's
  # cohort = 1 is original
  if(scenarios$design[k]=="cc" & scenarios$tp[k]==2){
    dff <- filter(mergedf[[i]], cohort==1 & (timept==0 | timept ==78))
    cat("cc, 2 tp \n")} 
  
  # cc, 3 tp's
  if(scenarios$design[k]=="cc" & scenarios$tp[k]==3){
    dff <- filter(mergedf[[i]], cohort==1 & (timept==0 | timept ==78 | timept==26)) 
    cat("cc, 3 tp \n")}
  
  # cc, 5 tp's
  if(scenarios$design[k]=="cc" & scenarios$tp[k]==5){
    dff <- filter(mergedf[[i]], cohort==1 & (timept==0 | timept ==78 | timept==26 | timept==13 | timept==52 )) 
    cat("cc, 5 tp \n")}
  
  # cc, cts tp's
  if(scenarios$design[k]=="cc" & scenarios$tp[k]==20){
    dff <- filter(mergedf[[i]], cohort==1 & (timept %in% c(0,3,6,9,13,16,19,23,26,29,32,35,38,41,44,52,58,65,71,78)))
    cat("cc, 20 tp \n")}
  
  # oc, 2 tp's
  if(scenarios$design[k]=="oc" & scenarios$tp[k]==2){
    dff <- filter(mergedf[[i]], timept==0 | timept == 78)
    cat("oc, 2 tp \n")}
  
  # oc, 3 tp's
  if(scenarios$design[k]=="oc" & scenarios$tp[k]==3){
    dff <- filter(mergedf[[i]], timept==0 | timept == 78 | timept==26)
    cat("oc, 3 tp \n")}
  
  # oc, 5 tp's
  if(scenarios$design[k]=="oc" & scenarios$tp[k]==5){
    dff <- filter(mergedf[[i]], timept==0 | timept == 78 | timept==26 | timept==13 | timept==52)
    cat("oc, 5 tp \n")}
  
  # oc, cts tp's
  if(scenarios$design[k]=="oc" & scenarios$tp[k]==20){
    dff <- filter(mergedf[[i]], timept %in% c(0,3,6,9,13,16,19,23,26,29,32,35,38,41,44,52,58,65,71,78))
    cat("oc, 20 tp \n")}
  
  # cs, 2 tp's
  if(scenarios$design[k]=="cs" & scenarios$tp[k]==2){
    dff <- filter(mergedf[[i]], timept==0 | timept == 78)
    dff$subject <- paste(dff$subject, dff$timept, sep=".") # Change repeated measures id's for the cross-sectional analyses
    cat("cs, 2 tp \n")}
  
  # cs, 3 tp's
  if(scenarios$design[k]=="cs" & scenarios$tp[k]==3){
    dff <- filter(mergedf[[i]], timept==0 | timept == 78 | timept == 26)
    dff$subject <- paste(dff$subject, dff$timept, sep=".") # Change repeated measures id's for the cross-sectional analyses
    cat("cs, 3 tp \n")}
  
  # cs, 5 tp's
  if(scenarios$design[k]=="cs" & scenarios$tp[k]==5){
    dff <- filter(mergedf[[i]], timept==0 | timept == 78 | timept==26 | timept==13 | timept==52)
    dff$subject <- paste(dff$subject, dff$timept, sep=".") # Change repeated measures id's for the cross-sectional analyses
    cat("cs, 5 tp \n")}
  
  # cs, 20 tp's
  if(scenarios$design[k]=="cs" & scenarios$tp[k]==20){
    dff <- filter(mergedf[[i]], timept %in% c(0,3,6,9,13,16,19,23,26,29,32,35,38,41,44,52,58,65,71,78))
    dff$subject <- paste(dff$subject, dff$timept, sep=".") # Change repeated measures id's for the cross-sectional analyses
    cat("cs, 20 tp \n")}
  
  ##### First do Kasza 1 timescale analysis model ####
  if(scenarios$design[k]!="cs"){
    # If the design is NOT cross-sectional we can include a subject RE
    
    tryCatch(
      {
        
        #dff$timept <- as.factor(dff$timept)   # Time is treated as continuous for OC estimand so don't run this line
        dff$trt<-as.factor(dff$trt)
        mod <- lmer(y ~ timept + trt:timept + (1|cluster) + (1|cluster:timept) + (1|cluster:subject), dff,  
                    control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
      },
      error=function(err){
        message('On iteration ',i, ' there was an error: ', err)
        flag <<-c(flag,i)
      },
      message=function(m){
        message('On iteration ',i, ' there was a message: ', m, ' ignoring this result')
        msgflag <<-c(msgflag,i)
        results[i,17:19] <<- NA   
      }
      
    )
    
  } else{
    
    tryCatch(
      {
        # If the design IS cross-sectional, drop subject RE
        
        #dff$timept <- as.factor(dff$timept)  # Time is treated as continuous for OC estimand so don't run this line
        dff$trt<-as.factor(dff$trt)
        mod <- lmer(y ~ timept + trt:timept + (1|cluster) + (1|cluster:timept), dff,
                    control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
      },
      error=function(err){
        message('On iteration ',i, ' there was an error: ', err)
        flag <<-c(flag,i)
      },
      message=function(m){
        message('On iteration ',i, ' there was a message: ', m, 'ignoring this result')
        msgflag <<-c(msgflag,i)
        results[i,17:19] <<- NA       
      }
      
    )  
    
  }
  
  #### Now run Kasza 2 timescale extended analysis model ####
  
  if(scenarios$design[k]=="cs" ){ 
    # If the design IS cross-sectional, drop subject RE
    
    tryCatch(
      {
        
        #dff$timept <- as.factor(dff$timept)  # Time is treated as continuous for OC estimand so don't run this line
        dff$trt<-as.factor(dff$trt)
        #dff$ind.time <-as.factor(dff$ord-1)  # Time is treated as continuous for OC estimand so don't run this line
        ind.time <- dff$ord - 1
        mod_lem <- lmer(y ~ timept + ind.time + trt:ind.time + trt:timept + (1|cluster) + (1|cluster:timept), dff,
                        control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5),check.rankX="silent.drop.cols"))
      },
      error=function(err){
        message('On iteration ',i, ' there was an error: ', err)
        flag2 <<-c(flag2,i)
      },
      message=function(m){
        message('On iteration ',i, ' there was a message: ', m, 'ignoring this result')
        msgflag2 <<-c(msgflag2,i)
        results[i,2:16] <<- NA            
      }
      
    )   
    
  } else if(scenarios$design[k]=="oc" ){ 
    
    tryCatch(
      {
        # If the design is NOT cross-sectional we can include a subject RE
        
        #dff$timept <- as.factor(dff$timept)  # Time is treated as continuous for OC estimand so don't run this line
        dff$trt<-as.factor(dff$trt)
        #dff$ind.time <-as.factor(dff$ord-1)  # Time is treated as continuous for OC estimand so don't run this line
        ind.time <- dff$ord - 1
        mod_lem <- lmer(y ~ timept + ind.time + trt:ind.time + trt:timept + (1|cluster) + (1|cluster:timept) + (1|cluster:subject) , dff,  
                        control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5),check.rankX="silent.drop.cols"))
      },
      error=function(err){
        message('On iteration ',i, ' there was an error: ', err)
        flag2 <<-c(flag2,i)
      },
      message=function(m){
        message('On iteration ',i, ' there was a message: ', m, 'ignoring this result')
        msgflag2 <<-c(msgflag2,i)
        results[i,2:16] <<- NA             
      }
      
    )  
    
  }
  
  # Extract elements for 1 timescale model
  if(!is.na(results$estimate[i])){
    sum<-summary(mod)
    mat <- data.frame(coef.name=dimnames(coef(sum))[[1]], coef.value = coef(sum)[,1], coef.se = coef(sum)[,2], 
                      coef.pval = coef(sum)[,5])
    results$estimate[i] <- mat$coef.value[mat$coef.name=="timept:trt1"] 
    results$se[i] <- mat$coef.se[mat$coef.name=="timept:trt1"] 
    results$sesq[i] <- (results$se[i])^2
    results$pval[i] <- mat$coef.pval[mat$coef.name=="timept:trt1"] 
  }
  
  # Extract elements for 2 timescale model - don't attempt for CC as doesn't have 2 timescales
  if((scenarios$design[k]=="oc" ) | (scenarios$design[k]=="cs")){ 
    if(!is.na(results$ch.estimate[i])){
      sum_lem<-summary(mod_lem)
      mat_lem <- data.frame(coef.name=dimnames(coef(sum_lem))[[1]], coef.value = coef(sum_lem)[,1],
                            coef.se = coef(sum_lem)[,2], coef.pval = coef(sum_lem)[,5])
      results$ind.estimate[i] <- ifelse(length(mat_lem$coef.value[mat_lem$coef.name=="ind.time:trt1"])!=0, mat_lem$coef.value[mat_lem$coef.name=="ind.time:trt1"],NA)
      results$ind.se[i] <- ifelse(length(mat_lem$coef.se[mat_lem$coef.name=="ind.time:trt1"])!=0, mat_lem$coef.se[mat_lem$coef.name=="ind.time:trt1"],NA) 
      results$ind.sesq[i] <- (results$ind.se[i])^2
      results$ind.pval[i] <- ifelse(length(mat_lem$coef.pval[mat_lem$coef.name=="ind.time:trt1"])!=0, mat_lem$coef.pval[mat_lem$coef.name=="ind.time:trt1"],NA) 
      results$ch.estimate[i] <- ifelse(length(mat_lem$coef.value[mat_lem$coef.name=="timept:trt1"])!=0, mat_lem$coef.value[mat_lem$coef.name=="timept:trt1"],NA) 
      results$ch.se[i] <- ifelse(length(mat_lem$coef.se[mat_lem$coef.name=="timept:trt1"])!=0, mat_lem$coef.se[mat_lem$coef.name=="timept:trt1"],NA) 
      results$ch.sesq[i] <- (results$ch.se[i])^2
      results$ch.pval[i] <- ifelse(length(mat_lem$coef.pval[mat_lem$coef.name=="timept:trt1"])!=0, mat_lem$coef.pval[mat_lem$coef.name=="timept:trt1"],NA) 
      results$total[i] <- sum(results$ind.estimate[i],results$ch.estimate[i],na.rm=TRUE)
      results$comb52[i] <- 78*results$ch.estimate[i] + 52*results$ind.estimate[i]
      results$comb26[i] <- 78*results$ch.estimate[i] + 26*results$ind.estimate[i]
      results$covariance[i] <- tryCatch({vcov(sum_lem)["timept:trt1","ind.time:trt1"]; vcov(sum_lem)["timept:trt1","ind.time:trt1"]}, error=function(e) NA)
      results$total_se[i] <- ifelse(is.na(results$ind.se[i]) | is.na(results$ch.se[i]) | is.na(results$covariance[i]), 
                                    sum(results$ind.se[i], results$ch.se[i], na.rm=TRUE), sqrt(results$ind.sesq[i] + results$ch.sesq[i] + 2*results$covariance[i]))
      results$total_sesq[i] <- (results$total_se[i])^2
      results$total_pval[i] <- 0 
    }
  }
  
  cat("i is")
  print(i)
}

singul <- length(msgflag)

########################
# 1 timescale model
########################

# Average estimate
finalresults$avg_estimate <- mean(results$estimate[!is.na(results$estimate)])

# Power (not used)
finalresults$power <- sum(results$pval[!is.na(results$estimate)]<0.05)/(Nsim-singul) 

# Average of model-based SE's for our parameter of interest in each simulation (not used - two different formulations)
finalresults$avg_se <- mean(results$se[!is.na(results$estimate)]) 
finalresults$avg_se2 <- sqrt(mean(results$sesq[!is.na(results$estimate)]))

# Empirical SE- equivalent to calculating sd(results$est)
finalresults$emp_se <- sqrt(1/(Nsim-singul-1)*sum((results$estimate[!is.na(results$estimate)]-finalresults$avg_estimate)^2))

# MC SE's for performance measures
finalresults$mc_se_bias <- finalresults$emp_se/sqrt(Nsim-singul)
finalresults$mc_se_power <- sqrt(finalresults$power*(1-finalresults$power)/(Nsim-singul))
finalresults$mc_se_empse <-  finalresults$emp_se/sqrt(2*(Nsim-singul-1))

# Errors
finalresults$error_count <- length(flag)
finalresults$msg_count <- length(msgflag)

# Record where ch.estimate is NA, these will be discounted in the denominator 
ch_na <- sum(is.na(results$ch.estimate))
print(ch_na)

########################
# 2 timescale model
########################

if((scenarios$design[k]=="oc") | (scenarios$design[k]=="cs" )){ 
  
# Average estimates  
finalresults$avg_ind_estimate <- mean(results$ind.estimate[!is.na(results$ch.estimate)])  # The is.na is for the ch.estimate - ONLY FOR OC ESTIMAND
finalresults$avg_ch_estimate <- mean(results$ch.estimate[!is.na(results$ch.estimate)]) # The is.na is for the ch.estimate - ONLY FOR OC ESTIMAND
finalresults$avg_total_estimate <-  mean(results$total[!is.na(results$ch.estimate)]) # Not interested in total for this estimand so can change this 

# Power (not used)
finalresults$ind_power <- sum(results$ind.pval[!is.na(results$ch.estimate)]<0.05)/(Nsim-ch_na)
finalresults$ch_power <- sum(results$ch.pval[!is.na(results$ch.estimate)]<0.05)/(Nsim-ch_na)

# Average of model based SE's (not used)
finalresults$avg_ind_se <- mean(results$ind.se[!is.na(results$ch.estimate)])
finalresults$avg_ind_se2 <- sqrt(mean(results$ind.sesq[!is.na(results$ch.estimate)]))
finalresults$avg_ch_se <- mean(results$ch.se[!is.na(results$ch.estimate)])
finalresults$avg_ch_se2 <- sqrt(mean(results$ch.sesq[!is.na(results$ch.estimate)]))
finalresults$avg_total_se <-  mean(results$total_se[!is.na(results$ch.estimate)])
finalresults$avg_total_se2 <-  sqrt(mean(results$total_sesq[!is.na(results$ch.estimate)]))

# Empirical SE's
finalresults$ind_emp_se <- sqrt(1/(Nsim-ch_na-1)*sum((results$ind.estimate[!is.na(results$ch.estimate)]-finalresults$avg_ind_estimate)^2))
finalresults$ch_emp_se <- sqrt(1/(Nsim-ch_na-1)*sum((results$ch.estimate[!is.na(results$ch.estimate)]-finalresults$avg_ch_estimate)^2))
finalresults$total_emp_se <- sqrt(1/(Nsim-ch_na-1)*sum((results$total[!is.na(results$ch.estimate)]-finalresults$avg_total_estimate)^2)) 
finalresults$oc52emp_se <- sd(results$comb52,na.rm=TRUE)
finalresults$oc26emp_se <- sd(results$comb26,na.rm=TRUE)

# MC SE's for performance measures
finalresults$ind_mc_se_bias <- finalresults$ind_emp_se/sqrt(Nsim-ch_na)
finalresults$ch_mc_se_bias <- finalresults$ch_emp_se/sqrt(Nsim-ch_na)
finalresults$total_mc_se_bias <- finalresults$total_emp_se/sqrt(Nsim-ch_na)

finalresults$ind_mc_se_power <- sqrt(finalresults$ind_power*(1-finalresults$ind_power)/(Nsim-ch_na))
finalresults$ch_mc_se_power <- sqrt(finalresults$ch_power*(1-finalresults$ch_power)/(Nsim-ch_na))
finalresults$ind_mc_se_empse <- finalresults$ind_emp_se/sqrt(2*(Nsim-ch_na-1))
finalresults$ch_mc_se_empse <-  finalresults$ch_emp_se/sqrt(2*(Nsim-ch_na-1))
finalresults$total_mc_se_empse <- finalresults$total_emp_se/sqrt(2*(Nsim-ch_na-1))
finalresults$oc52_mc_se_empse <- finalresults$oc52emp_se/sqrt(2*(Nsim-ch_na-1))	
finalresults$oc26_mc_se_empse <- finalresults$oc26emp_se/sqrt(2*(Nsim-ch_na-1))
finalresults$oc52_mc_se_bias <- finalresults$oc52emp_se/sqrt(Nsim-ch_na)	
finalresults$oc26_mc_se_bias <- finalresults$oc26emp_se/sqrt(Nsim-ch_na)

# Count number of errors and non-convergences
finalresults$error_count2 <- length(flag2)
finalresults$msg_count2 <- length(msgflag2)
} 

#############################################################
# Save each line of finalresults for a single value of k
#############################################################
saveRDS(finalresults, file=paste("ocfinalresults",k,".Rds",sep=""))
