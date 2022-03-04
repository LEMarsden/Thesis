library(foreign); library(haven)

# Create final results data frame which will have 72 rows for the 72 scenarios of Chapter 6
finalresults <- data.frame(scenario = numeric(1), avg_ind_estimate = numeric(1),
                           avg_ind_se = numeric(1), avg_ind_se2 = numeric(1), ind_power = numeric(1),
                           avg_ch_estimate = numeric(1),
                           avg_ch_se = numeric(1), avg_ch_se2 = numeric(1), ch_power = numeric(1),
                           avg_total_estimate = numeric(1), 
                           ind_emp_se = numeric(1), ch_emp_se = numeric(1), total_emp_se = numeric(1),
                           oc52emp_se = numeric(1), oc26emp_se = numeric(1),
                           ind_mc_se_bias = numeric(1), ch_mc_se_bias = numeric(1), total_mc_se_bias = numeric(1),
                           error_count2 = numeric(1),msg_count2 = numeric(1),
                           ind_mc_se_power = numeric(1), ch_mc_se_power = numeric(1),
                           ind_mc_se_empse = numeric(1), ch_mc_se_empse = numeric(1), total_mc_se_empse = numeric(1),
                           oc52_mc_se_empse = numeric(1), oc26_mc_se_empse = numeric(1), 
                           avg_estimate = numeric(1), power = numeric(1),
                           avg_se = numeric(1), avg_se2 = numeric(1), emp_se = numeric(1),
                           mc_se_bias = numeric(1), error_count = numeric(1),msg_count = numeric(1),
                           mc_se_power = numeric(1), mc_se_empse = numeric(1),
                           mc_se_conv = numeric(1), mc_se_conv2 = numeric(1))

for (i in 1:72){ 
file <- paste(paste("statapostfile",i,sep=""),".dta",sep="") 
results <- read_dta(file)

Nsim <- 650

# Manipulation of the results dataframe
results$sesq <- (results$se)^2
results$ind_sesq <- (results$ind_se)^2
results$ch_sesq <- (results$ch_se)^2
results$total <- results$ch_est + results$ind_est 
results$comb52 <- 78*results$ch_est + 52*results$ind_est
results$comb26 <- 78*results$ch_est + 26*results$ind_est
results$total_pval <- 0 # Not used
results$total_sesq <- (results$total_se)^2

# Count number of errors and non-convergences
finalresults$scenario <- i
finalresults$error_count <- length(which(results$err!=0))/650*100 # Err is for 1 timescale model
finalresults$error_count2 <- length(which(results$err2!=0))/650*100 # Err2 is for 2 timescale model
finalresults$nonconv <- length(which(results$conv==0))/650*100   # 0 is given in STATA if not converged
finalresults$nonconv2 <- length(which(results$conv2==0))/650*100 # 0 is given in STATA if not converged
singul <- length(which(results$err!=0)) + length(which(results$conv==0)) # Conv is for 1 timescale model
singul2 <- length(which(results$err2!=0)) + length(which(results$conv2==0)) # Conv2 is for 1 timescale model
finalresults$conv <- singul/650
finalresults$conv2 <- singul2/650


########################
# 1 timescale model
########################

# Average estimate
finalresults$avg_estimate <- mean(results$est[!is.na(results$est)])

# Power (not used)
finalresults$power <- sum(results$pval[!is.na(results$est)]<0.05)/(Nsim-singul)

# Average of model-based SE's for our parameter of interest in each simulation (not used - two different formulations)
finalresults$avg_se <- mean(results$se[!is.na(results$est)]) 
finalresults$avg_se2 <- sqrt(mean(results$sesq[!is.na(results$est)])) 

# Empirical SE- equivalent to calculating sd(results$est)
finalresults$emp_se <- sqrt(1/(Nsim-singul-1)*sum((results$est[!is.na(results$est)]-finalresults$avg_estimate)^2))  

# MC SE's for performance measures
finalresults$mc_se_bias <- finalresults$emp_se/sqrt(Nsim-singul)
finalresults$mc_se_power <- sqrt(finalresults$power*(1-finalresults$power)/(Nsim-singul))
finalresults$mc_se_empse <-  finalresults$emp_se/sqrt(2*(Nsim-singul-1))
finalresults$mc_se_conv <- sqrt(finalresults$conv*(1-finalresults$conv)/Nsim)

########################
# 2 timescale model
########################

# Average estimates
finalresults$avg_ind_estimate <- mean(results$ind_est[!is.na(results$ch_est)])  
finalresults$avg_ch_estimate <- mean(results$ch_est[!is.na(results$ch_est)]) 
finalresults$avg_total_estimate <-  mean(results$total[!is.na(results$ch_est)]) 

# Power (not used)
finalresults$ind_power <- sum(results$ind_pval[!is.na(results$ch_est)]<0.05)/(Nsim-singul2)
finalresults$ch_power <- sum(results$ch_pval[!is.na(results$ch_est)]<0.05)/(Nsim-singul2)

# Average of model based SE's (not used)
finalresults$avg_ind_se <- mean(results$ind_se[!is.na(results$ch_est)])
finalresults$avg_ind_se2 <- sqrt(mean(results$ind_sesq[!is.na(results$ch_est)]))
finalresults$avg_ch_se <- mean(results$ch_se[!is.na(results$ch_est)])
finalresults$avg_ch_se2 <- sqrt(mean(results$ch_sesq[!is.na(results$ch_est)]))
finalresults$avg_total_se <-  mean(results$total_se[!is.na(results$ch_est)])
finalresults$avg_total_se2 <-  sqrt(mean(results$total_sesq[!is.na(results$ch_est)]))

# Empirical SE's
finalresults$ind_emp_se <- sqrt(1/(Nsim-singul2-1)*sum((results$ind_est[!is.na(results$ch_est)]-finalresults$avg_ind_estimate)^2))
finalresults$ch_emp_se <- sqrt(1/(Nsim-singul2-1)*sum((results$ch_est[!is.na(results$ch_est)]-finalresults$avg_ch_estimate)^2))
finalresults$total_emp_se <- sqrt(1/(Nsim-singul2-1)*sum((results$total[!is.na(results$ch_est)]-finalresults$avg_total_estimate)^2)) 
finalresults$oc52emp_se <- sd(results$comb52, na.rm=TRUE)
finalresults$oc26emp_se <- sd(results$comb26, na.rm=TRUE)

# MC SE's for performance measures
finalresults$ind_mc_se_bias <- finalresults$ind_emp_se/sqrt(Nsim-singul2)
finalresults$ch_mc_se_bias <- finalresults$ch_emp_se/sqrt(Nsim-singul2)
finalresults$total_mc_se_bias <- finalresults$total_emp_se/sqrt(Nsim-singul2)
finalresults$ind_mc_se_power <- sqrt(finalresults$ind_power*(1-finalresults$ind_power)/(Nsim-singul2))
finalresults$ch_mc_se_power <- sqrt(finalresults$ch_power*(1-finalresults$ch_power)/(Nsim-singul2))
finalresults$ind_mc_se_empse <- finalresults$ind_emp_se/sqrt(2*(Nsim-singul2-1))
finalresults$ch_mc_se_empse <-  finalresults$ch_emp_se/sqrt(2*(Nsim-singul2-1))
finalresults$total_mc_se_empse <- finalresults$total_emp_se/sqrt(2*(Nsim-singul2-1))
finalresults$oc52_mc_se_empse <- finalresults$oc52emp_se/sqrt(2*(Nsim-singul2-1))	
finalresults$oc26_mc_se_empse <- finalresults$oc26emp_se/sqrt(2*(Nsim-singul2-1)) 
finalresults$oc52_mc_se_bias <- finalresults$oc52emp_se/sqrt(Nsim-singul2)
finalresults$oc26_mc_se_bias <- finalresults$oc26emp_se/sqrt(Nsim-singul2)
finalresults$mc_se_conv2 <- sqrt(finalresults$conv2*(1-finalresults$conv2)/Nsim)

# Emp SE's - multiply by 78 - for 1 ts
finalresults$emp_se78 <- finalresults$emp_se*78

# Emp SE's - multiply by 78 - for 2 ts
finalresults$emp_se78_2 <- finalresults$total_emp_se*78

# MC SE's for bias - multiply by 78 - for 1 ts
finalresults$mc_se_bias78 <- finalresults$mc_se_bias*78

# MC SE's for bias - multiply by 78 - for 2 ts
finalresults$total_mc_se_bias78 <- finalresults$total_mc_se_bias*78

# MC SE's for empSE - multiply by 78 - for 1 ts
finalresults$mc_se_empse78 <- finalresults$mc_se_empse*78

# MC SE's for empSE - multiply by 78 - for 2 ts
finalresults$total_mc_se_empse78 <- finalresults$total_mc_se_empse*78


saveRDS(finalresults, file=paste("jointfinalresults",i,".Rds",sep=""))
print(i)
}