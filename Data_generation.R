library(dplyr)
library(ggplot2)
library(MixMatrix)
library(mvtnorm)

################################################################################
# SIMULATE A UNIVERSAL DATASET UNDER A VARIETY OF SCENARIOS 
################################################################################

univ <- function(no_clusters, clust_size_m, clust_pop_k, md, maxfup, ind_coef, 
                 ch_coef, icc, cac, iac, sigma, beta1, beta2, beta3, design, 
                 tp, nonlin)
  
  # no_clusters = number of clusters; 
  # clust_size_m = cluster sample size;
  # clust_pop_k = cluster population size;
  # md = missing data mechanism;
  # maxfup = maximum length of follow-up;
  # ind_coef = coefficient for ind.time by trt interaction;
  # ch_coef = coefficient for cr.time by trt interaction;
  # icc = intraclass correlation coefficient;
  # cac = cluster autocorrelation;
  # iac = individual autocorrelation;
  # sigma = 
  # beta1, beta2, beta3 = coefficients for tuning the missing data mechanism;
  # design = design (oc, cs or cc);
  # tp = number of time points (2, 3, 5 or cts);
  # nonlin = 1 if the cluster intervention effect rate is constant, 2 for non-constant
  
{

  #####################################################
  # GENERAL - DOESN'T DEPEND ON MISSING DATA MECHANISM
  #####################################################
  
  if(clust_size_m==clust_pop_k){
    sample <- "full"
  } else{
      sample <- "sub"
  }
  
  # If full-sample, don't need to do anything else, loops run as normal 
  # If sub-sample, do other stuff later to sample from population depending on the design
  if(sample=="full"){
    cluster_size <- clust_size_m 
  } else if(sample=="sub"){
      cluster_size <- clust_pop_k
  }  
  
  # Total number of individuals     
  N <- no_clusters*cluster_size
  
  #######################################################
  # FUNCTION TO GENERATE THE LONGITUDINAL MEASUREMENTS
  #######################################################
  
  # Note that this fn has cluster_size as a general term which can be either clust_size_m or clust_pop_k depending on whether full/sub-sample
    simlong <- function(alltimepts, N, no_clusters, cluster_size, icc, cac, iac, sigma, ind_coef, ch_coef)
    {
      # Variables needed for dataframe including time
      timept <- rep(alltimepts,N)
      obs_pp <- length(alltimepts) 
      ord <- rep(1:obs_pp,N)
      obs <- seq(1,N*obs_pp)
      
      # Create cluster and subject numbers - control first half, trt second
      # clust_size_m for full, clust_pop_k for sub
      cluster <- rep(1:no_clusters, each=obs_pp*cluster_size)  
      subject <- rep(1:N,each=obs_pp)
      trt <- c(rep(0,length(obs)/2),rep(1,length(obs)/2))
      
      # CALCULATE VARIANCE COMPONENTS 
      # CLUSTER LEVEL
      var_chi <- icc*sigma^2
      
      # SUBJECT LEVEL
      var_eps <- sigma^2 - var_chi
      
      # CLUSTER-SPECIFIC RANDOM EFFECTS VECTORS - ONE OF THESE FOR EACH CLUSTER 
      # Times going across, clusters going down 
      cs_clust <- CSgenerate(obs_pp,cac)
      chi <- rmvnorm(n=no_clusters, mean=rep(0,times=obs_pp),sigma= var_chi*cs_clust)
      
      # SUBJECT-SPECIFIC RANDOM EFFECTS VECTORS - ONE OF THESE FOR EACH INDIVIDUAL
      # Times going across, individuals going down
      cs_indiv <- CSgenerate(obs_pp,iac)
      eps <- rmvnorm(n=N, mean=rep(0,times=obs_pp),sigma= var_eps*cs_indiv) 

      # RANDOM CLUSTER EFFECT FOR DROP-OUT MODEL, CALL IT B
      # clust_size_m for full, clust_pop_k for sub 
      b <- rep(rnorm(n=no_clusters, mean=0, sd=1),each=obs_pp*cluster_size) 
      
      # Turn the multivariate matrices into one long vector 
      # clust_size_m for full, clust_pop_k for sub 
      rows <- rep(1:no_clusters,each=cluster_size) 
      chi_long <- chi[rows,]
      clust_re=as.vector(t(chi_long))
      indiv_re=as.vector(t(eps))
      
      # Sum the individual and cluster intervention effect rates
      total_coef <- ind_coef + ch_coef
      
      # Create the longitudinal measurements, y
      if(nonlin==1){
        y <- 0.01*timept + total_coef*trt*timept + clust_re + indiv_re
      } else if(nonlin==2 & ch_coef==-0.003855){
        y <- 0.01*timept + ind_coef*trt*timept + (timept<8)*0*trt + (timept>=8 & timept <78)*(0.30095*exp(-0.1*(timept-8))-0.30095)*trt + (timept==78)*(-0.30069)*trt + clust_re + indiv_re
      } else if(nonlin==2 & ch_coef==-0.001285){
        y <- 0.01*timept + ind_coef*trt*timept + (timept<8)*0*trt + (timept>=8 & timept <78)*(0.10032*exp(-0.1*(timept-8))-0.10032)*trt + (timept==78)*(-0.10023)*trt + clust_re + indiv_re
      }
      
      # Create dataframe
      longdf <- data.frame(cluster, trt, timept, ord, subject, clust_re,indiv_re, y, b)
      return(longdf)
    }
    
    #######################################################
    # FUNCTION TO GENERATE THE TIMES TO EVENT
    #######################################################

    rexp = function(beta,x,time.axis,b) 
    {
      # Number of time points
      ni = length(time.axis) 
      Delta.t = time.axis[2:ni] - time.axis[1:(ni-1)]
      
      # This is because for MCAR we have a constant only and can't do matrix multiplication
      if (is.null(ncol(x))){
        exbeta = exp(c(x*beta + b))
      } else{
          exbeta = exp(c(x%*%beta + b)) 
      }

      # Simulate random uniform variable
      u = runif(1)
      ls = -log(u)

      # Find the appropriate intervals as in inverse CDF theorem 
      # Time intervals * exp(x*beta)
      cut = c(0,cumsum(Delta.t*exbeta[1:(ni-1)])) 
      
      # Find which of the above intervals the ls variable is in
      int = findInterval(ls,cut) 
      
      # cut[int] is the lower bound of the interval ls is in. Find how far into 
      # the interval it is, and divide by the value of exbeta in this interval 
      # (see formula). Then add on the lower bound timepoint of the interval
      sampl =  (ls - cut[int])/exbeta[int] + time.axis[int]
      
      return(sampl)
    }

  ##############################################################
  # MISSING DATA MECHANISMS FOR ORIGINAL COHORT
  ##############################################################

    if(md=="MAR")
    {
    beta_w <- c(beta1,beta2)
    alltimepts <- seq(0,maxfup,by=1) 
    alltimepts <- alltimepts 
    mardf <- simlong(alltimepts, N, no_clusters, cluster_size, icc, cac, iac, sigma, ind_coef, ch_coef)
    clust_re <- mardf$clust_re
    
    obs_pp <- length(alltimepts)
    
    # Vector of times
    time.axis_w <- mardf$timept[mardf$subject==1] 

    surv = rep(NA,N)

    for (i in 1:N) 
    {
      # Survival times for each individual
      surv[i] = rexp(beta_w,x = cbind(1,mardf$y[mardf$subject==i]),time.axis_w, 
                     mardf$b[mardf$subject==i]) 
    }
    
    # Dataframe with the first obs of every subject
    data.id = mardf[mardf$ord==1,]

    # The observed survival time is the minimum of true survival time and maxfup
    data.id$time = pmin(max(alltimepts),surv) 
    data.id$event = 1*(surv<max(alltimepts))

    # Create missing indicator
    mardf$mis = F

    # Now delete the appropriate cases from the longitudinal dataframe

    for (i in which(surv<max(alltimepts)))
    {
      # For each patient who failed, set observations after the surv time as missing
      mardf$mis[mardf$subject==i]=mardf$timept[mardf$subject==i]>surv[i]
    }

    # Observed data
    data.obs = mardf[mardf$mis==FALSE,]
    data.id[data.id$time==max(alltimepts),"time"] = max(alltimepts)+0.01

    data.obs$time = data.id$time[data.obs$subject]
    data.obs$event = data.id$event[data.obs$subject]
    }

   else if(md=="MNAR")
    {
    beta_w <- c(beta1,beta2,beta3)
    alltimepts <- seq(0,maxfup+1,by=1)

    alltimepts <- alltimepts
    mnardf <- simlong(alltimepts, N, no_clusters, cluster_size, icc, cac, iac, sigma, ind_coef, ch_coef)
    clust_re <- mnardf$clust_re
    
    obs_pp <- length(alltimepts)
    
    # Vector of times
    time.axis_w <- mnardf$timept[mnardf$subject==1 & mnardf$ord!=obs_pp] 

    surv = rep(NA,N)

    for (i in 1:N)  
    {
      # Survival times for each individual
      surv[i] = rexp(beta_w,x = cbind(1,mnardf$y[mnardf$subject==i & mnardf$ord != obs_pp],
                                      mnardf$y[mnardf$subject==i & mnardf$ord != 1]), time.axis_w,
                                      mnardf$b[mnardf$subject==i & mnardf$ord!=obs_pp]) 
    }

    # Dataframe with the first obs of every subject
    data.id = mnardf[mnardf$ord==1,]

    # The observed survival time is the minimum of true survival time and maxfup
    data.id$time = pmin(max(time.axis_w),surv) 
    data.id$event = 1*(surv<max(time.axis_w))   

    # Create missing indicator
    mnardf$mis = F

    # Now delete the appropriate cases from the longitudinal dataframe
    e.times = rep(surv,each = obs_pp)
    mnardf$mis = mnardf$timept > e.times 
    
    # Always set the missingness indicator to true for the last msmt time
    mnardf$mis[mnardf$ord == obs_pp] = T 

    # Observed data
    data.obs = mnardf[mnardf$mis==FALSE,]
    data.id[data.id$time==max(time.axis_w),"time"] = max(time.axis_w)+0.01

    data.obs$time = data.id$time[data.obs$subject]
    data.obs$event = data.id$event[data.obs$subject]
    }

    else if(md=="MCAR")
    {
    beta_w <- beta1
    alltimepts <- seq(0,maxfup,by=1)
    mcardf <- simlong(alltimepts, N, no_clusters, cluster_size, icc, cac, iac, sigma, ind_coef, ch_coef)
    clust_re <- mcardf$clust_re
    
    obs_pp <- length(alltimepts)
    secarg <- rep(1,length(alltimepts))
 
    surv = rep(NA,N)
    
    for(i in 1:N)
    {
      # Note that x is just a column of ones for the MCAR case because we don't include previous or next values in exponent
      surv[i] = rexp(beta_w,x = secarg,alltimepts, mcardf$b[mcardf$subject==i])
    }

    # Dataframe with the first obs of every subject
    data.id = mcardf[mcardf$ord==1,]

    # The observed survival time is the minimum of true survival time and maxfup
    data.id$time = pmin(max(alltimepts),surv) 
    data.id$event = 1*(surv<max(alltimepts))

    # Create missing indicator
    mcardf$mis = F

    # Now delete the appropriate cases from the longitudinal dataframe
    for (i in which(surv<max(alltimepts)))
    {
      # For each patient who failed, set observations after the surv time as missing
      mcardf$mis[mcardf$subject==i]=mcardf$timept[mcardf$subject==i]>surv[i]
    }

    # Observed data
    data.obs = mcardf[mcardf$mis==FALSE,]
    data.id[data.id$time==max(alltimepts),"time"] = max(alltimepts)+0.01

    data.obs$time = data.id$time[data.obs$subject]
    data.obs$event = data.id$event[data.obs$subject]
    }

    
    ##########################################
    # STEADY STATE ASSUMPTION - REPLACEMENTS
    ##########################################

    # This function will keep adding replacements for subject i until the maxfup time
    ss <- function(df,N)
    {
      # If the last timepoint in the df is equal to maxfup already, stop:
      if(max(df$timept) == maxfup){return(df)} 
      
      # p = 1 is the first replacement etc
      p <- 1 
      
      # fulldf is the df we keep adding the new replacements to as p increases
      fulldf<-df 
      
      while(max(fulldf$timept)<(maxfup))
      {
        # tempdf is a df for the replacement person
        # First look at how many remaining time points there are for this person
        # To do this look at their entry time 
        # Convert the original person's leavetime into days, round to days, round up, convert back to weeks
        entrytmp <- ceiling(max(fulldf$time)*7)/7

        # ceiling(entrytmp) is entry time. entrytmp is in weeks 
        # NOTE: entrytmp is an exact (integer) day but it is in weeks so non-integer in weeks 
        # Find next measurement point after their entry time and extend to maxfup
        if(ceiling(entrytmp)> maxfup){
          break
        } else if(ceiling(entrytmp)== maxfup){  
          tempdf <- data.frame(timept = maxfup)
        } else if(ceiling(entrytmp)< maxfup){
          tempdf <- data.frame(timept = seq(((ceiling(entrytmp))),maxfup,by=1)) 
        }

        tempdf$gaptime <- 0
        
        # The (max possible) number of msmts they will have in vector
        tsize <- length(tempdf$timept) 
        tempdf$entrytime <- rep(entrytmp,tsize)

        # Variables that stay the same for every replacement as same cluster
        tempdf$b <- rep(df$b[1],tsize)
        tempdf$cluster <- rep(df$cluster[1],tsize)
        tempdf$trt <- rep(df$trt[1],tsize)
        tempdf$bed <- rep(df$bed[1],tsize)
        tempdf$cohort <- 0
        
        # Make sure to bring through correct cluster times
        tempdf$clust_re <- clust_re[seq(from=ceiling(entrytmp)+1,length=tsize,by=1)] 
        
        # Calculate these later so leave in correct format for now
        tempdf$mis <- rep('NA',tsize)
        tempdf$time <- rep(99,tsize)
        tempdf$event <- rep(99,tsize)

        # Variables that will change for each person
        tempdf$bed_pers <- rep(max(fulldf$bed_pers)+1,tsize)
        tempdf$subject <- rep(df$subject[1] + p*N,tsize)  
        tempdf$ord <- 1:tsize
        
        # CALCULATE VARIANCE COMPONENTS FOR REPLACEMENTS
        # CLUSTER LEVEL
        var_chi <- icc*sigma^2
        
        # SUBJECT LEVEL
        var_eps <- sigma^2 - var_chi
        
        # Need to do this otherwise if tsize is 1 ie. one entry only, we can't create a CS matrix 
        if(tsize==1){
          eps <- rnorm(n=1,mean=0,sd=sqrt(var_eps))
          tempdf$indiv_re=eps
        } else{
          # SUBJECT-SPECIFIC RANDOM EFFECTS VECTORS - ONE OF THESE FOR EACH INDIVIDUAL
          cs_indiv <- CSgenerate(tsize,iac)
          
          # Times going across, individuals going down
          eps <- rmvnorm(n=1, mean=rep(0,times=tsize),sigma= var_eps*cs_indiv) 
          tempdf$indiv_re=as.vector(t(eps)) 
        }

        # Now for the case where there are 77 measurements in the dataset already and
        # their first measurement is at time=78. Calculate their first y value but 
        # don't do any of the MD stuff as breaks the while loop
        if(ceiling(entrytmp)== maxfup){ 
          tempdf$mis <- FALSE
          tempdf$time <- maxfup+0.01
          tempdf$event <- 0
          
          # In the below, ind.time is ord - 1 and cr.time is timept
          if(nonlin==1){
            tempdf$y <- 0*tempdf$timept + 0.01*(tempdf$ord-1) + ch_coef*tempdf$trt*tempdf$timept + ind_coef*tempdf$trt*(tempdf$ord-1) + tempdf$clust_re + tempdf$indiv_re
          } else if(nonlin==2 & ch_coef==-0.003855){
            tempdf$y <- 0*tempdf$timept + 0.01*(tempdf$ord-1) + ind_coef*tempdf$trt*(tempdf$ord-1) + (tempdf$timept<8)*0*tempdf$trt + 
              (tempdf$timept>=8 & tempdf$timept <78)*(0.30095*exp(-0.1*(tempdf$timept-8))-0.30095)*tempdf$trt + (tempdf$timept>=78)*(-0.30069)*tempdf$trt + tempdf$clust_re + tempdf$indiv_re
          } else if(nonlin==2 & ch_coef==-0.001285){
            tempdf$y <- 0*tempdf$timept + 0.01*(tempdf$ord-1) + ind_coef*tempdf$trt*(tempdf$ord-1) + (tempdf$timept<8)*0*tempdf$trt + 
              (tempdf$timept>=8 & tempdf$timept <78)*(0.10032*exp(-0.1*(tempdf$timept-8))-0.10032)*tempdf$trt + (tempdf$timept>=78)*(-0.10023)*tempdf$trt + tempdf$clust_re + tempdf$indiv_re
          }
          break
        }

        # Calculate y for this replacement person using replacement formula
        if(nonlin==1){
          tempdf$y <- 0*tempdf$timept + 0.01*(tempdf$ord-1) + ch_coef*tempdf$trt*tempdf$timept + ind_coef*tempdf$trt*(tempdf$ord-1) + tempdf$clust_re + tempdf$indiv_re
        } else if(nonlin==2 & ch_coef==-0.003855){
          tempdf$y <- 0*tempdf$timept + 0.01*(tempdf$ord-1) + ind_coef*tempdf$trt*(tempdf$ord-1) + (tempdf$timept<8)*0*tempdf$trt + 
            (tempdf$timept>=8 & tempdf$timept <78)*(0.30095*exp(-0.1*(tempdf$timept-8))-0.30095)*tempdf$trt + (tempdf$timept>=78)*(-0.30069)*tempdf$trt + tempdf$clust_re + tempdf$indiv_re
        } else if(nonlin==2 & ch_coef==-0.001285){
          tempdf$y <- 0*tempdf$timept + 0.01*(tempdf$ord-1) + ind_coef*tempdf$trt*(tempdf$ord-1) + (tempdf$timept<8)*0*tempdf$trt + 
            (tempdf$timept>=8 & tempdf$timept <78)*(0.10032*exp(-0.1*(tempdf$timept-8))-0.10032)*tempdf$trt + (tempdf$timept>=78)*(-0.10023)*tempdf$trt + tempdf$clust_re + tempdf$indiv_re
        }

        #########################################
        # MD MECHANISMS FOR REPLACEMENTS
        #########################################

            if(md=="MCAR")
            {
              beta_w <- beta1 
              
              # Vector of times
              time.axis_w <- tempdf$timept 

              # Note that x is just a column of ones for the MCAR case 
              tempsurv <- rexp(beta_w,x = rep(1,length(time.axis_w)),time.axis_w,tempdf$b )
              
              # Dataframe with the first obs of every subject
              # In this case it's just the first obs for our temp replacement
              data.id = tempdf[tempdf$ord==1,]

              # The observed survival time is the minimum of true survival time and maxfup
              data.id$time = pmin(max(time.axis_w),tempsurv) 
              data.id$event = 1*(tempsurv<max(time.axis_w))

              # Create missing indicator
              tempdf$mis = F

              # Delete the appropriate cases from the longitudinal dataframe
              # Set missing to true if timept > surv time
              tempdf$mis=tempdf$timept>tempsurv

              # Observed data
              data.obs = tempdf[tempdf$mis==FALSE,]
              data.id[data.id$time==max(time.axis_w),"time"] = max(time.axis_w)+0.01

              data.obs$time = data.id$time
              data.obs$event = data.id$event
            }

        else if(md=="MAR")
        {
          beta_w <- c(beta1,beta2)

          # Vector of times
          time.axis_w <- tempdf$timept 

          tempsurv <- rexp(beta_w,x = cbind(1,tempdf$y),time.axis_w,tempdf$b)

          # Dataframe with the first obs of every subject
          # In this case it's just the first obs for our temp replacement
          data.id = tempdf[tempdf$ord==1,]

          # The observed survival time is the minimum of true survival time and maxfup
          data.id$time = pmin(max(time.axis_w),tempsurv) 
          data.id$event = 1*(tempsurv<max(time.axis_w))

          # Create missing indicator
          tempdf$mis = F

          # Delete the appropriate cases from the longitudinal dataframe
          # Set missing to true if timept > surv time
          tempdf$mis=tempdf$timept>tempsurv

          # Observed data
          data.obs = tempdf[tempdf$mis==FALSE,]
          data.id[data.id$time==max(time.axis_w),"time"] = max(time.axis_w)+0.01

          data.obs$time = data.id$time
          data.obs$event = data.id$event
        }

        else if(md=="MNAR")
        {
          beta_w <- c(beta1,beta2,beta3)

          # IMPORTANT - for MNAR need to generate one further msmt row t = 79
          # Take the last row first and adapt indiv_re, timept, ord, clust_re and y only
          j <- filter(tempdf, timept== maxfup)
          j$timept <- maxfup+1
          j$ord <- j$ord + 1
          j$indiv_re <- rnorm(n=1,mean=0,sd=sqrt(var_eps))
          
          # The only way to get the t+1'th clust_re is to generate 80 chi variables upfront and then in the other md scenarios 
          # just don't use it
          if(nonlin==1){
            j$y <- 0*j$timept + 0.01*(j$ord-1) + ch_coef*j$trt*j$timept + ind_coef*j$trt*(j$ord-1) + clust_re[80] + j$indiv_re 
          } else if(nonlin==2 & ch_coef==-0.003855){
            j$y <- 0*j$timept + 0.01*(j$ord-1) + ind_coef*j$trt*(j$ord-1) + (j$timept<8)*0*j$trt + (j$timept>=8 & j$timept <78)*(0.30095*exp(-0.1*(j$timept-8))-0.30095)*j$trt + (j$timept>=78)*(-0.30069)*j$trt + clust_re[80] + j$indiv_re
          } else if(nonlin==2 & ch_coef==-0.001285){
            j$y <- 0*j$timept + 0.01*(j$ord-1) + ind_coef*j$trt*(j$ord-1) + (j$timept<8)*0*j$trt + (j$timept>=8 & j$timept <78)*(0.10032*exp(-0.1*(j$timept-8))-0.10032)*j$trt + (j$timept>=78)*(-0.10023)*j$trt + clust_re[80] + j$indiv_re
          }
          
          tempdf <- rbind(tempdf,j)

          # Vector of times (NOT including 79)
          time.axis_w <- tempdf$timept[tempdf$timept!=(maxfup +1)] 

          tempsurv <- rexp(beta_w,x = cbind(1,tempdf$y[tempdf$timept != (maxfup +1)],tempdf$y[tempdf$ord != 1]),
                           time.axis_w,tempdf$b[tempdf$timept != (maxfup+1)])

          # Dataframe with the first obs of every subject
          # In this case it's just the first obs for our temp replacement
          data.id = tempdf[tempdf$ord==1,]

          # The observed survival time is the minimum of true survival time and maxfup
          data.id$time = pmin(max(time.axis_w),tempsurv) 
          data.id$event = 1*(tempsurv<max(time.axis_w))

          # Create missing indicator
          tempdf$mis = F

          # Delete the appropriate cases from the longitudinal dataframe
          # Set missing to true if timept > surv time
          tempdf$mis=tempdf$timept>tempsurv

          # Always set the missingness indicator to true for the last msmt time
          tempdf$mis[tempdf$timept == (maxfup+1)] = T 

          # Observed data
          data.obs = tempdf[tempdf$mis==FALSE,]
          data.id[data.id$time==max(time.axis_w),"time"] = max(time.axis_w)+0.01

          data.obs$time = data.id$time
          data.obs$event = data.id$event
        }

        # Finally, add these observations from this replacement to the overall df for this bed
        fulldf <- rbind(fulldf,data.obs)
        p <- p + 1
      }

      # The while loop has broken for the below conditions so need to still bind the previous
      # dataframe with the tempdf with one measurement in it
      if(ceiling(entrytmp) == maxfup){
        fulldf <- rbind(fulldf,tempdf)
      }

      return(fulldf)

    }

    # Applies the steady state (ss) function over the whole dataframe
    fn <- function(df)
    {
      # Before running, need to put bed, bed person, entry and leave time in data.obs df
      # If cohort = 1 this is for the original cohort, 0 = additional cohort
      df <- mutate(df, bed_pers = 1, bed = subject, entrytime = 0, gaptime = 0, cohort=1) 
      
      # Pre-allocate a list
      listofdfs <- vector("list", max(df$subject)) 
    

      for (i in 1:length(unique(df$subject)))
      {
        # Create temp df of one specific subject in the original cohort
        tmp <- filter(df,subject== unique(df$subject)[i])  
        
        # Based on this person now fill with replacements
        # And remove any entries where the timept is > maxfup
        listofdfs[[i]] <- filter(ss(tmp,N),timept<=maxfup) 
      }
      
      # Combine all the dataframes (rows of subjects)
      return(do.call("rbind",listofdfs)) 
    }

    ##################################
    # CHOICE OF DESIGN
    ##################################

    # Function to filter by different timepoints
    tps <- function(df, num)
    {
      if (num == 2)
      {
        return(filter(df, timept==0 | as.character(timept) ==maxfup))
      }
      else if (num == 3)
      {
        return(filter(df, timept==0 | as.character(timept) == maxfup/3 | as.character(timept) ==maxfup)) 
      }
      else if (num == 5)
      {
        return(filter(df, as.character(timept)==0 | as.character(timept) == maxfup/6 | as.character(timept) == maxfup/3 | as.character(timept) == maxfup*(2/3) | as.character(timept) ==maxfup))
      }
      else if (num == 20)
      {
        return(filter(df, as.character(timept)==0 | as.character(timept) == 3| as.character(timept) == 6 | as.character(timept) == 9 | as.character(timept) == maxfup/6 
                      | as.character(timept) == 16 | as.character(timept) == 19 | as.character(timept) == 23 | as.character(timept) == maxfup/3
                      | as.character(timept) == 29 | as.character(timept) == 32 | as.character(timept) == 35 | as.character(timept) == 38
                      | as.character(timept) == 41 | as.character(timept) == 44 | as.character(timept) == 52 | as.character(timept) == 58 
                      | as.character(timept) == 65  | as.character(timept) == 71  | as.character(timept) == maxfup))
      }
      else if (num == "cts")
      {
        return(df)
      }
    }

    
    # If taking a full sample, just need to run the tps function for final output
      if(sample=="full"){
      if(design=="cc"){
        z <- tps(data.obs,tp)
        return(z)
      } else{
        z <- tps(fn(data.obs),tp)
        return(z)  
      }}
    
    
    # If taking a subsample, this is more complicated
    # Next code makes sure that the sample at timept = 0 is the same for all three designs
  
    if(sample=="sub"){ 
      
      # For this section if we have MNAR need to tell program that obs_pp isn't 80 any more but 79
      if(md=="MNAR"){obs_pp <- 79} 
      
      # Make array of who to sample at each time point - m
      Nsample <- array(NA,dim=c(no_clusters,clust_size_m,obs_pp)) 
      
      # Array of who is present at each time point - k
      # I have multipled by 1.4 here to accommodate for the largest possible k size we could have 
      who <- array(NA,dim=c(no_clusters,round(clust_pop_k*1.4),obs_pp)) 
      
      # Pre-allocate a list to fill 
      popnframe <- z <- vector("list", obs_pp) 
      
      # Everyone in the popn at timept j 
      popnframe[[1]] <- filter(data.obs,timept==0) 
      
      # Only do this for j = 1 for timept = 0
      for(i in 1:no_clusters){ 
        # List of who is in the population at time point j
        n <- unique(popnframe[[1]]$subject[popnframe[[1]]$cluster==i])
        
        # Who is there
        who[i,1:length(n),1] <- n 
        
        # Take a sample size m out of who is there
        Nsample[i,,1] <- sample(who[i,1:length(n),1], size=clust_size_m, replace=F) 
      } 
      
      # This is all the timept 0 measurements from those chosen at time 0
      z[[1]] <- filter(popnframe[[1]], subject %in% Nsample[,,1]) 
      z[[1]] <- mutate(z[[1]], bed_pers = 1, bed = subject, entrytime = 0, gaptime = 0, cohort=1)
      
      if(design=="cc"){   
        # For CC only need those present at timept 0 but need rest of their values
        z <- tps(filter(data.obs,subject %in% Nsample[,,1]),tp)  
        return(z)
      }
      
      if(design=="oc"){ 
        # The initial sample has already been taken, so just apply fn to it and return that
        z <- tps( fn(filter(data.obs,subject %in% Nsample[,,1])), tp)
        return(z)
      }
      
      if(design=="cs"){
        fullpop <- fn(data.obs) 
        
        for (j in 2:obs_pp){
          for(i in 1:no_clusters){
            # Everyone in the popn at timept j
            popnframe[[j]] <- filter(fullpop,timept==(j-1)) 
            
            # List of who is in the population at time point j
            n <- unique(popnframe[[j]]$subject[popnframe[[j]]$cluster==i]) 
            
            who[i,1:length(n),j] <- n 
            
            # The cs design takes random samples at each tp
            Nsample[i,,j] <- sample(who[i,1:length(n),j], size=clust_size_m, replace=F) 
            z[[j]] <- filter(popnframe[[j]], subject %in% Nsample[,,j])
          }} 
        return(tps(do.call("rbind",z),tp))
        
      } 
    }
}
