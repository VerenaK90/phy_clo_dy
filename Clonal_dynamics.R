
####################  
## Estimate the turnover rate during tumor growth
####################

## Define a mutation rate and a target tumor size; then plot the relationship between cellular death and number of mutations per tumor cell
## Provide the average number of mutations from your data
  
  mut.rate <- 2.66*10^-9*3.3*10^9
  tumor.size <- 10^9
  data.median <- 2300
  data.quantiles <- c(1900, 2900)

  plot(seq(0, 1000000, 1), apply(rbind(0, 1-log(tumor.size)/(seq(0,1000000, 1)/mut.rate )), 2, max), type="l", xlab = "# mutations", ylab = "delta_0", log="x",
       cex = 1.5, cex.lab = 1.5, cex.axis = 1.5)
  lines(seq(0, 1000000, 1), apply(rbind(0, 1-9*log(10)/(seq(0,1000000, 1)/mut.rate/4 )), 2, max), lty=2)
  legend(x = "topleft", lty = c(1,2), legend = c("somatic mutation rate", "4x somatic mutation rate"))
  abline(v = data, col = "red", lty=1)
  polygon(x = c(data.quantiles, rev(data.quantiles)),
          y=c(-0.2, -0.2, 1.2, 1.2), col=rgb(1,0,0, alpha=0.2), border = NA)
  
  ## obtain the turnover rate for your data
  
  1-log(tumor.size)/(data.median/mut.rate )
  1-log(tumor.size)/(data.quantiles/mut.rate )

  
  
####################  
## fit a two stage model with fast turnover followed by TERT promoter mutation
####################
  
  library(rootSolve)
  ## P2 computes the probability for at least 1 TERT mutation after n divisions (1 - (1-mut.rate)^n_div)
  ## The number of divisions, n_div, is 1/(1-delta)*(exp((1-delta)*time)-1), 
  ## where delta is the death rate relative to the proliferation rate, lambda; mut.rate is the single-base-pair substitution rate per division x 2 (2 loci in the TERT promoter); and unif is a random number from a uniform distribution

  ## Function to solve with root-solve: when is the probability of at least 1 mutation equal to unif (unif is a random number)?
  ## in this sense we sample a random timepoint of the TERT promoter mutation by sampling unif and assessing when the probability of at least one TERT promoter muation is equal to unif
  
  P2 <- function(time, delta, unif, mut.rate){
     1-(1-mut.rate)^(1/(1-delta)*(exp((1-delta)*time)-1)) - unif
  }
  

  ## Tumor growth of a tumor consisting of a TERT-WT (delta) and a TERT-mutated population (delat2, starts growing after start2).
  ## Function to solve with root-solve: When has the tumor grown to a tumor size of 10^9 cells?
  tumor.size2 <- function(time, delta, delta2, start2, tumor.size = 10^9){
    exp((1-delta)*time) + exp((1-delta2)*(time - start2))-tumor.size
  }
  
  ## Function to compute the relative fraction of a tumor sample that carries a TERT promoter mutation over time
  ## M: mutated fraction (grows at delta_tert after start_tert (timepoint of mutation acquisition)), S: non-mutated fraction (grows at delta)
  relative.tert.proportion <- function(time, delta, delta_tert, start_tert){
    M <- ceiling(exp((1-delta_tert)*(time - start_tert)))
    S <- ceiling(exp((1-delta)*time))
    ## simulate sampling
    res <- apply(rbind(M,S), 2, function(x){
      rbinom(n=1, size = sum(x), prob=x[1]/sum(x))/sum(x)
    })
    ## no TERT mutated tumor fraction if delta_tert is larger than 1 (death > division) or before the mutation
    res[delta_tert > 1 | start_tert > time | delta_tert < 0] <- 0
    res
  }
  
  
  ## We estimate the error between simulated and measured mean and variance of the TERT mutated tumor fraction. To this end, try different death rates for wt and mutated fractions and record the error at each parameter set
  ## At each parameter set we do 1000 simulations
  
  ## delta2 stores the death rate of the stem population
  delta2 <- c()
  ## delta_tert stores the death rate of the TERT mutated population
  delta_tert <- c()
  ## delta0 stores the average death rate if approximating growth with a single exponential expansion
  delta0 <- c()
  ## delta_stem stores the death rate of the stem population (like delta0; but since not all values of delta_stem will give valid fits, we need delta0 to store only the working parameters, while delta_stem stores all tested parameters)
  delta_stem <- seq(0.6, 0.999, 0.001)
  ## single-base substitution rate times 2 (2 TERT promoter mutations possible)
  mut.rate <- 2*2.66*10^-9
  
  ## mean and var store the mean and the variance from the simulation
  mean <- c()
  var <- c()
  ## i: average delta if modeling growth as a single exponential expansion (we do this approximation to solve for Tres, the end timepoint/resection)
  for(i in delta_stem - 0.001){
    print(i)
    ## j: death rate of stem population (TERT-wt). j must be larger than i, because i is the average of the death rate of the TERT-wt and mutated fraction and the delta_tert <= delta2
    for(j in seq(i+0.001, 0.999, 0.001)){
      
      ## randomly sample 1000 numbers between 0 and 1. With these random numbers, assess the time of the TERT mutation in the founder population (death rate = j)
      unifs <- runif(1000, 0, 1)
      unifs <- sapply(unifs, function(x){uniroot(P2, delta=j, unif=x, interval = c(0, 100000), mut.rate = 2*2.66*10^-9)$root})
      ## determine tres (end point /timepoint of resection) from approximation by single-exponential expansion (see above)
      tres <- 9*log(10)/(1-i)

      ## delta: death rate of TERT mutated population; can now be determined from tres and using the expected time of the TERT mutation as the average of the simulated times (unifs).
      delta <- 1 - log(exp((1-i)*tres) - exp((1-j)*tres))/(tres - mean(unifs))
      ## death rate mustn't exceed 1 (by definition tumors grow, thus lambda > delta)
      if(delta > 1){next}
      ## store i, the approximate death rate if modeling tumor growth with a single expansion
      delta0 <- c(delta0, i)
      
      ## now determine tres exactly for each evaluation; with the help of j (death rate in the stem) and delta (death rate in the TERT mutated fraction)
      tres.ind <- sapply(unifs, function(x){
        uniroot(tumor.size2, delta = j, delta2 = delta, start2 = x, interval = c(0, x+9*log(10)/(1-j)))$root})
      ## determine the relative tert proportion at tres
      rel.tert <- relative.tert.proportion(tres.ind, delta = j, delta_tert=delta, start_tert=unifs)
      
      ## mean and variance of the TERT mutated tumor fraction at death rates j and delta
      mean <- c(mean, mean(rel.tert))
      var <- c(var, var(rel.tert))
      
      ## store the death rate of the stem population
      delta2 <- c(delta2, j)
      ## and the death rate of the TERT mutated fraction
      delta_tert <- c(delta_tert, mean(delta))
      ## only accept death rates bewteen 0 and 1
      delta_tert[delta_tert>1 | delta_tert < 0] <- NA

    }
  }
  
  
  
  ## with this we obtain simulated values of the mean and variance of the TERT mutated tumor fraction at different parameters
  

  ## We now want to compare the simulations to the data
  
  ## provide data (inferred TERT mutated tumor fraction). We excluded samples with non-canonical TERT promoter mutations and ATRX mutations since it is unclear
  ## whether the same selective processes act in these tumors
  measured.data <- c(0.6845086, 0.8947414, 0.7876193, 1, 0.808011, 1, 0.8044239, 1, 1, 0.6701654, 1, 1 , 1, 1, 0.8795581, 0.8415703, 1, 1, 0.8063186)
  
  
  ## variance of mean and variance of the measured data by bootstrapping
  
  bootstrapped.means <- c()
  bootstrapped.var <- c()
  
  for(i in 1:10000){
    tmp <- sample(measured.data, length(measured.data), replace = T)
    bootstrapped.means <- c(bootstrapped.means, mean(tmp))
    bootstrapped.var <- c(bootstrapped.var, var(tmp))
  }
  
  ## least squares when comparing measured data and simulation
  RSS <- ((mean - mean(measured.data))/sd(bootstrapped.means))^2 + ((var - var(measured.data))/sd(bootstrapped.var))^2
  RSS <- RSS - min(RSS, na.rm=T)
  

  library(scales)
  library(ggplot2)

  library(scales)
  library(ggplot2)
  ## plot the death rate of the TERT-wt (delta2) and the TERT-mutated (delta_tert) tumor fractions; color based on RSS
  data <- data.frame(x=delta2, y=delta_tert, RSS=(RSS))[mean!=0 & !is.na(delta_tert) & !is.na(var),]
  ggplot(data, aes(x=x, y=y, col=(RSS), na.rm=T)) + 
    geom_point(size=0.5) +
    scale_color_gradientn(colors=c("firebrick","orange",  "yellow"), values = rescale(c(0,5.99,max(data[,3]))), ## Chi2 of 5.99 corresponds to the 95% confidence interval with 2 free parameters (delta1, delta_tert)
                          guide = "colorbar", limits=c(0,max(data[,3]))) + scale_x_continuous(name = "delta_0") + 
    scale_y_continuous(name = "delta_TERT") +
    theme_bw(base_size = 14) +
    geom_abline(slope = 1, intercept = 0, col = "grey", linetype=2) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_hline(yintercept=0.69)  ## lower boundary of 69% turnover as estimated in Suppl. Fig. 9c
  
  
  ## best parameters; delta_tert >= 0.92 comes from the overall analysis of cellular death in the tumor sample (see above)
  delta2. <- delta2[delta_tert >= 0.92 ][which(RSS[delta_tert>=0.92]<=5.99)[1]]
  delta0. <- delta0[delta_tert >= 0.92 ][which(RSS[delta_tert>=0.92]<=5.99)[1]]
  delta_tert. <- delta_tert[delta_tert >= 0.92 ][which(RSS[delta_tert>=0.92]<=5.99)[1]]
  
  ## simulate distribution
  unifs <- runif(1000, 0, 1)
  unifs <- sapply(unifs, function(x){uniroot(P2, delta=delta2., mut.rate = 2*2.66*10^-9, unif=x, interval = c(0, 1000000))$root})
  tres.ind <- sapply(unifs, function(x){
    uniroot(tumor.size2, delta = delta2., delta2 = delta_tert., start2 = x, interval = c(0, x+9*log(10)/(1-delta2.)))$root})
  rel.tert <- relative.tert.proportion(tres.ind, delta = delta2., delta2=delta_tert., start2=unifs)
  
  to.plot <- data.frame(x=sort(measured.data), y=sapply(sort(measured.data), function(x){sum(measured.data <= x)})/21)
  to.plot <- rbind(c(0,0), to.plot)
  to.plot.model <- data.frame(x=sort(rel.tert), y=sapply(sort(rel.tert), function(x){sum(rel.tert <= x)})/1000)
  to.plot.model <- rbind(c(0,0), to.plot.model)
  
  ggplot(to.plot) + geom_step(aes(x=x, y=y), col="grey")+ theme_grey(base_size = 14) + scale_x_continuous(name = "TERT mutated tumor fraction", limits=c(0, 1)) +
    scale_y_continuous(name = " % Tumors", breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "25", "50", "75", "100")) +
    theme_bw(base_size = 14) + geom_step(data=to.plot.model, aes(x=x, y=y), col="red")+
    scale_fill_manual(values=c(Data = "darkgrey")) + scale_color_manual(values=c(Model="red"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  legend("topleft", lty = c(1, NA), fill = c(NA, "grey"), col=c("red", "grey"), legend=c("model", "data"))
  
  
  ## repeat the same procedure for the four-fold mutation rate
  
  ## stores the death rate of the stem population
  delta2 <- c()
  ## stores the death rate of the TERT mutated population
  delta_tert <- c()
  ## stores the average death rate if approximating growth with a single exponential expansion
  delta0 <- c()
  ## same as delta0; but since not all delta_stems will give valid fits, we need delta0 to store only the working parameters
  delta_stem <- seq(0.6, 0.999, 0.001)
  ## 4x mutation rate
  mut.rate <- 2*2.66*10^-9*4
  
  ## mean and var store the mean and the variance from the simulation
  mean <- c()
  var <- c()
  ## i: average delta if modeling growth as a single exponential expansion (we do this approximation to solve for Tres, the end timepoint/resection)
  for(i in delta_stem - 0.001){
    print(i)
    ## j: death rate of stem population (TERT-wt). j must be larger than i, because i is the average of the death rate of the TERT-wt and mutated fraction and the delta_tert <= delta2
    for(j in seq(i+0.001, 0.999, 0.001)){
      
      ## randomly sample 1000 numbers between 0 and 1. With these random numbers, assess the time of the TERT mutation in the founder population (death rate = j)
      unifs <- runif(1000, 0, 1)
      unifs <- sapply(unifs, function(x){uniroot(P2, delta=j, unif=x, interval = c(0, 100000), mut.rate = 2*2.66*10^-9)$root})
      ## determine tres (end point /timepoint of resection) from approximation by single-exponential expansion (see above)
      tres <- 9*log(10)/(1-i)
      
      ## delta: death rate of TERT mutated population; can now be determined from tres and using the expected time of the TERT mutation as the average of the simulated times (unifs).
      delta <- 1 - log(exp((1-i)*tres) - exp((1-j)*tres))/(tres - mean(unifs))
      ## death rate mustn't exceed 1 (by definition tumors grow, thus lambda > delta)
      if(delta > 1){next}
      ## store i, the approximate death rate if modeling tumor growth with a single expansion
      delta0 <- c(delta0, i)
      
      ## now determine tres exactly for each evaluation; with the help of j (death rate in the stem) and delta (death rate in the TERT mutated fraction)
      tres.ind <- sapply(unifs, function(x){
        uniroot(tumor.size2, delta = j, delta2 = delta, start2 = x, interval = c(0, x+9*log(10)/(1-j)))$root})
      ## determine the relative tert proportion at tres
      rel.tert <- relative.tert.proportion(tres.ind, delta = j, delta_tert=delta, start_tert=unifs)
      
      ## mean and variance of the TERT mutated tumor fraction at death rates j and delta
      mean <- c(mean, mean(rel.tert))
      var <- c(var, var(rel.tert))
      
      ## store the death rate of the stem population
      delta2 <- c(delta2, j)
      ## and the death rate of the TERT mutated fraction
      delta_tert <- c(delta_tert, mean(delta))
      ## only accept death rates bewteen 0 and 1
      delta_tert[delta_tert>1 | delta_tert < 0] <- NA
      
    }
  }
  
  
  ## variance of mean and variance of the measured data by bootstrapping
  
  bootstrapped.means <- c()
  bootstrapped.var <- c()
  
  for(i in 1:10000){
    tmp <- sample(measured.data, length(measured.data), replace = T)
    bootstrapped.means <- c(bootstrapped.means, mean(tmp))
    bootstrapped.var <- c(bootstrapped.var, var(tmp))
  }
  
  RSS <- ((mean - mean(measured.data))/sd(bootstrapped.means))^2 + ((var - var(measured.data))/sd(bootstrapped.var))^2
  RSS <- RSS - min(RSS, na.rm=T)
  

  data <- data.frame(x=delta2, y=delta_tert, RSS=(RSS))[mean!=0 & !is.na(delta_tert) & !is.na(var),]
  ggplot(data, aes(x=x, y=y, col=(RSS), na.rm=T)) + 
    geom_point(size=0.5) +
    scale_color_gradientn(colors=c("firebrick","orange",  "yellow"), values = rescale(c(0,5.99,max(data[,3]))),
                          guide = "colorbar", limits=c(0,max(data[,3]))) + scale_x_continuous(name = "delta_0") + 
    scale_y_continuous(name = "delta_TERT") +
    theme_bw(base_size = 14) +
    geom_abline(slope = 1, intercept = 0, col = "grey", linetype=2) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_hline(yintercept=0.69)  ## lower boundary of 69% turnover as estimated in Suppl. Fig. 9c

  
  
  ## best parameters: delta_tert >= 0.69 comes from the overall analysis of cellular death in the tumor sample (see above)
  delta2. <- delta2[delta_tert >= 0.69 ][which(RSS[delta_tert>=0.69]<=5.99)[1]]
  delta0. <- delta0[delta_tert >= 0.69 ][which(RSS[delta_tert>=0.69]<=5.99)[1]]
  delta_tert. <- delta_tert[delta_tert >= 0.69 ][which(RSS[delta_tert>=0.69]<=5.99)[1]]
  
  ## simulate distribution
  unifs <- runif(1000, 0, 1)
  unifs <- sapply(unifs, function(x){uniroot(P2, delta=delta2., mut.rate = 2*4*2.66*10^-9, unif=x, interval = c(0, 1000000))$root})
  tres.ind <- sapply(unifs, function(x){
    uniroot(tumor.size2, delta = delta2., delta2 = delta_tert., start2 = x, interval = c(0, x+9*log(10)/(1-delta2.)))$root})
  rel.tert <- relative.tert.proportion(tres.ind, delta = delta2., delta2=delta_tert., start2=unifs)
  
  to.plot <- data.frame(x=sort(measured.data), y=sapply(sort(measured.data), function(x){sum(measured.data <= x)})/21)
  to.plot <- rbind(c(0,0), to.plot)
  to.plot.model <- data.frame(x=sort(rel.tert), y=sapply(sort(rel.tert), function(x){sum(rel.tert <= x)})/1000)
  to.plot.model <- rbind(c(0,0), to.plot.model)
  
  ggplot(to.plot) + geom_step(aes(x=x, y=y), col="grey")+ theme_grey(base_size = 14) + scale_x_continuous(name = "TERT mutated tumor fraction", limits=c(0, 1)) +
    scale_y_continuous(name = " % Tumors", breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "25", "50", "75", "100")) +
    theme_bw(base_size = 14) + geom_step(data=to.plot.model, aes(x=x, y=y), col="red")+
    scale_fill_manual(values=c(Data = "darkgrey")) + scale_color_manual(values=c(Model="red"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  legend("topleft", lty = c(1, NA), fill = c(NA, "grey"), col=c("red", "grey"), legend=c("model", "data"))
  

 