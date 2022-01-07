library(tidyverse)

nbacc_dose_response_data <- read.csv("nbacc_sero_fever.csv")

wa1_data <- filter(nbacc_dose_response_data,variant=='WA-1')
gamma_data <- filter(nbacc_dose_response_data,variant=='Gamma')

wa1_model <- glm(fever ~ logdose,family=binomial(link="probit"),data=wa1_data)
gamma_model <- glm(fever ~ logdose,family=binomial(link="probit"),data=gamma_data)

wa1_pars <- coef(wa1_model)
gamma_pars <- coef(gamma_model)

wa1_logD50_mle <- -wa1_pars[1]/wa1_pars[2]
gamma_logD50_mle <- -gamma_pars[1]/gamma_pars[2]


hypothetical_logD50_diff <- data.frame()
for(logD50_diff in seq(0,-1,-0.05)){
  
  hypothetical_gamma_pars <- c(-wa1_pars[2]*(wa1_logD50_mle + logD50_diff),wa1_pars[2])
  
  res <- do.call("rbind",lapply(1:500,function(i){
    
    hypothetical_gamma_res <- rbinom(nrow(gamma_data),1,prob=pnorm(hypothetical_gamma_pars[1] + hypothetical_gamma_pars[2]*gamma_data$logdose))
    
    hypothetical_combined <- rbind(
      data.frame(variant='WA-1',logdose=wa1_data$logdose,fever=wa1_data$fever),
      data.frame(variant='Gamma',logdose=gamma_data$logdose,fever=hypothetical_gamma_res)
    )
    
    pooled_probit <- glm(fever ~ logdose,family=binomial(link="probit"),data=hypothetical_combined)
    variant_probit <- glm(fever ~ variant + logdose,family=binomial(link="probit"),data=hypothetical_combined)
    
    teststat <- -2 * (as.numeric(logLik(pooled_probit))-as.numeric(logLik(variant_probit)))
    p.val <- pchisq(teststat, df = 1, lower.tail = FALSE)
    
    return(c(logD50_diff,teststat,p.val))
  }))
  
  hypothetical_logD50_diff <- rbind(hypothetical_logD50_diff,res)
}

power_alpha_05 <- hypothetical_logD50_diff %>% group_by(V1) %>% summarise(power= mean(V3<=0.05))

plot(-power_alpha_05$V1,power_alpha_05$power,xlab="Log10(D50) Difference",ylab="Power")
abline(v=(wa1_logD50_mle - gamma_logD50_mle),col='red')
