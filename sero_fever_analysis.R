library(tidyverse)
library(VGAM)
library(lmtest)
nbacc_dose_response_data <- read.csv("nbacc_sero_fever.csv")

nbacc_dose_response_data[nbacc_dose_response_data == '-'] <- 0 
nbacc_dose_response_data[nbacc_dose_response_data == '+'] <- 1 
nbacc_dose_response_data[nbacc_dose_response_data == 'ND'] <- NA 

nbacc_dose_response_data$PRNT50_binary <- as.numeric(nbacc_dose_response_data$PRNT50_binary)
nbacc_dose_response_data$PRNT80_binary <- as.numeric(nbacc_dose_response_data$PRNT80_binary)
nbacc_dose_response_data$IgG_ELISA_kit <- as.numeric(nbacc_dose_response_data$IgG_ELISA_kit)
nbacc_dose_response_data$IgG_ELISA_rebase <- as.numeric(nbacc_dose_response_data$IgG_ELISA_rebase)
nbacc_dose_response_data$Fever_binary <- as.numeric(nbacc_dose_response_data$Fever_binary)
nbacc_dose_response_data$logdose <- as.numeric(nbacc_dose_response_data$logdose)

### GLM based modeling

### Bivariate VGLM based modeling

wa1_biv_probit <- vglm( cbind(PRNT50_binary,Fever_binary) ~  logdose,
                              family=binom2.or(lmu="probitlink",zero = "or",exchangeable=F),
                              data=filter(nbacc_dose_response_data,variant=='WA-1'),
                              control=vglm.control(stepsize = 0.5,maxit = 1000))

wa1_biv_probit_2 <- vglm( cbind(PRNT50_binary,Fever_binary) ~  logdose,
                        family=binom2.or(lmu="probitlink",zero = "or",exchangeable=T),
                        data=filter(nbacc_dose_response_data,variant=='WA-1'),
                        control=vglm.control(stepsize = 0.5,maxit = 1000))


pchisq(-2*(logLik(wa1_biv_probit_2) - logLik(wa1_biv_probit)), df = 2, lower.tail = FALSE)

summary(wa1_biv_probit,score0=T)
score.stat(wa1_biv_probit)

#### Parameter estimates
wa1_pars <- coef(wa1_biv_probit)

wa1_sero_logD50_mle <- -wa1_pars[1]/wa1_pars[4]
wa1_sero_logD50_mle
10^wa1_sero_logD50_mle


wa1_fever_logD50_mle <- -wa1_pars[2]/wa1_pars[5]
wa1_fever_logD50_mle
10^wa1_fever_logD50_mle

### Gamma variant
gamma_biv_probit <- vglm( cbind(PRNT50_binary,Fever_binary) ~  logdose,
                        family=binom2.or(lmu="probitlink",zero = "or",exchangeable=F),
                        data=filter(nbacc_dose_response_data,variant=='Gamma'),
                        control=vglm.control(stepsize = 0.5,maxit = 1000))

gamma_biv_probit_2 <- vglm( cbind(PRNT50_binary,Fever_binary) ~  logdose,
                          family=binom2.or(lmu="probitlink",zero = "or",exchangeable=T),
                          data=filter(nbacc_dose_response_data,variant=='Gamma'),
                          control=vglm.control(stepsize = 0.5,maxit = 1000))

pchisq(-2*(logLik(gamma_biv_probit_2) - logLik(gamma_biv_probit)), df = 2, lower.tail = FALSE)

summary(wa1_biv_probit,score0=T)
score.stat(wa1_biv_probit)

### Parameter
gamma_pars <- coef(gamma_biv_probit)

gamma_sero_logD50_mle <- -gamma_pars[1]/gamma_pars[4]
gamma_sero_logD50_mle
10^gamma_sero_logD50_mle

gamma_fever_logD50_mle <- -gamma_pars[2]/gamma_pars[5]
gamma_fever_logD50_mle
10^gamma_fever_logD50_mle

###

#Gamma vs WA-1 fever
pooled_probit_fever <- glm(Fever_binary ~ logdose,family=binomial(link="probit"),data=nbacc_dose_response_data)
variant_probit_fever <- glm(Fever_binary ~ variant + logdose,family=binomial(link="probit"),data=nbacc_dose_response_data)

wa1_probit_fever <- glm(Fever_binary ~ logdose,family=binomial(link="probit"),data=filter(nbacc_dose_response_data,variant=='WA-1'))
gamma_probit_fever <- glm(Fever_binary ~  logdose,family=binomial(link="probit"),data=filter(nbacc_dose_response_data,variant=='Gamma'))

variant_fever_probit_pars <- coef(variant_probit_fever)
wa1_fever_probit_pars <- coef(wa1_probit_fever)
gamma_fever_probit_pars <- coef(gamma_probit_fever)

-wa1_fever_probit_pars[1]/wa1_fever_probit_pars[2]
-gamma_fever_probit_pars[1]/gamma_fever_probit_pars[2]

-(variant_fever_probit_pars[1] + variant_fever_probit_pars[2])/variant_fever_probit_pars[3]
-variant_fever_probit_pars[1]/variant_fever_probit_pars[3]

teststat <- -2 * (as.numeric(logLik(pooled_probit_fever))-as.numeric(logLik(variant_probit_fever)))
p.val <- pchisq(teststat, df = 1, lower.tail = FALSE)
p.val

#Gamma vs WA-1 sero
pooled_probit_sero <- glm(PRNT50_binary ~ logdose,family=binomial(link="probit"),data=nbacc_dose_response_data)
variant_probit_sero <- glm(PRNT50_binary ~ variant + logdose,family=binomial(link="probit"),data=nbacc_dose_response_data)

wa1_probit_sero <- glm(PRNT50_binary ~ logdose,family=binomial(link="probit"),data=filter(nbacc_dose_response_data,variant=='WA-1'))
gamma_probit_sero <- glm(PRNT50_binary ~  logdose,family=binomial(link="probit"),data=filter(nbacc_dose_response_data,variant=='Gamma'))

variant_sero_probit_pars <- coef(variant_probit_sero)
wa1_sero_probit_pars <- coef(wa1_probit_sero)
gamma_sero_probit_pars <- coef(gamma_probit_sero)

-wa1_sero_probit_pars[1]/wa1_sero_probit_pars[2]
-gamma_sero_probit_pars[1]/gamma_sero_probit_pars[2]

-(variant_sero_probit_pars[1] + variant_sero_probit_pars[2])/variant_sero_probit_pars[3]
-variant_sero_probit_pars[1]/variant_sero_probit_pars[3]

teststat <- -2 * (as.numeric(logLik(pooled_probit_sero))-as.numeric(logLik(variant_probit_sero)))
p.val <- pchisq(teststat, df = 1, lower.tail = FALSE)
p.val



### Parking lot
variant_bivariate_probit <- vglm( cbind(PRNT50_binary,Fever_binary) ~  variant + logdose,
                                  family=binom2.or(lmu="probitlink",zero='or',exchangeable = F),
                                  data=nbacc_dose_response_data,
                                  control=vglm.control(stepsize = 0.1,maxit = 1000))

pooled_bivariate_probit <- vglm( cbind(PRNT50_binary,Fever_binary) ~  logdose,
                                family=binom2.or(lmu="probitlink",zero='or',exchangeable = F),
                                data=nbacc_dose_response_data,
                                control=vglm.control(stepsize = 0.1,maxit = 1000))


teststat <- -2 * (as.numeric(logLik(pooled_bivariate_probit))-as.numeric(logLik(variant_bivariate_probit)))
p.val <- pchisq(teststat, df = length(coef(variant_bivariate_probit)) - length(coef(pooled_bivariate_probit)) , lower.tail = FALSE)
p.val 
 

variant_pars <- coef(variant_bivariate_probit)

wa1_sero_logD50_mle <- -(variant_pars[1] + variant_pars[4])/(variant_pars[6])
gamma_sero_logD50_mle <- -(variant_pars[1])/variant_pars[6]

10^wa1_sero_logD50_mle
10^gamma_sero_logD50_mle

wa1_fever_logD50_mle <- -(variant_pars[2] + variant_pars[5])/(variant_pars[7])
gamma_fever_logD50_mle <- -(variant_pars[2])/variant_pars[7]

10^wa1_fever_logD50_mle
10^gamma_fever_logD50_mle
