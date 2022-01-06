library(tidyverse)
library(rstanarm)

nbacc_dose_response_data <- read.csv("nbacc_sero_fever.csv")

pooled_probit <- glm(fever ~ logdose,family=binomial(link="probit"),data=nbacc_dose_response_data)
variant_probit <- glm(fever ~ variant + logdose,family=binomial(link="probit"),data=nbacc_dose_response_data)

teststat <- -2 * (as.numeric(logLik(pooled_probit))-as.numeric(logLik(variant_probit)))
p.val <- pchisq(teststat, df = 1, lower.tail = FALSE)

p <- ggplot(nbacc_dose_response_data,aes(logdose,fever,color=variant)) + 
  geom_point() + 
  geom_smooth(method = "glm", fullrange=T,
              se=F, method.args = list(family = binomial(link="probit"))) + 
  theme_bw()
p


variant_probit_bayes <- stan_glm(fever ~ variant + logdose,family=binomial(link="probit"),data=nbacc_dose_response_data)

pplot<-plot(variant_probit_bayes, "areas", prob = 0.95, prob_outer = 1)
pplot+ geom_vline(xintercept = 0)
