library(medflex)
library(mediation)

# load data set
data(framing)

# drop empty factor levels to avoid trouble with sandwich estimator in medflex later on
framing <- droplevels(framing)

# more details on the data set
?framing

# MSM for causal effect of framing
summary(glm(immigr ~ tone * eth, data = framing)) 
# accounting for the 2 x 2 factorial design
# only negatively framed news story that involves Latino immigrants invokes considerably more negative attitude towards immigration
# => focus on contrast between this condition and three other conditions
fitY <- glm(immigr ~ treat, data = framing)
summary(fitY)

#----

# inspecting the mediator variable
with(framing, prop.table(table(treat, anx), margin = 1)) 
# anxiety seems to be reverse coded (as it has a negative effect on immigr)
# recode as numeric variable such that lower values correspond to lower levels of anxiety
framing$anxiety <- as.numeric(framing$anx)

## fit nuisance working models
# outcome working model (treating immigr as a continuous normally distributed variable)
fitYM <- glm(immigr ~ treat + anxiety + gender + age + educ + income, data = framing)

# mediator working model (treating anxiety as a continuous normally distributed variable)
fitM <- glm(anxiety ~ treat + gender + age + educ + income, data = framing)

## preparing an expanded data set
# first create ID variable
framing <- data.frame(ID = 1:nrow(framing), framing)

# duplicate the data set and create auxiliary variables
expData <- data.frame(replicate = rep(1:2, times = nrow(framing)),
                      framing[rep(framing$ID, each = 2), ])


## weighting-based approach ----
weightData <- within(expData, 
                     {
                       treat0 <- treat
                       treat1 <- ifelse(replicate == 1, treat, 1-treat)
                     })
head(weightData)
# calculate weights based on generalized propensity scores
num <- with(weightData, dnorm(anxiety, mean = predict(fitM, newdata = within(weightData, treat <- treat1)), sd = sqrt(summary(fitM)$dispersion)))
denom <- with(weightData, dnorm(anxiety, mean = predict(fitM, newdata = within(weightData, treat <- treat0)), sd = sqrt(summary(fitM)$dispersion)))
w <- num/denom
head(data.frame(weightData, w))
# fit natural effect model
glm(immigr ~ treat0 + treat1, data = weightData, weights = w)
# corresponding standard errors called by summary() function will generally be invalid
# as clustering and uncertainty related to estimating weights is not taken into account

# in medflex
weightData <- neWeight(anxiety ~ factor(treat) + gender + age + educ + income, data = framing)
# important to explicitly code treat as a factor in order to get correct data expansion (otherwise treated as a continuous variable)
neModW <- neModel(immigr ~ treat0 + treat1, expData = weightData, se = "robust")
summary(neModW)
# valid standard errors for neModel objects => based on robust sandwich variance-covariance matrix
# does take into account clustering and uncertainty related to estimating weights
# same holds for other functions adapted to neModel objects
?`neModel-methods`

neEffdecomp(neModW)
plot(neModW)


## imputation-based approach ----
impData <- within(expData, 
                  {
                    treat1 <- treat
                    treat <- treat0 <- ifelse(replicate == 1, treat, 1-treat)
                  })
head(impData)
# impute counterfactual outcomes
impData$immigr <- predict(fitYM, newdata = impData)
# fit natural effect model
glm(immigr ~ treat0 + treat1, data = impData)

# in medflex
impData <- neImpute(immigr ~ factor(treat) + anxiety + gender + age + educ + income, data = framing)
neModI <- neModel(immigr ~ treat0 + treat1, expData = impData, se = "robust")
summary(neModI)
neEffdecomp(neModI)
plot(neModI)


## equivalence with traditional LSEM estimators under strictly linear working models

neModIss <- neModel(immigr ~ treat0 + treat1 + gender + age + educ + income, expData = impData, se = "robust")
# stratum-specific natural effect model

# (controlled) direct effect
theta <- as.list(coef(fitYM))
theta$treat
coef(neModI)["treat01"]
coef(neModIss)["treat01"]

# difference-in-coefficients estimator for indirect effect
gamma <- as.list(coef(fitY))
gamma$treat - theta$treat
coef(neModI)["treat11"]

fitYC <- glm(immigr ~ treat + gender + age + educ + income, data = framing)
# yields different estimator for framing effect (both stratum-specific and population-average, under assumption of no effect modification)
coef(fitYC)["treat"] - theta$treat
coef(neModIss)["treat11"]

# product-of-coefficients estimator for indirect effect
beta <- as.list(coef(fitM))
beta$treat * theta$anxiety
coef(neModIss)["treat11"]


## non-linearities ----

## allow for mediated interaction ----
neModIia <- neModel(immigr ~ treat0 * treat1, expData = impData, se = "robust")
summary(neModIia)
neEffdecomp(neModIia)
# interaction term attenuates to zero if imputation model is not adapted accordingly

impData <- neImpute(immigr ~ factor(treat) * anxiety + gender + age + educ + income, data = framing)
neModIia <- neModel(immigr ~ treat0 * treat1, expData = impData, se = "robust")
summary(neModIia)
neEffdecomp(neModIia)
plot(neModIia)


## stratum-specific effects and effect modification ----
neModIiass <- neModel(immigr ~ treat0 * treat1 + gender + age + educ + income, expData = impData, se = "robust")
summary(neModIiass)
neEffdecomp(neModIiass)
neEffdecomp(neModIiass, covLev = c(gender = "female", age = 50, educ = "high school", income = 10))
# no difference because no effect modification allowed for in natural effect parameterization

## direct application of mediation formula
fitYM <- glm(immigr ~ treat * anxiety + gender + age + educ + income, data = framing)
theta <- as.list(coef(fitYM))

c(theta$treat + theta$`treat:anxiety`*beta$`(Intercept)`, # theta1 + theta3*beta0
  theta$anxiety*beta$treat, # theta2*beta1
  theta$`treat:anxiety`*beta$treat) # theta3*beta1
# stratum-specific natural effect model that is compatible with parameterization resulting from working model specifications
neModIiass <- neModel(immigr ~ treat0 * treat1 + treat0 * (gender + age + educ + income), expData = impData, se = "robust")
coef(neModIiass)[c("treat01", "treat11", "treat01:treat11")]


# test for effect modification by gender
impData <- neImpute(immigr ~ factor(treat) * anxiety * gender + age + educ + income, data = framing)
# refit imputation model explicitly allowing for effect modification by gender
neModI_G <- neModel(immigr ~ treat0 * treat1 * gender, expData = impData, se = "robust")
summary(neModI_G)
neEffdecomp(neModI_G)
neEffdecomp(neModI_G, covLev = c(gender = "female"))

# hypothesis captured by multiple model parameters => use neLht function to test linear hypothesis
# total effect modified by gender?
testTE_G <- neLht(neModI_G, linfct = c("`treat01:genderfemale` + `treat11:genderfemale` + `treat01:treat11:genderfemale` == 0"))
summary(testTE_G)
# NDE(0) modified by gender?
testNDE0_G <- neLht(neModI_G, linfct = c("`treat01:genderfemale` == 0"))
summary(testNDE0_G)
# NDE(1) modified by gender?
testNDE1_G <- neLht(neModI_G, linfct = c("`treat01:genderfemale` + `treat01:treat11:genderfemale` == 0"))
summary(testNDE1_G)
# NIE(0) modified by gender?
testNIE0_G <- neLht(neModI_G, linfct = c("`treat11:genderfemale` == 0"))
summary(testNIE0_G)
# NIE(1) modified by gender?
testNIE1_G <- neLht(neModI_G, linfct = c("`treat11:genderfemale` + `treat01:treat11:genderfemale` == 0"))
summary(testNIE1_G)
# omnibus test for effect modification of natural direct and indirect effects by gender
testNE_G <- neLht(neModI_G, linfct = c("`treat01:genderfemale` = 0",
                                       "`treat11:genderfemale` = 0",
                                       "`treat01:treat11:genderfemale` == 0"))
summary(testNE_G, test = Chisqtest())


# test for effect modification by education level
impData <- neImpute(immigr ~ factor(treat) * anxiety * educ + age + gender + income, data = framing)
# refit imputation model explicitly allowing for effect modification by education level
neModI_E <- neModel(immigr ~ treat0 * treat1 * educ, expData = impData, se = "robust")
summary(neModI_E)
plot(neEffdecomp(neModI_E, covLev = c(educ = "less than high school")))
plot(neEffdecomp(neModI_E, covLev = c(educ = "high school")))
plot(neEffdecomp(neModI_E, covLev = c(educ = "some college")))
plot(neEffdecomp(neModI_E, covLev = c(educ = "bachelor's degree or higher")))

# total effect modified by education level?
testTE_E <- neLht(neModI_E, linfct = c("`treat01:educhigh school` + `treat11:educhigh school` + `treat01:treat11:educhigh school` == 0",
                                       "`treat01:educsome college` + `treat11:educsome college` + `treat01:treat11:educsome college` == 0",
                                       "`treat01:educbachelor's degree or higher` + `treat11:educbachelor's degree or higher` + `treat01:treat11:educbachelor's degree or higher` == 0"))
summary(testTE_E, test = Chisqtest())
# NDE(0) modified by education level?
testNDE0_E <- neLht(neModI_E, linfct = c("`treat01:educhigh school` == 0",
                                         "`treat01:educsome college` == 0",
                                         "`treat01:educbachelor's degree or higher` == 0"))
summary(testNDE0_E, test = Chisqtest())
# NDE(1) modified by education level?
testNDE1_E <- neLht(neModI_E, linfct = c("`treat01:educhigh school` + `treat01:treat11:educhigh school` == 0",
                                         "`treat01:educsome college` + `treat01:treat11:educsome college` == 0",
                                         "`treat01:educbachelor's degree or higher` + `treat01:treat11:educbachelor's degree or higher` == 0"))
summary(testNDE1_E, test = Chisqtest())
# NIE(0) modified by education level?
testNIE0_E <- neLht(neModI_E, linfct = c("`treat11:educhigh school` == 0",
                                         "`treat11:educsome college` == 0",
                                         "`treat11:educbachelor's degree or higher` == 0"))
summary(testNIE0_E, test = Chisqtest())
# NIE(1) modified by education level?
testNIE1_E <- neLht(neModI_E, linfct = c("`treat11:educhigh school` + `treat01:treat11:educhigh school` == 0",
                                         "`treat11:educsome college` + `treat01:treat11:educsome college` == 0",
                                         "`treat11:educbachelor's degree or higher` + `treat01:treat11:educbachelor's degree or higher` == 0"))
summary(testNIE1_E, test = Chisqtest())
# omnibus test for effect modification of natural direct and indirect effects by education level
testNE_E <- neLht(neModI_E, linfct = c("`treat01:educhigh school` = 0",
                                       "`treat01:educsome college` = 0",
                                       "`treat01:educbachelor's degree or higher` = 0",
                                       "`treat11:educhigh school` = 0",
                                       "`treat11:educsome college` = 0",
                                       "`treat11:educbachelor's degree or higher` = 0",
                                       "`treat01:treat11:educhigh school` = 0",
                                       "`treat01:treat11:educsome college` = 0",
                                       "`treat01:treat11:educbachelor's degree or higher` = 0"))
summary(testNE_E, test = Chisqtest())

# effect modification of NDE(0) and NIE(0) can also easily be verified by Anova function
library(car)
Anova(neModI_E, type = 3)


# If possible, try to aim for a rich imputation model that makes few parametric assumptions
impData <- neImpute(immigr ~ factor(treat) * factor(anxiety) * (gender + income + age) + educ, data = framing)
neModIia <- neModel(immigr ~ treat0 * treat1, expData = impData, se = "robust")
summary(neModIia)

## inverse-odds weighting
fitA <- glm(treat ~ anxiety + gender + age + educ + income, family = binomial, data = framing)
num <- with(weightData, dbinom(as.numeric(treat1) - 1, size = 1, prob = predict(fitA, newdata = weightData, type = "response")))
denom <- with(weightData, dbinom(as.numeric(treat0) - 1, size = 1, prob = predict(fitA, newdata = weightData, type = "response")))
w <- num/denom
head(data.frame(weightData, w))
# fit natural effect model
glm(immigr ~ treat0 * treat1, data = weightData, weights = w)



## linear and logistic natural effect models for a binary outcome ----
framing <- within(framing, immigr_bin <- ifelse(immigr >= 3, 1, 0))

impData <- neImpute(immigr_bin ~ factor(treat) * factor(anxiety) * (gender + income + age) + educ, 
                    family = binomial, data = framing)

# risk difference scale (warranted because natural effect model saturated)
neModRD <- neModel(immigr_bin ~ treat0 * treat1, expData = impData, se = "robust")
summary(neModRD)
plot(neModRD)

# odds ratio scale
neModOR <- neModel(immigr_bin ~ treat0 * treat1, family = binomial, expData = impData, se = "robust")
summary(neModOR)
plot(neModOR, transf = exp)

# effect modification by gender (on risk difference and odds ratio scales)
neModRD_G <- neModel(immigr_bin ~ treat0 * treat1 * gender, expData = impData, se = "robust")
neModOR_G <- neModel(immigr_bin ~ treat0 * treat1 * gender, family = binomial, expData = impData, se = "robust")

summary(neModRD_G)
plot(neEffdecomp(neModRD_G, covLev = c(gender = "male")))
plot(neEffdecomp(neModRD_G, covLev = c(gender = "female")))
# omnibus test
testNE_G <- neLht(neModRD_G, linfct = c("`treat01:genderfemale` = 0",
                                        "`treat11:genderfemale` = 0",
                                        "`treat01:treat11:genderfemale` == 0"))
summary(testNE_G, test = Chisqtest())


summary(neModOR_G)
plot(neEffdecomp(neModOR_G, covLev = c(gender = "male")), transf = exp)
plot(neEffdecomp(neModOR_G, covLev = c(gender = "female")), transf = exp)
# omnibus test
testNE_G <- neLht(neModOR_G, linfct = c("`treat01:genderfemale` = 0",
                                        "`treat11:genderfemale` = 0",
                                        "`treat01:treat11:genderfemale` == 0"))
summary(testNE_G, test = Chisqtest())


# effect modification by age
neModRD_A <- neModel(immigr_bin ~ treat0 * treat1 * age, expData = impData, se = "robust")
neModOR_A <- neModel(immigr_bin ~ treat0 * treat1 * age, family = binomial, expData = impData, se = "robust")

summary(neModRD_A)
plot(neEffdecomp(neModRD_A, covLev = c(age = 20)))
plot(neEffdecomp(neModRD_A, covLev = c(age = 50)))
plot(neEffdecomp(neModRD_A, covLev = c(age = 80)))
range(predict(neModRD_A$neModelFit)) # still within 0-1 range (risk for extrapolation!)
# omnibus test
testNE_A <- neLht(neModRD_A, linfct = c("`treat01:age` = 0",
                                        "`treat11:age` = 0",
                                        "`treat01:treat11:age` = 0"))
summary(testNE_A, test = Chisqtest())

summary(neModOR_A)
plot(neEffdecomp(neModOR_A, covLev = c(age = 20)), transf = exp)
plot(neEffdecomp(neModOR_A, covLev = c(age = 50)), transf = exp)
plot(neEffdecomp(neModOR_A, covLev = c(age = 80)), transf = exp)
# omnibus test
testNE_A <- neLht(neModOR_A, linfct = c("`treat01:age` = 0",
                                        "`treat11:age` = 0",
                                        "`treat01:treat11:age` = 0"))
summary(testNE_A, test = Chisqtest())


## treatment-induced mediator-outcome confounding (multiple interdependent mediators) ---- 

# perceived harm (belief rather than emotion) as alternative mediator
# independent of anxiety (conditional on treatment and baseline covariates)?
library(dagitty)
g <- downloadGraph("dagitty.net/myYTGFa")
plot(g)
framing2 <- within(framing, {
                              gender <- as.numeric(gender)
                              educ <- as.numeric(educ)
                            })
r <- localTests(g, framing2)
r$p.value <- p.adjust(r$p.value)
par(mar = c(5, 20, 4, 2) + 0.1)
plotLocalTestResults(r)
par(mar = c(5, 4, 4, 2) + 0.1)
# data seems to indicate that anxiety and perceived harm are associated conditional on treatment and baseline covariates
# Following Imai et al (2013) we assume that perceived harm may affect anxiety rather than vice versa (and that there is no unmeasured confounding between these two mediators)


## sequential approach to recover partial indirect effect 
# i.e. the component that is solely mediated by anxiety (or the effect mediated by treatment-induced changes in anxiety as far as these aren't induced by preceding changes in perceived harm)

# joint mediation
impDataLM <- neImpute(immigr ~ treat * p_harm * anxiety + gender + age + educ + income, nMed = 2, data = framing) 
neModLM <- neModel(immigr ~ treat0 * treat1, expData = impDataLM, se = "robust")
summary(neModLM)
neEffdecomp(neModLM)

# mediation by perceived economic harm
impDataL <- neImpute(immigr ~ treat * p_harm + gender + age + educ + income, data = framing) 
neModL <- neModel(immigr ~ treat0 * treat1, expData = impDataL, se = "robust")
summary(neModL)
neEffdecomp(neModL)

# eta3 + eta5
coef(neEffdecomp(neModLM))["total indirect effect"] - coef(neEffdecomp(neModL))["total indirect effect"]
coef(neEffdecomp(neModL))["pure direct effect"] - coef(neEffdecomp(neModLM))["pure direct effect"]

# eta3 + eta6
coef(neEffdecomp(neModL))["total direct effect"] - coef(neEffdecomp(neModLM))["total direct effect"]
coef(neEffdecomp(neModLM))["pure indirect effect"] - coef(neEffdecomp(neModL))["pure indirect effect"]



## direct approach to recover partial indirect effect (option A: weigh for L)
# working models
fitYLM <- glm(immigr ~ treat * p_harm * anxiety + gender + age + educ + income, data = framing)
fitL <- glm(p_harm ~ treat + gender + age + educ + income, data = framing)

# duplicate the data set and create auxiliary variables
expData <- data.frame(replicate = rep(1:4, times = nrow(framing)),
                      framing[rep(framing$ID, each = 4), ])

expData <- within(expData, 
                  {
                    treat0 <- ifelse(replicate %in% c(1,3), treat, 1-treat)
                    treat1 <- ifelse(replicate %in% c(1,2), treat, 1-treat)
                    treat2 <- treat
                    treat <- treat0
                  })
head(expData[, c("treat0", "treat1", "treat2")], 8)

# impute counterfactuals
expData$immigr <- predict(fitYLM, newdata = expData)

# calculate weights based on generalized propensity scores
num <- with(expData, dnorm(p_harm, mean = predict(fitL, newdata = within(expData, treat <- treat1)), sd = sqrt(summary(fitL)$dispersion)))
denom <- with(expData, dnorm(p_harm, mean = predict(fitL, newdata = within(expData, treat <- treat2)), sd = sqrt(summary(fitL)$dispersion)))
w <- num/denom
head(data.frame(expData, w))

# natural effect model
neMod3A <- glm(immigr ~ treat0 * treat1 * treat2, data = expData, weights = w)
neMod3A



## direct approach to recover partial indirect effect (option B: weigh for M)
# working models
fitYLM <- glm(immigr ~ treat * p_harm * anxiety + gender + age + educ + income, data = framing)
fitM <- glm(anxiety ~ factor(treat) + p_harm + gender + age + educ + income, data = framing)

# duplicate the data set and create auxiliary variables
expData <- data.frame(replicate = rep(1:4, times = nrow(framing)),
                      framing[rep(framing$ID, each = 4), ])

expData <- within(expData, 
                  {
                    treat0 <- ifelse(replicate %in% c(1,3), treat, 1-treat)
                    treat2 <- ifelse(replicate %in% c(1,2), treat, 1-treat)
                    treat1 <- treat
                    treat <- treat0
                  })
head(expData[, c("treat0", "treat1", "treat2")], 8)

# impute counterfactuals
expData$immigr <- predict(fitYLM, newdata = expData)

# calculate weights based on generalized propensity scores
num <- with(expData, dnorm(anxiety, mean = predict(fitM, newdata = within(expData, treat <- treat2)), sd = sqrt(summary(fitM)$dispersion)))
denom <- with(expData, dnorm(anxiety, mean = predict(fitM, newdata = within(expData, treat <- treat1)), sd = sqrt(summary(fitM)$dispersion)))
w <- num/denom
head(data.frame(expData, w))

# natural effect model
neMod3B <- glm(immigr ~ treat0 * treat1 * treat2, data = expData, weights = w)
neMod3B





#### weighting-based estimator for pure and total partial indirect effect

# eta3 + eta5
sum(coef(neMod3A)[c("treat2", "treat0:treat2")])
sum(coef(neMod3B)[c("treat2", "treat0:treat2")])

# eta3 + eta6
sum(coef(neMod3A)[c("treat2", "treat1:treat2")])
sum(coef(neMod3B)[c("treat2", "treat1:treat2")])

# pure partial indirect effect (eta3)
coef(neMod3A)["treat2"]
coef(neMod3B)["treat2"]

# total partial indirect effect (eta3 + eta5 + eta6 + eta7)
sum(coef(neMod3A)[c("treat2", "treat0:treat2", "treat1:treat2", "treat0:treat1:treat2")])
sum(coef(neMod3B)[c("treat2", "treat0:treat2", "treat1:treat2", "treat0:treat1:treat2")])

# simpler NE model
weightData2 <- neWeight(anxiety ~ factor(treat) + p_harm + gender + age + educ + income, data = framing)
neModW2 <- neModel(immigr ~ treat0 * treat1, expData = weightData2, se = "robust")
summary(neModW2)

# pure partial indirect effect (beta2)
coef(neModW2)["treat11"]

# total partial indirect effect (beta2 + beta3)
sum(coef(neModW2)[c("treat11", "treat01:treat11")])




# obtain bootstrap standard errors
library(boot)

bootFun <- function(data, index) {
  dat <- data[index, ]
  
  fitYLM <- glm(immigr ~ treat * p_harm * anxiety + gender + age + educ + income, data = dat)
  fitL <- glm(p_harm ~ treat + gender + age + educ + income, data = dat)
  
  dat$ID <- 1:nrow(dat)
  
  expData <- data.frame(replicate = rep(1:4, times = nrow(dat)),
                        dat[rep(dat$ID, each = 4), ])
  
  expData <- within(expData, 
                    {
                      treat0 <- ifelse(replicate %in% c(1,3), treat, 1-treat)
                      treat1 <- ifelse(replicate %in% c(1,2), treat, 1-treat)
                      treat2 <- treat
                      treat <- treat0
                    })
  
  expData$immigr <- predict(fitYLM, newdata = expData)
  
  num <- with(expData, dnorm(p_harm, mean = predict(fitL, newdata = within(expData, treat <- treat1)), sd = sqrt(summary(fitL)$dispersion)))
  denom <- with(expData, dnorm(p_harm, mean = predict(fitL, newdata = within(expData, treat <- treat2)), sd = sqrt(summary(fitL)$dispersion)))
  w <- num/denom
  
  neMod <- glm(immigr ~ treat0 * treat1 * treat2, data = expData, weights = w)
  return(coef(neMod))
}

bootSE <- boot(data = framing, statistic = bootFun, R = 1000)

confint <- sapply(1:length(bootSE$t0), FUN = function(x) boot.ci(bootSE, conf = 0.95, type = "norm", index = x)$normal[2:3])
data.frame(estimate = coef(neMod3), LCL = confint[1, ], UCL = confint[2, ])

# function to calculate 95% confidence intervals 
# for linear combination of parameter estimates
linfunCI <- function(boot.out, L, conf) {
  est <- sum(L %*% boot.out$t0)
  se <- diag(sqrt(t(L) %*% var(boot.out$t) %*% L))
  bias <- est - L %*% colMeans(boot.out$t)
  CI <- est + bias + c(-1,1) * qnorm(1-(1-conf)/2) * se
  return(c("LCL" = CI[1], "UCL" = CI[2])) 
}

# eta3 + eta5
L <- c(0, 0, 0, 1, 0, 1, 0, 0)
c(L %*% bootSE$t0, linfunCI(bootSE, L, 0.95))

# test whether different three-way decompositions are significantly different (i.e. interaction terms differ from zero)

# function for Wald-type Chisquare test
# based on the bootstrap covariance matrix 
bootChisq <- function(boot.out, L) {
  chisq <- t(boot.out$t0) %*% t(L) %*% solve(L %*% var(boot.out$t)
                                             %*% t(L)) %*% L %*% boot.out$t0
  df <- dim(L)[1]
  p <- pchisq(q = chisq, df, lower.tail = FALSE)
  return(c("Chisq" = chisq, "df" = df, "p" = p))
}

L <- matrix(c(0, 0, 0, 0, 1, 0, 0, 0,
              0, 0, 0, 0, 0, 1, 0, 0,
              0, 0, 0, 0, 0, 0, 1, 0,
              0, 0, 0, 0, 0, 0, 0, 1),
            nrow = 4, byrow = TRUE)
bootChisq(bootSE, L)
# no evidence for significant differences => fit more parsimonious model that simplifies interpretation
neMod3 <- glm(immigr ~ treat0 + treat1 + treat2, data = expData, weights = w)
neMod3