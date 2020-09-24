library(tmle)
library(SuperLearner)
library(gam)
library(polspline)
library(ggplot2)
# ----------------------------
# simulation functions
#----------------------------
# function to bound probabilities
boundx <- function(x=c(1, 2, 3), new_min=3, new_max=9){
  old_min <- min(x)
  old_max <- max(x)
  stretch_me <- (new_max-new_min)/(old_max-old_min)
  stretch_me*(x-old_min) + new_min
}

# data generating process
simple_simdat <- function(rseed, n, true_ate){
  set.seed(rseed)
  W1 <- rnorm(n)
  W2 <- rnorm(n)
  pA <- boundx(W1+W2+W1*W2, new_min=0.05, new_max=0.95)
  A <- rbinom(n, 1, pA)
  epsilon<- rnorm(n, 0, .2)
  Y <- true_ate*A - W1 - W2 - 2*W1*W2 + epsilon
  data.frame(W1=W1, W2=W2, pA=pA, A=A, Y=Y)
}

# randomized data generating process
simple_random_simdat <- function(rseed, n, true_ate){
  set.seed(rseed)
  W1 <- rnorm(n)
  W2 <- rnorm(n)
  pA <- sample(c(0,1), size = n, replace = TRUE)
  A <- rbinom(n, 1, pA)
  epsilon<- rnorm(n, 0, .2)
  Y <- true_ate*A - W1 - W2 - 2*W1*W2 + epsilon
  data.frame(W1=W1, W2=W2, pA=pA, A=A, Y=Y)
}
#---------------------
# observational
#---------------------
dat <- simple_simdat(rseed=2020, n=5000, true_ate=0.25)
mean_diff = mean(dat$Y[dat$A==1]) - mean(dat$Y[dat$A==0])
ttest_out <- t.test(x=dat$Y[dat$A==1], y=dat$Y[dat$A==0])


glm_out <- glm(Y~A+W1+W2, data=dat)
glm_coef = summary(glm_out)$coefficients
glm_ci = confint(glm_out)


# glm works if you specify model correctly
glm_out_correct <- glm(Y~ A+W1+W2 + W1*W2, data=dat)
glm_coef_correct = summary(glm_out_correct)$coefficients
glm_ci_correct = confint(glm_out_correct)


tml_incorrect_exposure_model <-
  tmle(Y=dat$Y
       ,A=dat$A
       ,W=dat[, c("W1", "W2")]
       ,Q.SL.library=c("SL.glm", "SL.glm.interaction", "SL.gam", "SL.polymars")
       ,gform="A~W1+W2"    #note: this is incorrect because A  =  W1 + W2 + W1*W2
       ,V=10
  )
tmle_incorrect_exposure_model_ATE = tml_incorrect_exposure_model$estimates$ATE$psi
tmle_incorrect_exposure_model_CI = tml_incorrect_exposure_model$estimates$ATE$CI



tml_incorrect_outcome_model <- tmle(
  Y=dat$Y
  ,A=dat$A
  ,W=dat[, c("W1", "W2")]
  ,g.SL.library=c("SL.glm", "SL.glm.interaction", "SL.gam", "SL.polymars")
  ,Qform="Y~W1+W2"
  ,V=10
)
tmle_incorrect_outcome_model_ATE = tml_incorrect_outcome_model$estimates$ATE$psi
tmle_incorrect_outcome_model_CI = tml_incorrect_outcome_model$estimates$ATE$CI

tmle_learn_both_models <- tmle(Y=dat$Y
                               ,A=dat$A
                               ,W=dat[, c("W1", "W2")]
                               ,Q.SL.library=c("SL.glm", "SL.glm.interaction", "SL.gam", "SL.polymars")
                               ,g.SL.library=c("SL.glm", "SL.glm.interaction", "SL.gam", "SL.polymars")
                               ,V=10
)
tmle_learn_both_models_ATE = tmle_learn_both_models$estimates$ATE$psi
tmle_learn_both_models_CI = tmle_learn_both_models$estimates$ATE$CI


# plot
est_vec = c(mean_diff, 
            glm_coef[2,1], 
            glm_coef_correct[2,1],
            tmle_incorrect_exposure_model_ATE, 
            tmle_incorrect_outcome_model_ATE, 
            tmle_learn_both_models_ATE)

ci_matrix = rbind(ttest_out$conf.int, 
                  glm_ci[2,], 
                  glm_ci_correct[2,],
                  tmle_incorrect_exposure_model_CI,
                  tmle_incorrect_outcome_model_CI, 
                  tmle_learn_both_models_CI)
results = data.frame(est_vec, ci_matrix)
names(results) = c('estimate', 'ci_low', 'ci_high')
results$type = factor(c('mean_difference', 'glm_wrong_model', 'glm_correct_model', 
                        'tmle_wrong_expos_model',
                        'tmle_wrong_outcome_model',
                        'tmle_learn_both_models'),
                      levels = c('mean_difference', 'glm_wrong_model', 'glm_correct_model',
                        'tmle_wrong_expos_model',
                        'tmle_wrong_outcome_model',
                        'tmle_learn_both_models'))
ggplot(results, aes(type, est_vec)) + geom_point() + 
      geom_errorbar(aes(ymin =ci_low, ymax = ci_high )) + 
      theme(axis.text.x = element_text(angle = 25,  hjust = 1)) +
      ggtitle("Estimators attempting to recover true ATE = 0.25") +
      ylab('') + xlab ('') +
      geom_hline(yintercept = 0.25, color = 'red', linetype = 'dashed')

#-----------------
# random exposure
#-----------------
dat <- simple_random_simdat(rseed=2021, n=5000, true_ate=0.25)
mean_diff = mean(dat$Y[dat$A==1]) - mean(dat$Y[dat$A==0])
ttest_out <- t.test(x=dat$Y[dat$A==1], y=dat$Y[dat$A==0])

glm_out <- glm(Y~A+W1+W2, data=dat)
glm_coef = summary(glm_out)$coefficients
glm_ci = confint(glm_out)

# glm works if you specify model correctly
glm_out_correct <- glm(Y~A+W1+W2 + W1*W2, data=dat)
glm_coef_correct = summary(glm_out_correct)$coefficients
glm_ci_correct = confint(glm_out_correct)


tml_incorrect_exposure_model <-
  tmle(Y=dat$Y
       ,A=dat$A
       ,W=dat[, c("W1", "W2")]
       ,Q.SL.library=c("SL.glm", "SL.glm.interaction", "SL.gam", "SL.polymars")
       ,gform="A~W1+W2"    #note: this is incorrect because A  =  W1 + W2 + W1*W2
       ,V=10
  )
tmle_incorrect_exposure_model_ATE = tml_incorrect_exposure_model$estimates$ATE$psi
tmle_incorrect_exposure_model_CI = tml_incorrect_exposure_model$estimates$ATE$CI



tml_incorrect_outcome_model <- tmle(
  Y=dat$Y
  ,A=dat$A
  ,W=dat[, c("W1", "W2")]
  ,g.SL.library=c("SL.glm", "SL.glm.interaction", "SL.gam", "SL.polymars")
  ,Qform="Y~W1+W2"
  ,V=10
)
tmle_incorrect_outcome_model_ATE = tml_incorrect_outcome_model$estimates$ATE$psi
tmle_incorrect_outcome_model_CI = tml_incorrect_outcome_model$estimates$ATE$CI



tmle_learn_both_models <- tmle(Y=dat$Y
                               ,A=dat$A
                               ,W=dat[, c("W1", "W2")]
                               ,Q.SL.library=c("SL.glm", "SL.glm.interaction", "SL.gam", "SL.polymars")
                               ,g.SL.library=c("SL.glm", "SL.glm.interaction", "SL.gam", "SL.polymars")
                               ,V=10
)
tmle_learn_both_models_ATE = tmle_learn_both_models$estimates$ATE$psi
tmle_learn_both_models_CI = tmle_learn_both_models$estimates$ATE$CI


# plot
est_vec = c(mean_diff, 
            glm_coef[2,1], 
            glm_coef_correct[2,1],
            tmle_incorrect_exposure_model_ATE, 
            tmle_incorrect_outcome_model_ATE, 
            tmle_learn_both_models_ATE)

ci_matrix = rbind(ttest_out$conf.int, 
                  glm_ci[2,], 
                  glm_ci_correct[2,],
                  tmle_incorrect_exposure_model_CI,
                  tmle_incorrect_outcome_model_CI, 
                  tmle_learn_both_models_CI)
results = data.frame(est_vec, ci_matrix)
names(results) = c('estimate', 'ci_low', 'ci_high')
results$type = factor(c('mean_difference', 'glm_wrong_model', 'glm_correct_model', 
                        'tmle_wrong_expos_model',
                        'tmle_wrong_outcome_model',
                        'tmle_learn_both_models'),
                      levels = c('mean_difference', 'glm_wrong_model', 'glm_correct_model',
                                 'tmle_wrong_expos_model',
                                 'tmle_wrong_outcome_model',
                                 'tmle_learn_both_models'))
ggplot(results, aes(type, est_vec)) + geom_point() + 
  geom_errorbar(aes(ymin =ci_low, ymax = ci_high )) + 
  theme(axis.text.x = element_text(angle = 25,  hjust = 1)) +
  ggtitle("Estimators attempting to recover true ATE = 0.25") +
  ylab('') + xlab ('') +
  geom_hline(yintercept = 0.25, color = 'red', linetype = 'dashed')
