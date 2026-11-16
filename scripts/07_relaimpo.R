# 07_relaimpo.R: relative importance of ILS, gene flow, and estimation error for gene tree
#                variation (LMG metric, relaimpo; 95% CIs from 100 bootstrap replicates)
# Input:  relative_contribution.csv (one row per node; Y_BSsum = 100 - gCF, X1_ILS = theta,
#         X2_Err = estimation error, X3_Hyb = reticulation index)
# Output: Buddleja_relative_contribution_with_CI.csv (Figure 4, Appendix S18)
library(relaimpo)

d   <- read.csv("relative_contribution.csv")
fit <- lm(Y_BSsum ~ X1_ILS + X2_Err + X3_Hyb, data = d)
cat("R2 =", round(100 * summary(fit)$r.squared, 2), "%\n")

set.seed(1)
ci <- booteval.relimp(boot.relimp(fit, b = 100, type = "lmg", rela = TRUE))
k  <- c(1, 3, 2)                       # ILS, gene flow, estimation error
write.csv(data.frame(Factor = c("ILS", "Gene flow", "Estimation error"),
                     Percentage = 100 * ci@lmg[k],
                     Lower = 100 * ci@lmg.lower[1, k], Upper = 100 * ci@lmg.upper[1, k]),
          "Buddleja_relative_contribution_with_CI.csv", row.names = FALSE)
