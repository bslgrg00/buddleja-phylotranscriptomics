# 05_MSCquartets.R: quartet tests under the multispecies coalescent (T1, ILS) and
#                   NANUQ distances (T3, gene flow), all 34 taxa
# Input: M20_genetrees_rooted.tre (04_concordance.sh), astral_M20.tre
# Output: tree_rejection_by_alpha.csv, simplex_T1_alpha*.pdf (Figure 2B),
#         nanuq_all34_alpha1e-04_beta0.95.nex (Neighbor-Net in SplitsTree, Figure 3A)
library(MSCquartets); library(ape)

genetrees <- read.tree("M20_genetrees_rooted.tre")
taxa      <- sort(unique(unlist(lapply(genetrees, function(t) t$tip.label))))
sp        <- read.tree("astral_M20.tre"); sp$node.label <- NULL; sp$edge.length <- NULL
ALPHAS    <- c(0.05, 0.01, 1e-3, 1e-4, 1e-6)

# Quartet counts; unresolved quartets split equally among the three resolutions
RQT <- quartetTableResolved(quartetTable(genetrees, taxonnames = taxa), omit = FALSE)

# 1. T1 test on the ASTRAL-III topology and T3 test
pT1 <- quartetStarTestInd(quartetTreeTestInd(RQT, model = "T1",
                                             speciestree = write.tree(unroot(sp))))
pT3 <- quartetStarTestInd(quartetTreeTestInd(RQT, model = "T3"))

# 2. Quartets rejecting each model, raw (reported) and Holm-corrected
rates <- function(pT, m) do.call(rbind, lapply(ALPHAS, function(a) {
  n <- nrow(pT); hb <- HolmBonferroni(pT, m, a)
  nh <- sum(hb[[ncol(hb)]] == "Y"); nr <- sum(pT[, paste0("p_", m)] < a)
  data.frame(set = "all34", model = m, n_quartets = n, alpha = a,
             n_reject_holm = nh, pct_reject_holm = round(100 * nh / n, 2),
             n_reject_raw = nr, pct_reject_raw = round(100 * nr / n, 2))
}))
write.csv(rbind(rates(pT1, "T1"), rates(pT3, "T3")), "tree_rejection_by_alpha.csv",
          row.names = FALSE)

# 3. Simplex plots of the T1 test (Figure 2B)
for (a in c(0.01, 1e-3, 1e-4, 1e-6)) {
  pdf(sprintf("simplex_T1_alpha%g.pdf", a), 7, 7)
  quartetTestPlot(pT1, test = "T1", alpha = a, beta = 0.95)
  dev.off()
}

# 4. NANUQ distances (alpha = 1e-4, beta = 0.95)
NANUQdist(pT3, outfile = "nanuq_all34", alpha = 1e-4, beta = 0.95, plot = FALSE)
