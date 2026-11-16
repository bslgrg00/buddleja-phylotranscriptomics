# 05_anomaly_zone.R: anomaly-zone internode pairs (Degnan and Rosenberg 2006) on the
#                    ASTRAL-III tree, branch lengths in coalescent units
# Output: anomaly_zone.csv (Appendix S8)
library(ape)

a_x <- function(x) log(2/3 + (3*exp(2*x) - 2) / (18*(exp(3*x) - exp(2*x))))

tr <- root(read.tree("astral_M20.tre"), resolve.root = TRUE,
           outgroup = c("Rehmannia_glutinosa_OWAS", "Lindenbergia_philippensis_WUZV"))
E <- tr$edge; L <- tr$edge.length; nt <- Ntip(tr)

# every pair of adjacent internal branches: parent (x) and descendant (y)
res <- do.call(rbind, lapply(which(E[, 2] > nt), function(i) {
  j <- which(E[, 1] == E[i, 2] & E[, 2] > nt)
  if (length(j) == 0) return(NULL)
  data.frame(parent_node = E[i, 1], child_node = E[i, 2], grandchild_node = E[j, 2],
             x_parent = L[i], y_descendant = L[j], a_of_x = a_x(L[i]),
             clade_below_y = sapply(E[j, 2], function(n)
               paste(sub(".*_", "", extract.clade(tr, n)$tip.label), collapse = " ")))
}))
res$margin <- res$a_of_x - res$y_descendant
res$in_anomaly_zone <- res$margin > 0
num <- c("x_parent", "y_descendant", "a_of_x", "margin")
res[num] <- lapply(res[num], round, 5)
write.csv(res[order(res$x_parent), ], "anomaly_zone.csv", row.names = FALSE)
cat(sum(res$in_anomaly_zone), "of", nrow(res), "internode pairs in the anomaly zone\n")
