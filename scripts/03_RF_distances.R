# 03_RF_distances.R: Robinson-Foulds distances among the 12 species trees
# Output: RF_normalized.csv, RF_raw_counts.csv, RF_distance_heatmap.pdf/.png
library(ape); library(phangorn); library(ggplot2); library(reshape2)

tiers <- c(10, 20, 30)
files <- c(sprintf("M%d_IQTREE.treefile", tiers), sprintf("M%d_RAxML.support", tiers),
           sprintf("astral_M%d.tre", tiers), sprintf("astral_pro3_M%d.tre", tiers))
labs  <- paste(rep(c("IQtree", "RAxML", "ASTRAL", "ASTRAL_pro3"), each = 3), tiers, sep = "_")

trees <- lapply(files, function(f) {
  t <- read.tree(f); t$tip.label <- gsub(" ", "_", t$tip.label); unroot(t) })
class(trees) <- "multiPhylo"
stopifnot(all(sapply(trees, function(t) setequal(t$tip.label, trees[[1]]$tip.label))))

nrf <- as.matrix(RF.dist(trees, normalize = TRUE));  dimnames(nrf) <- list(labs, labs)
rf  <- as.matrix(RF.dist(trees, normalize = FALSE)); dimnames(rf)  <- list(labs, labs)
off_mean <- function(m) { diag(m) <- NA; round(rowMeans(m, na.rm = TRUE), 4) }
write.csv(cbind(round(nrf, 4), Mean_nRF = off_mean(nrf)), "RF_normalized.csv")
write.csv(cbind(rf, Mean_distance = off_mean(rf)), "RF_raw_counts.csv")

d <- melt(nrf, varnames = c("Tree1", "Tree2"), value.name = "RF")
d$Tree1 <- factor(d$Tree1, levels = labs); d$Tree2 <- factor(d$Tree2, levels = rev(labs))
p <- ggplot(d, aes(Tree1, Tree2, fill = RF)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(RF == 0, "", sprintf("%.2f", RF))), size = 3, fontface = "bold") +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0.10,
                       limits = c(0, 1), name = "Normalized\nRF distance") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank(),
        panel.grid = element_blank())
ggsave("RF_distance_heatmap.pdf", p, width = 8, height = 7)
ggsave("RF_distance_heatmap.png", p, width = 8, height = 7, dpi = 600)
