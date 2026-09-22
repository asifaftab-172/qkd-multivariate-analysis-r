# M5 Assignment 5.1 - Multivariate Analysis of QKD Coexistence Data
# Student: Asif Akhtab Ronggon
# Dataset: Horgan et al. (2026), Zenodo DOI: 10.5281/zenodo.21792844
#
# This script downloads five measurement-level CSV files from the 0-km
# coexistence campaign, preprocesses the data, and performs EDA, PCA,
# k-means clustering, hierarchical clustering, and LDA.

options(stringsAsFactors = FALSE)
set.seed(2026)

# -----------------------------
# 0. Project directories
# -----------------------------
dirs <- c("data/raw", "data/processed", "figures", "results")
for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

base_url <- "https://zenodo.org/records/21792844/files/"
files <- c(
  "baseline_no_signal_no_added_fibre_full_sweep.csv",
  "roadm_3dbm_injected_signal_no_added_fibre_full_sweep.csv",
  "roadm_7dbm_injected_signal_no_added_fibre_full_sweep.csv",
  "roadm_9dbm_injected_signal_no_added_fibre_full_sweep.csv",
  "roadm_12dbm_injected_signal_no_added_fibre_full_sweep.csv"
)

# -----------------------------
# 1. Data acquisition and preprocessing
# -----------------------------
for (f in files) {
  dest <- file.path("data/raw", f)
  if (!file.exists(dest)) {
    message("Downloading ", f)
    download.file(paste0(base_url, f, "?download=1"), destfile = dest,
                  mode = "wb", quiet = FALSE)
  }
}

raw_list <- lapply(file.path("data/raw", files), read.csv, check.names = FALSE)
dat <- do.call(rbind, raw_list)
rownames(dat) <- NULL

# Validate the analysis schema before accessing or transforming columns.
required_columns <- c("Sweep_km", "Noise_dBm", "VOA_dB", "QBER",
                      "SecureKeyRate_bps")
missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns) > 0) {
  stop("Missing required source columns: ", paste(missing_columns, collapse = ", "))
}
if (nrow(dat) == 0) stop("The source files contain no observations.")
numeric_columns <- c("Sweep_km", "VOA_dB", "QBER", "SecureKeyRate_bps")
if (!all(vapply(dat[numeric_columns], is.numeric, logical(1)))) {
  stop("Sweep_km, VOA_dB, QBER, and SecureKeyRate_bps must be numeric.")
}
if (!all(vapply(dat[numeric_columns],
                function(x) all(is.finite(x)), logical(1)))) {
  stop("Missing or non-finite values detected in required numeric columns.")
}
if (any(dat$QBER < 0 | dat$QBER > 1) ||
    any(dat$SecureKeyRate_bps < 0)) {
  stop("QBER must lie in [0, 1] and secure key rate must be nonnegative.")
}
if (anyNA(dat$Noise_dBm) ||
    any(!as.character(dat$Noise_dBm) %in% c("none", "3", "7", "9", "12"))) {
  stop("Noise_dBm must contain only none, 3, 7, 9, or 12.")
}

# Keep only the back-to-back (0 km) measurements. The selected source files
# already satisfy this condition; the assertion protects reproducibility.
stopifnot(all(dat$Sweep_km == 0))

# Create readable categorical condition labels.
dat$NoiseCondition <- ifelse(dat$Noise_dBm == "none", "Baseline",
                             paste0(dat$Noise_dBm, " dBm"))
dat$NoiseCondition <- factor(dat$NoiseCondition,
  levels = c("Baseline", "3 dBm", "7 dBm", "9 dBm", "12 dBm"))
if (anyNA(dat$NoiseCondition) || any(table(dat$NoiseCondition) < 3)) {
  stop("Each of the five recognized conditions needs at least three observations.")
}

# Secure key rate spans orders of magnitude and contains zeros near collapse.
# log10(SKR + 1) preserves zero-SKR measurements while reducing skew.
dat$LogSKR <- log10(dat$SecureKeyRate_bps + 1)

# Count missing values across all source and derived fields. Required analysis
# fields have already been checked; unused metadata may still contain NA.
missing_total <- sum(is.na(dat))

write.csv(dat, "data/processed/QKD_Coexistence_0km_Combined.csv", row.names = FALSE)

# Quantitative variables used for multivariate analysis.
X <- dat[, c("VOA_dB", "QBER", "LogSKR")]
predictor_sd <- vapply(X, sd, numeric(1))
if (any(!is.finite(predictor_sd) | predictor_sd <= 0)) {
  stop("Every analysis predictor must have a finite, positive standard deviation.")
}
X_scaled <- scale(X)

# Descriptive statistics helper.
desc_fun <- function(x) c(
  N = length(x), Mean = mean(x), SD = sd(x), Median = median(x),
  Min = min(x), Max = max(x)
)

overall_desc <- do.call(rbind, lapply(X, desc_fun))
overall_desc <- data.frame(Variable = rownames(overall_desc), overall_desc,
                           row.names = NULL, check.names = FALSE)
write.csv(overall_desc, "results/descriptive_statistics_overall.csv", row.names = FALSE)

cond_levels <- levels(dat$NoiseCondition)
cond_desc <- do.call(rbind, lapply(cond_levels, function(g) {
  d <- dat[dat$NoiseCondition == g, ]
  data.frame(
    NoiseCondition = g,
    N = nrow(d),
    Mean_VOA_dB = mean(d$VOA_dB), SD_VOA_dB = sd(d$VOA_dB),
    Mean_QBER = mean(d$QBER), SD_QBER = sd(d$QBER),
    Median_SKR_bps = median(d$SecureKeyRate_bps),
    Mean_LogSKR = mean(d$LogSKR), SD_LogSKR = sd(d$LogSKR),
    Zero_SKR_pct = 100 * mean(d$SecureKeyRate_bps == 0)
  )
}))
write.csv(cond_desc, "results/descriptive_statistics_by_condition.csv", row.names = FALSE)

cor_mat <- cor(X)
write.csv(cor_mat, "results/correlation_matrix.csv")

# -----------------------------
# 2. EDA figures
# -----------------------------
condition_cols <- c("#222222", "#0072B2", "#009E73", "#E69F00", "#D55E00")
condition_pch <- c(16, 17, 15, 18, 8)
cond_idx <- as.integer(dat$NoiseCondition)

pdf("figures/01_scatterplot_matrix.pdf", width = 7.2, height = 7.2)
pairs(X,
      col = adjustcolor(condition_cols[cond_idx], alpha.f = 0.35),
      pch = 16, cex = 0.45,
      oma = c(6, 4, 4, 2),
      main = "Scatterplot Matrix of QKD Performance Variables")
# Draw the key in its own strip, outside the scatterplot panels.
par(fig = c(0, 1, 0, 0.12), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend("center", legend = cond_levels, col = condition_cols, pch = 16,
       cex = 0.85, bty = "n", horiz = TRUE)
dev.off()

pdf("figures/02_correlation_matrix.pdf", width = 6.3, height = 5.5)
par(mar = c(5, 5, 3, 2))
image(1:ncol(cor_mat), 1:ncol(cor_mat), t(cor_mat[nrow(cor_mat):1, ]),
      axes = FALSE, xlab = "", ylab = "",
      col = hcl.colors(50, "Blue-Red 3"), zlim = c(-1, 1),
      main = "Correlation Matrix")
axis(1, at = 1:ncol(cor_mat), labels = colnames(cor_mat), las = 2)
axis(2, at = 1:ncol(cor_mat), labels = rev(rownames(cor_mat)), las = 2)
for (i in seq_len(nrow(cor_mat))) {
  for (j in seq_len(ncol(cor_mat))) {
    text(j, nrow(cor_mat) - i + 1, sprintf("%.2f", cor_mat[i, j]),
         cex = 0.95, col = if (abs(cor_mat[i, j]) > 0.6) "white" else "black")
  }
}
box()
dev.off()

# Condition-wise trends with semi-transparent raw points and mean lines.
pdf("figures/03_qber_vs_attenuation.pdf", width = 7.2, height = 5.2)
plot(dat$VOA_dB, dat$QBER, type = "n", xlab = "Optical attenuation (VOA, dB)",
     ylab = "Quantum Bit Error Rate (QBER)",
     main = "QBER Across Attenuation and Noise Conditions")
for (i in seq_along(cond_levels)) {
  g <- cond_levels[i]; d <- dat[dat$NoiseCondition == g, ]
  points(d$VOA_dB, d$QBER, pch = 16, cex = 0.4,
         col = adjustcolor(condition_cols[i], alpha.f = 0.22))
  m <- aggregate(QBER ~ VOA_dB, d, mean)
  lines(m$VOA_dB, m$QBER, col = condition_cols[i], lwd = 2)
}
legend("topleft", legend = cond_levels, col = condition_cols, lwd = 2,
       cex = 0.8, bty = "n")
dev.off()

pdf("figures/04_logskr_vs_attenuation.pdf", width = 7.2, height = 5.2)
plot(dat$VOA_dB, dat$LogSKR, type = "n", xlab = "Optical attenuation (VOA, dB)",
     ylab = expression(log[10](SKR + 1)),
     main = "Secure Key Rate Across Attenuation and Noise Conditions")
for (i in seq_along(cond_levels)) {
  g <- cond_levels[i]; d <- dat[dat$NoiseCondition == g, ]
  points(d$VOA_dB, d$LogSKR, pch = 16, cex = 0.4,
         col = adjustcolor(condition_cols[i], alpha.f = 0.22))
  m <- aggregate(LogSKR ~ VOA_dB, d, mean)
  lines(m$VOA_dB, m$LogSKR, col = condition_cols[i], lwd = 2)
}
legend("bottomleft", legend = cond_levels, col = condition_cols, lwd = 2,
       cex = 0.8, bty = "n")
dev.off()

# -----------------------------
# 3. Principal Component Analysis
# -----------------------------
pca <- prcomp(X, center = TRUE, scale. = TRUE)
pca_var <- (pca$sdev^2) / sum(pca$sdev^2)
pca_cum <- cumsum(pca_var)
n_pc_80 <- which(pca_cum >= 0.80)[1]

pca_variance <- data.frame(
  Component = paste0("PC", seq_along(pca_var)),
  Eigenvalue = pca$sdev^2,
  VarianceExplained = pca_var,
  CumulativeVariance = pca_cum
)
write.csv(pca_variance, "results/pca_variance.csv", row.names = FALSE)
write.csv(pca$rotation, "results/pca_loadings.csv")

pdf("figures/05_pca_scree.pdf", width = 6.8, height = 4.8)
scree_x <- barplot(100 * pca_var, names.arg = paste0("PC", seq_along(pca_var)),
        ylab = "Variance explained (%)", xlab = "Principal component",
        main = "PCA Scree Plot", ylim = c(0, max(100 * pca_var) * 1.18))
lines(scree_x, 100 * pca_var, type = "b", pch = 16, lwd = 1.5)
abline(h = 100 / ncol(X), lty = 2)
dev.off()

scores <- as.data.frame(pca$x)
scores$NoiseCondition <- dat$NoiseCondition

pdf("figures/06_pca_biplot.pdf", width = 7.2, height = 5.8)
load <- pca$rotation[, 1:2, drop = FALSE]
arrow_scale <- 0.78 * min(diff(range(scores$PC1)) / diff(range(load[,1])),
                          diff(range(scores$PC2)) / diff(range(load[,2])))
# Include loading arrows and their labels as well as observation scores.
biplot_xlim <- extendrange(c(scores$PC1, load[,1] * arrow_scale * 1.18))
biplot_ylim <- extendrange(c(scores$PC2, load[,2] * arrow_scale * 1.18))
plot(scores$PC1, scores$PC2, type = "n",
     xlim = biplot_xlim, ylim = biplot_ylim,
     xlab = sprintf("PC1 (%.1f%%)", 100 * pca_var[1]),
     ylab = sprintf("PC2 (%.1f%%)", 100 * pca_var[2]),
     main = "PCA Biplot of QKD Coexistence Measurements")
for (i in seq_along(cond_levels)) {
  sel <- dat$NoiseCondition == cond_levels[i]
  points(scores$PC1[sel], scores$PC2[sel], pch = condition_pch[i],
         col = adjustcolor(condition_cols[i], alpha.f = 0.45), cex = 0.65)
}
# Loading arrows, scaled to score space.
arrows(0, 0, load[,1] * arrow_scale, load[,2] * arrow_scale,
       length = 0.08, lwd = 1.6)
text(load[,1] * arrow_scale * 1.08, load[,2] * arrow_scale * 1.08,
     labels = rownames(load), cex = 0.9, font = 2)
legend("topright", legend = cond_levels, col = condition_cols,
       pch = condition_pch, cex = 0.78, bty = "n")
abline(h = 0, v = 0, col = "grey80", lty = 3)
dev.off()

# -----------------------------
# 4. Cluster analysis
# -----------------------------
# K-means is applied to standardized original variables.
k_range <- 1:10
wss <- sapply(k_range, function(k) {
  kmeans(X_scaled, centers = k, nstart = 50, iter.max = 100)$tot.withinss
})

# Automated elbow estimate: point farthest from the straight line joining k=1 and k=10.
x <- k_range; y <- wss
x1 <- x[1]; y1 <- y[1]; x2 <- x[length(x)]; y2 <- y[length(y)]
dist_line <- abs((y2-y1)*x - (x2-x1)*y + x2*y1 - y2*x1) /
             sqrt((y2-y1)^2 + (x2-x1)^2)
optimal_k <- x[which.max(dist_line)]
optimal_k <- max(2, optimal_k)

write.csv(data.frame(k = k_range, WSS = wss, DistanceFromChord = dist_line),
          "results/kmeans_elbow.csv", row.names = FALSE)

pdf("figures/07_elbow_method.pdf", width = 6.6, height = 4.8)
plot(k_range, wss, type = "b", pch = 16, lwd = 2,
     xlab = "Number of clusters (k)", ylab = "Total within-cluster sum of squares",
     main = "Elbow Method for K-means")
abline(v = optimal_k, lty = 2, lwd = 1.5)
text(optimal_k, max(wss) * 0.9, labels = paste("Selected k =", optimal_k), pos = 4)
dev.off()

km <- kmeans(X_scaled, centers = optimal_k, nstart = 100, iter.max = 200)
write.csv(data.frame(Observation = seq_len(nrow(dat)), Cluster = km$cluster,
                     NoiseCondition = dat$NoiseCondition),
          "results/kmeans_membership.csv", row.names = FALSE)
km_cross <- table(Cluster = km$cluster, NoiseCondition = dat$NoiseCondition)
write.csv(as.data.frame.matrix(km_cross), "results/kmeans_vs_condition.csv")

pdf("figures/08_kmeans_clusters_pca.pdf", width = 7.0, height = 5.3)
plot(scores$PC1, scores$PC2, col = km$cluster, pch = 16, cex = 0.62,
     xlab = sprintf("PC1 (%.1f%%)", 100 * pca_var[1]),
     ylab = sprintf("PC2 (%.1f%%)", 100 * pca_var[2]),
     main = sprintf("K-means Clusters in PCA Space (k = %d)", optimal_k))
legend("topright", legend = paste("Cluster", seq_len(optimal_k)),
       col = seq_len(optimal_k), pch = 16, cex = 0.8, bty = "n")
dev.off()

# Hierarchical clustering with Ward's minimum-variance method.
dmat <- dist(X_scaled, method = "euclidean")
hc <- hclust(dmat, method = "ward.D2")
hc_cluster <- cutree(hc, k = optimal_k)
hc_cross <- table(Cluster = hc_cluster, NoiseCondition = dat$NoiseCondition)
write.csv(as.data.frame.matrix(hc_cross), "results/hierarchical_vs_condition.csv")

pdf("figures/09_hierarchical_dendrogram.pdf", width = 8.4, height = 5.4)
plot(hc, labels = FALSE, hang = -1, main = "Hierarchical Clustering Dendrogram",
     xlab = "QKD observations", ylab = "Ward linkage height")
rect.hclust(hc, k = optimal_k, border = 2:(optimal_k + 1))
dev.off()

pdf("figures/10_hierarchical_clusters_pca.pdf", width = 7.0, height = 5.3)
plot(scores$PC1, scores$PC2, col = hc_cluster, pch = 16, cex = 0.62,
     xlab = sprintf("PC1 (%.1f%%)", 100 * pca_var[1]),
     ylab = sprintf("PC2 (%.1f%%)", 100 * pca_var[2]),
     main = sprintf("Hierarchical Clusters in PCA Space (k = %d)", optimal_k))
legend("topright", legend = paste("Cluster", seq_len(optimal_k)),
       col = seq_len(optimal_k), pch = 16, cex = 0.8, bty = "n")
dev.off()

# -----------------------------
# 5. Linear Discriminant Analysis
# -----------------------------
if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package 'MASS' is required. It is included with standard R distributions.")
}

# Stratified 70/30 split evaluates held-out rows from the same campaign.
# Shared experimental settings may occur in both sets, so this does not
# estimate generalization to independent runs or previously unseen settings.
train_idx <- unlist(lapply(cond_levels, function(g) {
  ids <- which(dat$NoiseCondition == g)
  sample(ids, size = floor(0.70 * length(ids)))
}))
train_idx <- sort(train_idx)
test_idx <- setdiff(seq_len(nrow(dat)), train_idx)

train <- dat[train_idx, c("NoiseCondition", "VOA_dB", "QBER", "LogSKR")]
test  <- dat[test_idx,  c("NoiseCondition", "VOA_dB", "QBER", "LogSKR")]

# Standardize predictors using training-set parameters to put coefficients
# on a common scale without using test-set information for preprocessing.
mu <- sapply(train[, -1], mean)
sig <- sapply(train[, -1], sd)
if (any(!is.finite(sig) | sig <= 0)) {
  stop("Every training predictor must have a finite, positive standard deviation.")
}
standardize_frame <- function(d) {
  out <- data.frame(NoiseCondition = d$NoiseCondition)
  for (v in names(mu)) out[[paste0(v, "_z")]] <- (d[[v]] - mu[v]) / sig[v]
  out
}
train_z <- standardize_frame(train)
test_z <- standardize_frame(test)
all_z <- standardize_frame(dat[, c("NoiseCondition", "VOA_dB", "QBER", "LogSKR")])

lda_fit <- MASS::lda(NoiseCondition ~ VOA_dB_z + QBER_z + LogSKR_z, data = train_z)
lda_scaling <- lda_fit$scaling
write.csv(lda_scaling, "results/lda_coefficients.csv")
write.csv(lda_fit$means, "results/lda_group_means_standardized.csv")

lda_test_pred <- predict(lda_fit, newdata = test_z)
conf <- table(Actual = test_z$NoiseCondition, Predicted = lda_test_pred$class)
accuracy <- sum(diag(conf)) / sum(conf)
write.csv(as.data.frame.matrix(conf), "results/lda_confusion_matrix.csv")
write.csv(data.frame(Metric = c("Held-out accuracy", "Training observations", "Test observations"),
                     Value = c(accuracy, nrow(train_z), nrow(test_z))),
          "results/lda_performance.csv", row.names = FALSE)

lda_all <- predict(lda_fit, newdata = all_z)
lda_scores <- as.data.frame(lda_all$x)
lda_scores$NoiseCondition <- dat$NoiseCondition
write.csv(lda_scores, "results/lda_scores.csv", row.names = FALSE)

pdf("figures/11_lda_separation.pdf", width = 7.0, height = 5.4)
if (ncol(lda_all$x) >= 2) {
  plot(lda_all$x[,1], lda_all$x[,2], type = "n",
       xlab = "LD1", ylab = "LD2",
       main = "LDA Separation of Noise Conditions")
  for (i in seq_along(cond_levels)) {
    sel <- dat$NoiseCondition == cond_levels[i]
    points(lda_all$x[sel,1], lda_all$x[sel,2], pch = condition_pch[i],
           col = adjustcolor(condition_cols[i], alpha.f = 0.45), cex = 0.65)
  }
  legend("topright", legend = cond_levels, col = condition_cols,
         pch = condition_pch, cex = 0.78, bty = "n")
} else {
  plot(lda_all$x[,1], jitter(rep(0, nrow(lda_all$x))), pch = 16,
       col = condition_cols[cond_idx], xlab = "LD1", ylab = "",
       yaxt = "n", main = "LDA Discriminant Scores")
}
dev.off()

pdf("figures/12_lda_confusion_matrix.pdf", width = 6.4, height = 5.6)
conf_num <- unclass(conf)
image(seq_len(ncol(conf_num)), seq_len(nrow(conf_num)), t(conf_num[nrow(conf_num):1, ]),
      axes = FALSE, xlab = "Predicted condition", ylab = "Actual condition",
      col = hcl.colors(30, "Blues 3"), main = "LDA Confusion Matrix (Held-out Test Set)")
axis(1, at = seq_len(ncol(conf_num)), labels = colnames(conf_num), las = 2, cex.axis = 0.8)
axis(2, at = seq_len(nrow(conf_num)), labels = rev(rownames(conf_num)), las = 2, cex.axis = 0.8)
for (i in seq_len(nrow(conf_num))) {
  for (j in seq_len(ncol(conf_num))) {
    text(j, nrow(conf_num) - i + 1, conf_num[i,j], font = 2)
  }
}
box()
dev.off()

# -----------------------------
# 6. Machine-readable and LaTeX-ready summaries
# -----------------------------
# Dominant variables for concise interpretation.
pc1_order <- names(sort(abs(pca$rotation[,1]), decreasing = TRUE))
pc2_order <- names(sort(abs(pca$rotation[,2]), decreasing = TRUE))
ld1_order <- rownames(lda_scaling)[order(abs(lda_scaling[,1]), decreasing = TRUE)]
ld2_order <- if (ncol(lda_scaling) >= 2) rownames(lda_scaling)[order(abs(lda_scaling[,2]), decreasing = TRUE)] else NA

summary_lines <- c(
  sprintf("Observations: %d", nrow(dat)),
  sprintf("Missing values: %d", missing_total),
  sprintf("PC1 variance explained: %.2f%%", 100*pca_var[1]),
  sprintf("PC2 variance explained: %.2f%%", 100*pca_var[2]),
  sprintf("PC1+PC2 cumulative variance: %.2f%%", 100*pca_cum[2]),
  sprintf("PCs needed for >=80%% variance: %d", n_pc_80),
  sprintf("Largest absolute PC1 loading: %s", pc1_order[1]),
  sprintf("Largest absolute PC2 loading: %s", pc2_order[1]),
  sprintf("Selected k by elbow geometry: %d", optimal_k),
  sprintf("LDA held-out accuracy: %.2f%%", 100*accuracy),
  sprintf("Largest absolute LD1 coefficient: %s", ld1_order[1]),
  sprintf("Largest absolute LD2 coefficient: %s", ifelse(is.na(ld2_order[1]), "N/A", ld2_order[1]))
)
writeLines(summary_lines, "results/analysis_summary.txt")

# Macros consumed by the LaTeX report.
escape_tex <- function(x) gsub("_", "\\_", x, fixed = TRUE)
macros <- c(
  sprintf("\\newcommand{\\NObs}{%d}", nrow(dat)),
  sprintf("\\newcommand{\\NMissing}{%d}", missing_total),
  sprintf("\\newcommand{\\PCOneVar}{%.1f\\%%}", 100*pca_var[1]),
  sprintf("\\newcommand{\\PCTwoVar}{%.1f\\%%}", 100*pca_var[2]),
  sprintf("\\newcommand{\\PCTwoCum}{%.1f\\%%}", 100*pca_cum[2]),
  sprintf("\\newcommand{\\NPCOptimal}{%d}", n_pc_80),
  sprintf("\\newcommand{\\PCOneTop}{%s}", escape_tex(pc1_order[1])),
  sprintf("\\newcommand{\\PCTwoTop}{%s}", escape_tex(pc2_order[1])),
  sprintf("\\newcommand{\\KOptimal}{%d}", optimal_k),
  sprintf("\\newcommand{\\LDAAccuracy}{%.1f\\%%}", 100*accuracy),
  sprintf("\\newcommand{\\TrainN}{%d}", nrow(train_z)),
  sprintf("\\newcommand{\\TestN}{%d}", nrow(test_z)),
  sprintf("\\newcommand{\\LDOneTop}{%s}", escape_tex(sub("_z$", "", ld1_order[1]))),
  sprintf("\\newcommand{\\LDTwoTop}{%s}", escape_tex(sub("_z$", "", ifelse(is.na(ld2_order[1]), "N/A", ld2_order[1]))))
)
writeLines(macros, "results/results.tex")

# Compact LaTeX tables generated from results.
fmt <- function(x, d=3) formatC(x, format = "f", digits = d)

# Descriptive table by condition.
con <- file("results/descriptive_table.tex", "w")
writeLines("\\begin{tabular}{lrrrr}", con)
writeLines("\\toprule", con)
writeLines("Condition & $n$ & Mean QBER & Median SKR (bps) & Zero-SKR (\\%) \\\\", con)
writeLines("\\midrule", con)
for (i in seq_len(nrow(cond_desc))) {
  r <- cond_desc[i,]
  writeLines(sprintf("%s & %d & %s & %s & %s \\\\",
    r$NoiseCondition, r$N, fmt(r$Mean_QBER,4),
    format(round(r$Median_SKR_bps), big.mark=",", scientific=FALSE), fmt(r$Zero_SKR_pct,1)), con)
}
writeLines("\\bottomrule", con); writeLines("\\end{tabular}", con); close(con)

# PCA table.
con <- file("results/pca_table.tex", "w")
writeLines("\\begin{tabular}{lrrr}", con); writeLines("\\toprule", con)
writeLines("Component & Eigenvalue & Variance (\\%) & Cumulative (\\%) \\\\", con)
writeLines("\\midrule", con)
for (i in seq_len(nrow(pca_variance))) {
  r <- pca_variance[i,]
  writeLines(sprintf("%s & %s & %s & %s \\\\", r$Component,
                     fmt(r$Eigenvalue,3), fmt(100*r$VarianceExplained,1),
                     fmt(100*r$CumulativeVariance,1)), con)
}
writeLines("\\bottomrule", con); writeLines("\\end{tabular}", con); close(con)

# PCA loadings table.
con <- file("results/pca_loadings_table.tex", "w")
writeLines("\\begin{tabular}{lrrr}", con); writeLines("\\toprule", con)
writeLines("Variable & PC1 & PC2 & PC3 \\\\", con); writeLines("\\midrule", con)
for (v in rownames(pca$rotation)) {
  vals <- pca$rotation[v,]
  writeLines(sprintf("%s & %s & %s & %s \\\\", escape_tex(v),
                     fmt(vals[1],3), fmt(vals[2],3), fmt(vals[3],3)), con)
}
writeLines("\\bottomrule", con); writeLines("\\end{tabular}", con); close(con)

# LDA coefficient table.
con <- file("results/lda_coefficients_table.tex", "w")
writeLines(paste0("\\begin{tabular}{l", paste(rep("r", ncol(lda_scaling)), collapse=""), "}"), con)
writeLines("\\toprule", con)
writeLines(paste(c("Predictor", colnames(lda_scaling)), collapse = " & ") |> paste0(" \\\\"), con)
writeLines("\\midrule", con)
for (v in rownames(lda_scaling)) {
  vals <- sapply(lda_scaling[v,], fmt, d=3)
  writeLines(paste(c(escape_tex(sub("_z$", "", v)), vals), collapse=" & ") |> paste0(" \\\\"), con)
}
writeLines("\\bottomrule", con); writeLines("\\end{tabular}", con); close(con)

# Confusion matrix LaTeX table.
con <- file("results/lda_confusion_table.tex", "w")
writeLines(paste0("\\begin{tabular}{l", paste(rep("r", ncol(conf)), collapse=""), "}"), con)
writeLines("\\toprule", con)
writeLines(paste(c("Actual $\\backslash$ Pred.", colnames(conf)), collapse=" & ") |> paste0(" \\\\"), con)
writeLines("\\midrule", con)
for (i in seq_len(nrow(conf))) {
  writeLines(paste(c(rownames(conf)[i], conf[i,]), collapse=" & ") |> paste0(" \\\\"), con)
}
writeLines("\\bottomrule", con); writeLines("\\end{tabular}", con); close(con)

message("Analysis complete. Results written to data/processed, results, and figures.")
message(paste(summary_lines, collapse = "\n"))
