# M5 Assignment 5.1 - Multivariate Analysis of QKD Coexistence Data
# Student: Asif Akhtab Ronggon
# Source: Horgan et al. (2026), https://doi.org/10.5281/zenodo.21792844
#
# Run from the repository root, or place this script and the submission CSV
# together and run Rscript "M5Assign 5.1 Data__Ronggon_Asif.R".
# All analyses describe the 0-km campaign; measurements are not independent runs.
# Dependencies are base/recommended R packages MASS and cluster.

options(stringsAsFactors = FALSE)
for (pkg in c("MASS", "cluster")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Required recommended R package is missing: ", pkg)
}
for (d in c("data/raw", "data/processed", "figures", "results"))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
seeds <- c(ElbowBase = 202601L, FinalKmeans = 202602L, LDASplit = 202603L)
write.csv(data.frame(Stage = names(seeds), Seed = unname(seeds)),
          "results/random_seeds.csv", row.names = FALSE)

# Shared formatting. Text is escaped before insertion into LaTeX.
escape_tex <- function(x) {
  x <- gsub("_", "\\_", x, fixed = TRUE)
  x <- gsub("%", "\\%", x, fixed = TRUE)
  x <- gsub("&", "\\&", x, fixed = TRUE)
  x
}
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
pct <- function(x, d = 1) paste0(fmt(100 * x, d), "\\%")
tex_p <- function(x) ifelse(is.na(x), "--",
                           ifelse(x < 0.001, "$<0.001$", fmt(x, 3)))
tex_table <- function(path, headers, rows, align = NULL) {
  rows <- as.matrix(rows)
  if (is.null(align)) align <- paste0("l", paste(rep("r", ncol(rows) - 1), collapse = ""))
  lines <- c(paste0("\\begin{tabular}{", align, "}"), "\\toprule",
             paste0(paste(headers, collapse = " & "), " \\\\"), "\\midrule")
  for (i in seq_len(nrow(rows)))
    lines <- c(lines, paste0(paste(rows[i, ], collapse = " & "), " \\\\"))
  writeLines(c(lines, "\\bottomrule", "\\end{tabular}"), file.path("results", path))
}
write_findings <- function(name, paragraphs)
  writeLines(paste(paragraphs, collapse = "\n\n"), file.path("results", paste0("findings_", name, ".tex")))

ink <- "#20334A"
muted <- "#65748B"
accent <- "#067A83"
condition_cols <- c("#243B65", "#008C95", "#3C8DCC", "#D78A25", "#A94374")
condition_pch <- c(16, 17, 15, 18, 8)
variable_labels <- c(VOA_dB = "Attenuation", QBER = "QBER", LogSKR = "log10(SKR + 1)")
open_plot <- function(name, width = 8.3, height = 5.4, legend = FALSE) {
  pdf(file.path("figures", name), width = width, height = height,
      family = "Helvetica", useDingbats = FALSE)
  if (legend) layout(matrix(c(1, 2), ncol = 1), heights = c(1, 0.12))
  par(mar = c(4.7, 4.8, 3.7, 1.1), mgp = c(2.8, 0.7, 0), tcl = -0.25,
      col = ink, col.axis = muted, col.lab = ink, col.main = ink,
      cex = 0.95, cex.main = 1.15, font.main = 2, las = 1, bty = "l")
}
plot_grid <- function() grid(col = "#E7EBF0", lty = 1)
legend_strip <- function(labels, cols, pch = 16) {
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = labels, col = cols, pch = pch,
         horiz = TRUE, bty = "n", cex = 0.9, x.intersp = 0.7)
}
subtitle <- function(text) mtext(text, side = 3, line = 0.4, cex = 0.78, col = muted)

# 1. Data acquisition, validation, and transparent preprocessing.
files <- c("baseline_no_signal_no_added_fibre_full_sweep.csv",
           "roadm_3dbm_injected_signal_no_added_fibre_full_sweep.csv",
           "roadm_7dbm_injected_signal_no_added_fibre_full_sweep.csv",
           "roadm_9dbm_injected_signal_no_added_fibre_full_sweep.csv",
           "roadm_12dbm_injected_signal_no_added_fibre_full_sweep.csv")
submission_csv <- "M5Assign 5.1 Data__Ronggon_Asif.csv"
candidates <- submission_csv
available <- candidates[file.exists(candidates)]
required_columns <- c("Noise_dBm", "Sweep_km", "VOA_dB", "QBER", "SecureKeyRate_bps")
read_source <- function(path) {
  d <- read.csv(path, check.names = FALSE)
  missing <- setdiff(required_columns, names(d))
  if (length(missing)) stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  if (!"SourceFile" %in% names(d)) d$SourceFile <- basename(path)
  d
}
if (length(available)) {
  source_paths <- available[1]
  dat <- read_source(source_paths)
  input_mode <- "Combined submission CSV"
} else {
  source_paths <- file.path("data/raw", files)
  for (i in seq_along(files)) {
    if (!file.exists(source_paths[i]))
      download.file(paste0("https://zenodo.org/records/21792844/files/", files[i], "?download=1"),
                    source_paths[i], mode = "wb")
  }
  raw_list <- lapply(source_paths, read_source)
  dat <- do.call(rbind, raw_list)
  input_mode <- "Five original Zenodo CSV files"
}
rownames(dat) <- NULL
if (!nrow(dat)) stop("No observations were found.")
num_cols <- c("Sweep_km", "VOA_dB", "QBER", "SecureKeyRate_bps")
if (!all(vapply(dat[num_cols], is.numeric, logical(1))))
  stop("Required measurement columns must be numeric.")
if (!all(vapply(dat[num_cols], function(x) all(is.finite(x)), logical(1))))
  stop("Required measurement columns contain missing or non-finite values.")
if (any(dat$Sweep_km != 0)) stop("This analysis requires only 0-km observations.")
if (any(dat$QBER < 0 | dat$QBER > 1) || any(dat$SecureKeyRate_bps < 0))
  stop("QBER must lie in [0,1] and secure key rate must be nonnegative.")
if (anyNA(dat$Noise_dBm) ||
    any(!as.character(dat$Noise_dBm) %in% c("none", "3", "7", "9", "12")))
  stop("Noise_dBm must contain only none, 3, 7, 9, or 12.")
cond_levels <- c("Baseline", "3 dBm", "7 dBm", "9 dBm", "12 dBm")
dat$NoiseCondition <- factor(ifelse(dat$Noise_dBm == "none", "Baseline",
                                  paste0(dat$Noise_dBm, " dBm")), levels = cond_levels)
if (any(table(dat$NoiseCondition) < 3))
  stop("Each condition needs at least three observations for stratified splitting.")
dat$LogSKR <- log10(dat$SecureKeyRate_bps + 1)
missing_total <- sum(is.na(dat))
duplicate_measurements <- sum(duplicated(dat[required_columns]))
# Equal measurement values are not proof of duplicate acquisition; retain them.
X <- dat[c("VOA_dB", "QBER", "LogSKR")]
predictor_sd <- vapply(X, sd, numeric(1))
if (any(!is.finite(predictor_sd) | predictor_sd <= 0))
  stop("All analysis predictors require finite, positive standard deviations.")
X_scaled <- scale(X)
cond_idx <- as.integer(dat$NoiseCondition)
write.csv(dat, "data/processed/QKD_Coexistence_0km_Combined.csv", row.names = FALSE)
write.csv(data.frame(Input = source_paths, MD5 = unname(tools::md5sum(source_paths))),
          "results/input_manifest.csv", row.names = FALSE)
writeLines(c(paste("Input mode:", input_mode),
             paste("Rows:", nrow(dat)), paste("Missing values across all fields:", missing_total),
             paste("Repeated five-field measurement tuples retained:", duplicate_measurements),
             "No imputation or outlier deletion was performed.",
             "Zero SKR retained; LogSKR = log10(SKR + 1); PCA/clusters use sample-SD standardization."),
           "results/preprocessing_log.txt")
desc_fun <- function(x) c(N = length(x), Mean = mean(x), SD = sd(x),
                          Median = median(x), Min = min(x), Max = max(x))
overall_desc <- data.frame(Variable = names(X), do.call(rbind, lapply(X, desc_fun)),
                           row.names = NULL)
write.csv(overall_desc, "results/descriptive_statistics_overall.csv", row.names = FALSE)
cond_desc <- do.call(rbind, lapply(cond_levels, function(g) {
  d <- dat[dat$NoiseCondition == g, ]
  data.frame(NoiseCondition = g, N = nrow(d), Mean_VOA_dB = mean(d$VOA_dB),
             SD_VOA_dB = sd(d$VOA_dB), Min_VOA_dB = min(d$VOA_dB), Max_VOA_dB = max(d$VOA_dB),
             Mean_QBER = mean(d$QBER), SD_QBER = sd(d$QBER),
             Median_SKR_bps = median(d$SecureKeyRate_bps),
             Mean_LogSKR = mean(d$LogSKR), SD_LogSKR = sd(d$LogSKR),
             Zero_SKR_pct = 100 * mean(d$SecureKeyRate_bps == 0))
}))
write.csv(cond_desc, "results/descriptive_statistics_by_condition.csv", row.names = FALSE)
cor_mat <- cor(X)
write.csv(cor_mat, "results/correlation_matrix.csv")

# 2. EDA: distribution geometry and observed attenuation profiles.
open_plot("01_scatterplot_matrix.pdf", 8, 7.4)
pairs(X, labels = unname(variable_labels), pch = 16, cex = 0.38,
      col = adjustcolor(condition_cols[cond_idx], 0.3), oma = c(5, 3, 4, 2),
      main = "Three measurements, one joint performance pattern",
      gap = 0.65, cex.labels = 1.05)
par(fig = c(0, 1, 0, 0.105), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend("center", legend = cond_levels, col = condition_cols, pch = 16,
       horiz = TRUE, bty = "n", cex = 0.9)
dev.off()

open_plot("02_correlation_matrix.pdf", 7.4, 5.7)
par(mar = c(5, 7.3, 4, 1.8))
pal_cor <- colorRampPalette(c("#A13F58", "#FAFAFC", "#176B80"))(201)
plot(NA, xlim = c(0.5, 3.5), ylim = c(0.5, 3.5), axes = FALSE,
     xlab = "", ylab = "", xaxs = "i", yaxs = "i", main = "Pearson correlation")
for (i in 1:3) for (j in 1:3) {
  value <- cor_mat[i, j]
  rect(j - 0.5, 3.5 - i, j + 0.5, 4.5 - i,
       col = pal_cor[round((value + 1) * 100) + 1], border = "white", lwd = 3)
  text(j, 4 - i, fmt(value, 2), cex = 1.35, font = 2,
       col = if (abs(value) > 0.65) "white" else ink)
}
axis(1, at = 1:3, labels = unname(variable_labels), tick = FALSE)
axis(2, at = 3:1, labels = unname(variable_labels), tick = FALSE, las = 1)
subtitle("Teal: positive association   |   Rose: negative association   |   Range: -1 to +1")
dev.off()

trend_plot <- function(response, filename, title_text, ylab_text) {
  open_plot(filename, legend = TRUE)
  plot(dat$VOA_dB, dat[[response]], type = "n", xlab = "Optical attenuation (dB)",
       ylab = ylab_text, main = title_text)
  plot_grid()
  for (i in seq_along(cond_levels)) {
    d <- dat[dat$NoiseCondition == cond_levels[i], ]
    points(d$VOA_dB, d[[response]], pch = 16, cex = 0.45,
           col = adjustcolor(condition_cols[i], 0.22))
    m <- aggregate(d[[response]], list(VOA = d$VOA_dB), mean)
    lines(m$VOA, m$x, col = condition_cols[i], lwd = 2.4)
  }
  subtitle("Faint points: measurements   |   Lines: means at observed attenuation settings")
  legend_strip(cond_levels, condition_cols)
  dev.off()
}
trend_plot("QBER", "03_qber_vs_attenuation.pdf", "Error rate rises across operating conditions", "QBER (proportion)")
trend_plot("LogSKR", "04_logskr_vs_attenuation.pdf", "Key throughput changes with attenuation", "log10(SKR + 1)")

# 3. Standardized PCA. Signs are fixed for presentation, not substantive inference.
pca <- prcomp(X, center = TRUE, scale. = TRUE)
for (j in seq_len(ncol(pca$rotation))) {
  anchor <- which.max(abs(pca$rotation[, j]))
  if (pca$rotation[anchor, j] < 0) {
    pca$rotation[, j] <- -pca$rotation[, j]
    pca$x[, j] <- -pca$x[, j]
  }
}
pca_var <- pca$sdev^2 / sum(pca$sdev^2)
pca_cum <- cumsum(pca_var)
n_pc_80 <- which(pca_cum >= 0.8)[1]
pca_variance <- data.frame(Component = paste0("PC", 1:3), Eigenvalue = pca$sdev^2,
                           VarianceExplained = pca_var, CumulativeVariance = pca_cum)
write.csv(pca_variance, "results/pca_variance.csv", row.names = FALSE)
write.csv(pca$rotation, "results/pca_loadings.csv")
write.csv(data.frame(pca$x, NoiseCondition = dat$NoiseCondition),
          "results/pca_scores.csv", row.names = FALSE)
scores <- as.data.frame(pca$x)
open_plot("05_pca_scree.pdf", 7.7, 5.1)
scree_x <- barplot(100 * pca_var, names.arg = paste0("PC", 1:3), col = c(accent, "#8CC8CF", "#DCEAF0"),
                   border = NA, ylim = c(0, 110), ylab = "Variance explained (%)",
                   main = "How many dimensions retain the information?")
text(scree_x, 100 * pca_var + 4, sprintf("%.1f%%", 100 * pca_var), font = 2)
lines(scree_x, 100 * pca_cum, type = "b", pch = 16, lwd = 1.7, col = "#A94374")
abline(h = 80, col = muted, lty = 3)
legend("right", legend = c("Individual component", "Cumulative variance", "80% criterion"),
       col = c(accent, "#A94374", muted), pch = c(15, 16, NA),
       lty = c(NA, 1, 3), bty = "n", cex = 0.82, bg = "white")
dev.off()

open_plot("06_pca_biplot.pdf", 8.3, 6.3, legend = TRUE)
load <- pca$rotation[, 1:2, drop = FALSE]
arrow_scale <- 0.65 * min(apply(abs(as.matrix(scores[1:2])), 2, max)) / max(abs(load))
arrow_xy <- load * arrow_scale
plot(scores$PC1, scores$PC2, type = "n",
     xlim = extendrange(c(scores$PC1, 1.3 * arrow_xy[, 1]), f = 0.09),
     ylim = extendrange(c(scores$PC2, 1.3 * arrow_xy[, 2]), f = 0.09),
     xlab = sprintf("PC1 (%.1f%%)", 100 * pca_var[1]),
     ylab = sprintf("PC2 (%.1f%%)", 100 * pca_var[2]),
     main = "PCA: common variation and secondary contrasts")
plot_grid()
abline(h = 0, v = 0, col = "#AAB5C3", lty = 3)
points(scores$PC1, scores$PC2, pch = condition_pch[cond_idx], cex = 0.53,
       col = adjustcolor(condition_cols[cond_idx], 0.38))
arrows(0, 0, arrow_xy[, 1], arrow_xy[, 2], length = 0.08, lwd = 2, col = ink)
text(arrow_xy[, 1] * 1.18, arrow_xy[, 2] * 1.18,
     labels = unname(variable_labels), cex = 0.9, font = 2, col = ink)
subtitle("Arrows are loadings multiplied by a common display factor; signs are conventional")
legend_strip(cond_levels, condition_cols, condition_pch)
dev.off()

# 4. K-means and Ward clustering, with evidence for and against the chosen cut.
k_range <- 1:min(10, nrow(unique(X_scaled)) - 1)
if (length(k_range) < 3) stop("Too few distinct measurements for clustering.")
km_candidates <- lapply(k_range, function(k) {
  set.seed(seeds["ElbowBase"] + k)
  kmeans(X_scaled, centers = k, nstart = 50, iter.max = 200)
})
wss <- vapply(km_candidates, function(x) x$tot.withinss, numeric(1))
# Maximum distance to the endpoint chord is a candidate-selection heuristic.
x <- k_range
dist_line <- abs((tail(wss, 1) - wss[1]) * x - (tail(x, 1) - x[1]) * wss +
                 tail(x, 1) * wss[1] - tail(wss, 1) * x[1]) /
             sqrt((tail(wss, 1) - wss[1])^2 + (tail(x, 1) - x[1])^2)
optimal_k <- max(2, k_range[which.max(dist_line)])
dmat <- dist(X_scaled)
silhouette_mean <- function(groups) mean(cluster::silhouette(groups, dmat)[, "sil_width"])
sil_km <- c(NA_real_, vapply(km_candidates[-1], function(m) silhouette_mean(m$cluster), numeric(1)))
hc <- hclust(dmat, method = "ward.D2")
sil_hc <- c(NA_real_, vapply(k_range[-1], function(k) silhouette_mean(cutree(hc, k)), numeric(1)))
validation <- data.frame(k = k_range, WSS = wss, DistanceFromChord = dist_line,
                         KmeansSilhouette = sil_km, WardSilhouette = sil_hc)
write.csv(validation, "results/kmeans_elbow.csv", row.names = FALSE)
write.csv(validation, "results/cluster_validation.csv", row.names = FALSE)
set.seed(seeds["FinalKmeans"])
km <- kmeans(X_scaled, centers = optimal_k, nstart = 100, iter.max = 200)
# Order cluster labels by mean log-SKR so profile labels have a stable meaning.
renumber <- function(groups) {
  order_ids <- as.integer(names(sort(tapply(dat$LogSKR, groups, mean))))
  match(groups, order_ids)
}
km_order <- as.integer(names(sort(tapply(dat$LogSKR, km$cluster, mean))))
km$cluster <- match(km$cluster, km_order)
km$centers <- km$centers[km_order, , drop = FALSE]
rownames(km$centers) <- as.character(seq_len(optimal_k))
km$withinss <- km$withinss[km_order]
km$size <- km$size[km_order]
hc_cluster <- renumber(cutree(hc, k = optimal_k))
km_sil <- silhouette_mean(km$cluster)
hc_sil <- silhouette_mean(hc_cluster)
choose2 <- function(x) x * (x - 1) / 2
cross_methods <- table(km$cluster, hc_cluster)
index <- sum(choose2(cross_methods))
expected <- sum(choose2(rowSums(cross_methods))) * sum(choose2(colSums(cross_methods))) / choose2(nrow(dat))
maximum <- (sum(choose2(rowSums(cross_methods))) + sum(choose2(colSums(cross_methods)))) / 2
ari <- if (maximum == expected) NA_real_ else (index - expected) / (maximum - expected)
write.csv(as.data.frame.matrix(cross_methods), "results/kmeans_vs_hierarchical.csv")
km_cross <- table(Cluster = km$cluster, NoiseCondition = dat$NoiseCondition)
hc_cross <- table(Cluster = hc_cluster, NoiseCondition = dat$NoiseCondition)
write.csv(as.data.frame.matrix(km_cross), "results/kmeans_vs_condition.csv")
write.csv(as.data.frame.matrix(hc_cross), "results/hierarchical_vs_condition.csv")
write.csv(data.frame(Observation = seq_len(nrow(dat)), Cluster = km$cluster,
                     WardCluster = hc_cluster, NoiseCondition = dat$NoiseCondition),
          "results/kmeans_membership.csv", row.names = FALSE)
cluster_profiles <- do.call(rbind, lapply(seq_len(optimal_k), function(k) {
  d <- dat[km$cluster == k, ]
  data.frame(Cluster = k, N = nrow(d), Mean_VOA_dB = mean(d$VOA_dB),
             Mean_QBER = mean(d$QBER), Mean_LogSKR = mean(d$LogSKR),
             Median_SKR_bps = median(d$SecureKeyRate_bps),
             Zero_SKR_pct = 100 * mean(d$SecureKeyRate_bps == 0))
}))
write.csv(cluster_profiles, "results/cluster_profiles.csv", row.names = FALSE)
cluster_cols <- hcl.colors(optimal_k, palette = "Dark 3")
open_plot("07_elbow_method.pdf", 9.1, 4.9)
par(mfrow = c(1, 2), mar = c(4.5, 4.7, 4, 1))
plot(k_range, wss, type = "n", xlab = "Number of clusters (k)", ylab = "Within-cluster sum of squares",
     main = "Elbow candidate")
plot_grid()
lines(k_range, wss, type = "b", pch = 16, col = accent, lwd = 2)
abline(v = optimal_k, lty = 2, col = "#A94374")
subtitle(paste("Chord-distance candidate: k =", optimal_k))
plot(k_range[-1], sil_km[-1], type = "n",
     ylim = range(c(sil_km, sil_hc), na.rm = TRUE), xlab = "Number of clusters (k)",
     ylab = "Mean silhouette width", main = "Check compactness and separation")
plot_grid()
lines(k_range[-1], sil_km[-1], type = "b", pch = 16, col = accent, lwd = 2)
lines(k_range[-1], sil_hc[-1], type = "b", pch = 17, col = "#A94374", lwd = 2)
abline(v = optimal_k, lty = 2, col = muted)
legend("bottomright", c("K-means", "Ward"), col = c(accent, "#A94374"),
       pch = c(16, 17), lty = 1, bty = "n", cex = 0.85)
dev.off()

cluster_projection <- function(groups, filename, title_text) {
  open_plot(filename, legend = TRUE)
  plot(scores$PC1, scores$PC2, type = "n", xlab = sprintf("PC1 (%.1f%%)", 100 * pca_var[1]),
       ylab = sprintf("PC2 (%.1f%%)", 100 * pca_var[2]), main = title_text)
  plot_grid()
  points(scores$PC1, scores$PC2, col = adjustcolor(cluster_cols[groups], 0.55), pch = 16, cex = 0.6)
  subtitle("Clustering uses all three standardized variables; the display uses two PCs")
  legend_strip(paste("Cluster", seq_len(optimal_k)), cluster_cols)
  dev.off()
}
cluster_projection(km$cluster, "08_kmeans_clusters_pca.pdf", sprintf("K-means: %d exploratory performance groups", optimal_k))
lower_cut <- hc$height[nrow(dat) - optimal_k]
upper_cut <- hc$height[nrow(dat) - optimal_k + 1]
cut_height <- mean(c(lower_cut, upper_cut))
open_plot("09_hierarchical_dendrogram.pdf", 9, 5.3)
plot(hc, labels = FALSE, hang = -1, main = "Ward hierarchy: where the tree is cut",
     xlab = sprintf("%d measurement rows; labels omitted for readability", nrow(dat)),
     ylab = "Ward merge height", sub = "", col = ink)
abline(h = cut_height, lty = 2, col = "#A94374", lwd = 1.6)
rect.hclust(hc, k = optimal_k, border = cluster_cols[unique(hc_cluster[hc$order])])
legend("topright", sprintf("k = %d; cut height %.1f", optimal_k, cut_height),
       lty = 2, col = "#A94374", bty = "n", cex = 0.9)
dev.off()
cluster_projection(hc_cluster, "10_hierarchical_clusters_pca.pdf", sprintf("Ward partition at the same k = %d", optimal_k))

# 5. LDA classification. Independent seed; training-only centering/scaling.
set.seed(seeds["LDASplit"])
train_idx <- sort(unlist(lapply(cond_levels, function(g) {
  ids <- which(dat$NoiseCondition == g)
  sample(ids, floor(0.70 * length(ids)))
})))
test_idx <- setdiff(seq_len(nrow(dat)), train_idx)
train <- dat[train_idx, c("NoiseCondition", names(X))]
test <- dat[test_idx, c("NoiseCondition", names(X))]
mu <- vapply(train[-1], mean, numeric(1))
sig <- vapply(train[-1], sd, numeric(1))
if (any(!is.finite(sig) | sig <= 0)) stop("Training predictors require positive finite SDs.")
standardize_frame <- function(d) {
  out <- data.frame(NoiseCondition = d$NoiseCondition)
  for (v in names(mu)) out[[paste0(v, "_z")]] <- (d[[v]] - mu[v]) / sig[v]
  out
}
train_z <- standardize_frame(train)
test_z <- standardize_frame(test)
all_z <- standardize_frame(dat)
lda_prior <- prop.table(table(train$NoiseCondition))
lda_fit <- MASS::lda(NoiseCondition ~ VOA_dB_z + QBER_z + LogSKR_z, data = train_z,
                     prior = as.numeric(lda_prior))
write.csv(data.frame(Condition = names(lda_prior), TrainingPrior = as.numeric(lda_prior)),
          "results/lda_class_priors.csv", row.names = FALSE)
lda_scaling <- lda_fit$scaling
lda_test_pred <- predict(lda_fit, newdata = test_z)
conf <- table(Actual = test_z$NoiseCondition,
              Predicted = factor(lda_test_pred$class, levels = cond_levels))
accuracy <- sum(diag(conf)) / sum(conf)
recall <- diag(conf) / rowSums(conf)
precision <- ifelse(colSums(conf) > 0, diag(conf) / colSums(conf), NA_real_)
macro_recall <- mean(recall)
majority_class <- names(which.max(table(train$NoiseCondition)))
majority_accuracy <- mean(test$NoiseCondition == majority_class)
class_metrics <- data.frame(Condition = cond_levels, Support = as.integer(rowSums(conf)),
                             PredictedN = as.integer(colSums(conf)),
                             Recall = as.numeric(recall), Precision = as.numeric(precision))
write.csv(class_metrics, "results/lda_class_metrics.csv", row.names = FALSE)
write.csv(data.frame(Metric = c("Held-out accuracy", "Training-majority baseline accuracy",
                               "Macro recall", "Training observations", "Test observations"),
                     Value = c(accuracy, majority_accuracy, macro_recall, nrow(train), nrow(test))),
          "results/lda_performance.csv", row.names = FALSE)
write.csv(as.data.frame.matrix(conf), "results/lda_confusion_matrix.csv")
write.csv(lda_scaling, "results/lda_coefficients.csv")
write.csv(lda_fit$means, "results/lda_group_means_standardized.csv")
write.csv(data.frame(Variable = names(mu), TrainingMean = mu, TrainingSD = sig),
          "results/lda_preprocessing.csv", row.names = FALSE)
split_provenance <- data.frame(Observation = seq_len(nrow(dat)), SourceFile = dat$SourceFile,
                               NoiseCondition = dat$NoiseCondition, VOA_dB = dat$VOA_dB,
                               Split = ifelse(seq_len(nrow(dat)) %in% train_idx, "Train", "Test"))
write.csv(split_provenance, "results/lda_split_provenance.csv", row.names = FALSE)
write.csv(data.frame(Observation = test_idx, Actual = test$NoiseCondition,
                     Predicted = lda_test_pred$class, lda_test_pred$posterior),
          "results/lda_test_predictions.csv", row.names = FALSE)
lda_all <- predict(lda_fit, newdata = all_z)
write.csv(data.frame(lda_all$x, NoiseCondition = dat$NoiseCondition, Split = split_provenance$Split),
          "results/lda_scores.csv", row.names = FALSE)
open_plot("11_lda_separation.pdf", legend = TRUE)
lda_y <- if (ncol(lda_all$x) >= 2) lda_all$x[, 2] else rep(0, nrow(dat))
plot(lda_all$x[, 1], lda_y, type = "n", xlab = "First discriminant (LD1)",
     ylab = if (ncol(lda_all$x) >= 2) "Second discriminant (LD2)" else "",
     main = "LDA scores: overlap remains among known conditions")
plot_grid()
points(lda_all$x[, 1], lda_y, pch = condition_pch[cond_idx], cex = 0.58,
       col = adjustcolor(condition_cols[cond_idx], 0.4))
subtitle("All observations projected using the training model; accuracy uses test rows only")
legend_strip(cond_levels, condition_cols, condition_pch)
dev.off()

# Counts and within-actual-class percentages make imbalance and failures visible.
open_plot("12_lda_confusion_matrix.pdf", 9, 7.5)
layout(matrix(c(1, 2), ncol = 1), heights = c(1, 0.18))
par(mar = c(4.7, 7.1, 4.6, 1.6), xaxs = "i", yaxs = "i")
row_prop <- sweep(as.matrix(conf), 1, rowSums(conf), "/")
heat_cols <- colorRampPalette(c("#F2F7FB", "#C3DDEA", "#589EBB", "#144766"))(101)
plot(NA, xlim = c(0.5, 5.5), ylim = c(0.5, 5.5), axes = FALSE,
     xlab = "Predicted condition", ylab = "", main = "Where the classifier succeeds and fails")
subtitle(sprintf("Held-out accuracy %.1f%%  |  Majority baseline %.1f%%  |  Macro recall %.1f%%",
                 100 * accuracy, 100 * majority_accuracy, 100 * macro_recall))
for (i in 1:5) for (j in 1:5) {
  y <- 6 - i
  rect(j - 0.5, y - 0.5, j + 0.5, y + 0.5, col = heat_cols[round(100 * row_prop[i, j]) + 1],
       border = "white", lwd = 2)
  txt_col <- if (row_prop[i, j] >= 0.68) "white" else ink
  text(j, y + 0.11, conf[i, j], cex = 1.18, font = 2, col = txt_col)
  text(j, y - 0.17, sprintf("%.1f%%", 100 * row_prop[i, j]), cex = 0.87, col = txt_col)
  if (i == j) rect(j - 0.46, y - 0.46, j + 0.46, y + 0.46,
                  border = "#C38B30", lwd = 2)
}
axis(1, at = 1:5, labels = cond_levels, tick = FALSE, cex.axis = 1)
axis(2, at = 5:1, labels = paste0(cond_levels, "\n(n = ", rowSums(conf), ")"),
     tick = FALSE, las = 1, cex.axis = 0.94)
mtext("Actual condition", side = 2, line = 5.8, las = 0, cex = 0.95)
par(mar = c(2.7, 7.1, 0.6, 1.6))
plot(NA, xlim = c(0, 100), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "",
     xaxs = "i", yaxs = "i")
for (j in 0:99) rect(j, 0.35, j + 1, 0.7, col = heat_cols[j + 1], border = NA)
axis(1, at = seq(0, 100, 20), labels = paste0(seq(0, 100, 20), "%"), cex.axis = 0.82)
mtext("Within-actual-class share; each row sums to 100%. Gold outlines mark correct predictions.",
      side = 1, line = 1.6, cex = 0.78, col = muted)
dev.off()

# Cluster profile plot, with both original-unit table and standardized geometry.
profile_z <- aggregate(as.data.frame(X_scaled), list(Cluster = km$cluster), mean)
open_plot("13_cluster_profiles.pdf", 8.4, 5.3, legend = TRUE)
plot(1:3, rep(0, 3), type = "n", xlim = c(0.85, 3.15),
     ylim = extendrange(as.matrix(profile_z[-1])), xaxt = "n",
     xlab = "", ylab = "Cluster mean (standardized units)",
     main = "What distinguishes the selected K-means groups?")
plot_grid()
abline(h = 0, lty = 2, col = muted)
axis(1, at = 1:3, labels = unname(variable_labels))
for (k in seq_len(optimal_k))
  lines(1:3, as.numeric(profile_z[k, -1]), type = "b", pch = 16, lwd = 2.3,
        cex = 1.2, col = cluster_cols[k])
subtitle("Zero is the full-sample mean; labels are ordered from lower to higher mean log-SKR")
legend_strip(paste0("Cluster ", seq_len(optimal_k), " (n=", cluster_profiles$N, ")"), cluster_cols)
dev.off()

# 6. Multivariate multiple regression: two continuous outcomes jointly modeled.
# Centering makes condition intercept differences refer to mean attenuation.
# These are descriptive, within-campaign associations, not causal estimates.
voa_center <- mean(dat$VOA_dB)
dat$VOA_c <- dat$VOA_dB - voa_center
reg_formula <- cbind(QBER, LogSKR) ~ VOA_c * NoiseCondition
reg_fit <- lm(reg_formula, data = dat)
if (reg_fit$rank < nrow(coef(reg_fit))) stop("The regression design is rank deficient.")
manova_fit <- manova(reg_formula, data = dat)
manova_stats <- summary(manova_fit, test = "Pillai")$stats
manova_stats <- manova_stats[rownames(manova_stats) != "Residuals", , drop = FALSE]
reg_tests <- data.frame(Term = rownames(manova_stats), manova_stats, row.names = NULL, check.names = FALSE)
write.csv(reg_tests, "results/regression_manova.csv", row.names = FALSE)
outcomes <- c("QBER", "LogSKR")
outcome_fits <- lapply(outcomes, function(y) lm(as.formula(paste(y, "~ VOA_c * NoiseCondition")), data = dat))
names(outcome_fits) <- outcomes
reg_summary <- do.call(rbind, lapply(outcomes, function(y) {
  s <- summary(outcome_fits[[y]])
  data.frame(Outcome = y, R_squared = s$r.squared, Adjusted_R_squared = s$adj.r.squared,
             RMSE = sqrt(mean(residuals(outcome_fits[[y]])^2)), Residual_df = df.residual(outcome_fits[[y]]))
}))
write.csv(reg_summary, "results/regression_fit.csv", row.names = FALSE)
reg_coef <- coef(reg_fit)
reg_se <- sapply(outcome_fits, function(m) summary(m)$coefficients[, "Std. Error"])
term_labels <- rownames(reg_coef)
term_labels[term_labels == "(Intercept)"] <- "Baseline at mean attenuation"
term_labels[term_labels == "VOA_c"] <- "Baseline slope per dB"
term_labels <- sub("^VOA_c:NoiseCondition", "Extra slope: ", term_labels)
term_labels <- sub("^NoiseCondition", "At center: ", term_labels)
coef_table <- data.frame(Term = term_labels, QBER_Estimate = reg_coef[, "QBER"],
                         QBER_SE = reg_se[, "QBER"], LogSKR_Estimate = reg_coef[, "LogSKR"],
                         LogSKR_SE = reg_se[, "LogSKR"], row.names = NULL)
write.csv(coef_table, "results/regression_coefficients.csv", row.names = FALSE)
reg_slopes <- do.call(rbind, lapply(cond_levels, function(g) {
  term <- paste0("VOA_c:NoiseCondition", g)
  effect <- reg_coef["VOA_c", ]
  if (g != "Baseline") effect <- effect + reg_coef[term, ]
  data.frame(Condition = g, QBER_per_dB = effect["QBER"], LogSKR_per_dB = effect["LogSKR"])
}))
write.csv(reg_slopes, "results/regression_condition_slopes.csv", row.names = FALSE)
open_plot("14_regression_fits.pdf", 10, 5.4)
par(mfrow = c(1, 2), mar = c(4.7, 4.5, 4.2, 1), oma = c(4.5, 0, 0, 0))
for (response in outcomes) {
  plot(dat$VOA_dB, dat[[response]], type = "n", xlab = "Attenuation (dB)",
       ylab = if (response == "QBER") "QBER" else "log10(SKR + 1)",
       main = if (response == "QBER") "Regression: error rate" else "Regression: transformed throughput")
  plot_grid()
  for (i in seq_along(cond_levels)) {
    d <- dat[dat$NoiseCondition == cond_levels[i], ]
    points(d$VOA_dB, d[[response]], pch = 16, cex = 0.4, col = adjustcolor(condition_cols[i], 0.2))
    new <- data.frame(VOA_dB = seq(min(d$VOA_dB), max(d$VOA_dB), length.out = 100),
                      NoiseCondition = factor(cond_levels[i], levels = cond_levels))
    new$VOA_c <- new$VOA_dB - voa_center
    lines(new$VOA_dB, predict(outcome_fits[[response]], newdata = new), col = condition_cols[i], lwd = 2.4)
  }
  subtitle(sprintf("Within-sample R-squared = %.3f; lines limited to each observed range",
                   reg_summary$R_squared[reg_summary$Outcome == response]))
}
par(fig = c(0, 1, 0, 0.10), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend("center", cond_levels, col = condition_cols, lty = 1, lwd = 2, bty = "n", horiz = TRUE, cex = 0.9)
dev.off()

open_plot("15_regression_diagnostics.pdf", 9.6, 7.1)
par(mfrow = c(2, 2), mar = c(4, 4.6, 3, 1), oma = c(1.7, 0, 0, 0))
for (response in outcomes) {
  m <- outcome_fits[[response]]
  fit <- fitted(m)
  res <- residuals(m)
  plot(fit, res, pch = 16, cex = 0.42, col = adjustcolor(accent, 0.28),
       xlab = paste("Fitted", response), ylab = "Residual",
       main = paste(response, ": residuals versus fitted"))
  abline(h = 0, col = muted, lty = 2)
  lines(lowess(fit, res), col = "#A94374", lwd = 2)
  qqnorm(res / summary(m)$sigma, pch = 16, cex = 0.4,
         col = adjustcolor(accent, 0.35), main = paste(response, ": normal Q-Q"),
         xlab = "Theoretical normal quantile", ylab = "Residual / residual SD")
  qqline(res / summary(m)$sigma, col = "#A94374", lwd = 2)
}
mtext("Curvature, changing spread, or tail departures limit the linear-Gaussian model; repeated rows may also be dependent.",
      outer = TRUE, side = 1, line = 0.3, cex = 0.78, col = muted)
dev.off()
outside_qber <- sum(fitted(reg_fit)[, "QBER"] < 0 | fitted(reg_fit)[, "QBER"] > 1)
outside_logskr <- sum(fitted(reg_fit)[, "LogSKR"] < 0)
reg_diagnostics <- do.call(rbind, lapply(outcomes, function(y) {
  m <- outcome_fits[[y]]
  residual <- residuals(m)
  centered <- residual - mean(residual)
  data.frame(Outcome = y,
             ResidualSkewness = mean(centered^3) / mean(centered^2)^1.5,
             AbsResidualFittedCorrelation = cor(abs(residual), fitted(m)),
             MaximumCooksDistance = max(cooks.distance(m)),
             FittedOutsidePhysicalRange = if (y == "QBER") outside_qber else outside_logskr)
}))
write.csv(reg_diagnostics, "results/regression_diagnostics.csv", row.names = FALSE)
writeLines(c("Model: cbind(QBER, LogSKR) ~ VOA_c * NoiseCondition",
             paste("Attenuation center (dB):", format(voa_center, digits = 12)),
             "Pillai tests are sequential (Type I) in the displayed formula order.",
             "Classical standard errors and reference p-values assume independent errors and suitable covariance structure.",
             "Repeated within-setting measurements may be dependent; no run-level identifiers support independent-run validation.",
             "QBER is bounded; log-SKR has a point mass at zero. Linear trends and Gaussian residual assumptions are approximations.",
             paste("Observed-row fitted QBER outside [0,1]:", outside_qber),
             paste("Observed-row fitted log-SKR below zero:", outside_logskr),
             "Residual/fitted and Q-Q plots are diagnostic evidence; no causal conclusions or population-level significance claims are made."),
           "results/regression_diagnostics.txt")

# 7. Machine-readable summaries and complete LaTeX tables.
pc1_order <- names(sort(abs(pca$rotation[, 1]), decreasing = TRUE))
pc2_order <- names(sort(abs(pca$rotation[, 2]), decreasing = TRUE))
ld1_order <- rownames(lda_scaling)[order(abs(lda_scaling[, 1]), decreasing = TRUE)]
ld2_order <- if (ncol(lda_scaling) >= 2) rownames(lda_scaling)[order(abs(lda_scaling[, 2]), decreasing = TRUE)] else "N/A"
macro_values <- c(NObs = as.character(nrow(dat)), NMissing = as.character(missing_total),
                  PCOneVar = pct(pca_var[1]), PCTwoVar = pct(pca_var[2]), PCTwoCum = pct(pca_cum[2]),
                  NPCOptimal = as.character(n_pc_80), PCOneTop = escape_tex(pc1_order[1]),
                  PCTwoTop = escape_tex(pc2_order[1]), KOptimal = as.character(optimal_k),
                  LDAAccuracy = pct(accuracy), TrainN = as.character(nrow(train)),
                  TestN = as.character(nrow(test)), LDOneTop = escape_tex(sub("_z$", "", ld1_order[1])),
                  LDTwoTop = escape_tex(sub("_z$", "", ld2_order[1])),
                  LDAMajorityAccuracy = pct(majority_accuracy), LDAMacroRecall = pct(macro_recall),
                  ClusterARI = fmt(ari), KmeansSilhouette = fmt(km_sil), WardSilhouette = fmt(hc_sil),
                  RegressionQberRsq = fmt(reg_summary$R_squared[1]),
                  RegressionLogSkrRsq = fmt(reg_summary$R_squared[2]),
                  AttenuationCenter = fmt(voa_center, 2))
writeLines(sprintf("\\newcommand{\\%s}{%s}", names(macro_values), macro_values), "results/results.tex")
tex_table("descriptive_table.tex",
          c("Condition", "$n$", "Mean QBER", "Median SKR (bps)", "Zero-SKR (\\%)"),
          cbind(cond_desc$NoiseCondition, cond_desc$N, fmt(cond_desc$Mean_QBER, 4),
                format(round(cond_desc$Median_SKR_bps), big.mark = ",", scientific = FALSE),
                fmt(cond_desc$Zero_SKR_pct, 1)))
tex_table("pca_table.tex", c("Component", "Eigenvalue", "Variance (\\%)", "Cumulative (\\%)"),
          cbind(pca_variance$Component, fmt(pca_variance$Eigenvalue),
                fmt(100 * pca_var, 1), fmt(100 * pca_cum, 1)))
tex_table("pca_loadings_table.tex", c("Variable", "PC1", "PC2", "PC3"),
          cbind(escape_tex(rownames(pca$rotation)), apply(pca$rotation, 2, fmt)))
tex_table("lda_coefficients_table.tex", c("Predictor", colnames(lda_scaling)),
          cbind(escape_tex(sub("_z$", "", rownames(lda_scaling))), apply(lda_scaling, 2, fmt)))
tex_table("lda_confusion_table.tex", c("Actual $\\backslash$ Pred.", colnames(conf)),
          cbind(rownames(conf), as.matrix(conf)))
tex_table("cluster_profiles_table.tex",
          c("Cluster", "$n$", "Mean VOA", "Mean QBER", "Mean log-SKR", "Zero SKR (\\%)"),
          cbind(cluster_profiles$Cluster, cluster_profiles$N, fmt(cluster_profiles$Mean_VOA_dB, 2),
                fmt(cluster_profiles$Mean_QBER, 4), fmt(cluster_profiles$Mean_LogSKR, 2),
                fmt(cluster_profiles$Zero_SKR_pct, 1)))
tex_table("lda_metrics_table.tex", c("Condition", "Test $n$", "Predicted $n$", "Recall (\\%)", "Precision (\\%)"),
          cbind(class_metrics$Condition, class_metrics$Support, class_metrics$PredictedN,
                fmt(100 * class_metrics$Recall, 1),
                ifelse(is.na(class_metrics$Precision), "--", fmt(100 * class_metrics$Precision, 1))))
test_labels <- c("Attenuation (centered)", "Noise condition", "Attenuation $\\times$ condition")
tex_table("regression_tests_table.tex",
          c("Sequential term", "Pillai", "Approx. $F$", "df", "Reference $p$"),
          cbind(test_labels, fmt(manova_stats[, "Pillai"]), fmt(manova_stats[, "approx F"], 2),
                paste0(manova_stats[, "num Df"], ", ", manova_stats[, "den Df"]),
                tex_p(manova_stats[, "Pr(>F)"])))
tex_table("regression_coefficients_table.tex",
          c("Term", "QBER est.", "QBER SE", "Log-SKR est.", "Log-SKR SE"),
          cbind(coef_table$Term, fmt(coef_table$QBER_Estimate, 4), fmt(coef_table$QBER_SE, 4),
                fmt(coef_table$LogSKR_Estimate, 3), fmt(coef_table$LogSKR_SE, 3)))
tex_table("regression_fit_table.tex", c("Outcome", "$R^2$", "Adjusted $R^2$", "RMSE", "Residual df"),
          cbind(c("QBER", "Log-SKR"), fmt(reg_summary$R_squared), fmt(reg_summary$Adjusted_R_squared),
                fmt(reg_summary$RMSE, 4), reg_summary$Residual_df))

# 8. Computed interpretation. Narratives are generated from the same objects as plots.
pair_idx <- which(lower.tri(cor_mat), arr.ind = TRUE)
strong_idx <- pair_idx[which.max(abs(cor_mat[pair_idx])), ]
write_findings("eda", c(
  sprintf("The analysis contains %d measurements across five conditions, with group sizes ranging from %d to %d. Attenuation spans %.1f--%.1f dB, QBER spans %.4f--%.4f, and raw SKR spans %.0f--%.0f bits/s. There are %d zero-SKR measurements (%s); these are retained, and the $\\log_{10}(\\mathrm{SKR}+1)$ transformation preserves their value at zero.",
          nrow(dat), min(cond_desc$N), max(cond_desc$N), min(dat$VOA_dB), max(dat$VOA_dB),
          min(dat$QBER), max(dat$QBER), min(dat$SecureKeyRate_bps), max(dat$SecureKeyRate_bps),
          sum(dat$SecureKeyRate_bps == 0), pct(mean(dat$SecureKeyRate_bps == 0))),
  sprintf("The strongest absolute pairwise correlation is between %s and %s ($r=%.3f$). In particular, attenuation correlates %.3f with QBER and %.3f with transformed SKR. These pooled associations combine within-condition trends and differences in the sampled condition/attenuation mix; they do not by themselves isolate a causal effect of injected power.",
          escape_tex(rownames(cor_mat)[strong_idx[1]]), escape_tex(colnames(cor_mat)[strong_idx[2]]),
          cor_mat[strong_idx[1], strong_idx[2]], cor_mat["VOA_dB", "QBER"], cor_mat["VOA_dB", "LogSKR"]),
  sprintf("For example, mean sampled attenuation is %.2f dB for baseline and %.2f dB for 12 dBm, while their median raw SKR values are %.0f and %.0f bits/s, respectively. Those medians summarize different attenuation distributions. A raw median comparison must therefore not be interpreted as a controlled injected-power effect.",
          cond_desc$Mean_VOA_dB[1], cond_desc$Mean_VOA_dB[5],
          cond_desc$Median_SKR_bps[1], cond_desc$Median_SKR_bps[5])))
loading_sentence <- function(j) paste(sprintf("%s $%+.3f$",
                                               escape_tex(rownames(pca$rotation)), pca$rotation[, j]),
                                      collapse = ", ")
component_contrast <- function(j) {
  loadings <- pca$rotation[, j]
  positive <- paste(unname(variable_labels[names(loadings)[loadings > 0]]), collapse = ", ")
  negative <- paste(unname(variable_labels[names(loadings)[loadings < 0]]), collapse = ", ")
  paste0("Higher PC", j, " scores combine higher ", positive,
         if (nzchar(negative)) paste0(" with lower ", negative) else "", ".")
}
write_findings("pca", c(
  sprintf("PC1 retains %s of standardized variance, PC2 retains %s, and their cumulative share is %s. The 80\\%% retention rule selects %d component(s); the two-dimensional display is still useful for inspecting the smaller second-direction contrast.",
          pct(pca_var[1]), pct(pca_var[2]), pct(pca_cum[2]), n_pc_80),
  sprintf("The displayed PC1 loadings are %s; PC2 loadings are %s. Variables with the same loading sign move together along a component, whereas opposite signs define a contrast. The largest absolute contributions are %s on PC1 and %s on PC2. Component signs are arbitrary and have been oriented so each largest loading is positive; reversing an entire component changes no statistical conclusion.",
          loading_sentence(1), loading_sentence(2), escape_tex(pc1_order[1]), escape_tex(pc2_order[1])),
  paste(component_contrast(1), component_contrast(2))))
best_sil_k <- k_range[which.max(replace(sil_km, is.na(sil_km), -Inf))]
best_hc_k <- k_range[which.max(replace(sil_hc, is.na(sil_hc), -Inf))]
profile_sentences <- vapply(seq_len(optimal_k), function(k) {
  r <- cluster_profiles[k, ]
  composition <- km_cross[k, ] / sum(km_cross[k, ])
  major <- which.max(composition)
  sprintf("Cluster %d ($n=%d$) has mean attenuation %.2f dB, mean QBER %.4f, mean log-SKR %.2f, and %.1f\\%% zero SKR; its largest condition share is %s (%s).",
          k, r$N, r$Mean_VOA_dB, r$Mean_QBER, r$Mean_LogSKR, r$Zero_SKR_pct,
          cond_levels[major], pct(composition[major]))
}, character(1))
write_findings("clusters", c(
  sprintf("The geometric elbow selects the candidate $k=%d$. At this cut, mean silhouette width is %.3f for the final K-means solution and %.3f for Ward clustering. Across the tested cuts, the largest mean silhouettes occur at $k=%d$ (K-means) and $k=%d$ (Ward). The elbow is a compact-summary heuristic, not proof of a uniquely optimal number of physical regimes.",
          optimal_k, km_sil, hc_sil, best_sil_k, best_hc_k),
  sprintf("The Ward dendrogram cut lies between merge heights %.2f and %.2f. The adjusted Rand index comparing K-means and Ward memberships is %.3f (1 denotes identical partitions; approximately 0 denotes chance-level agreement). Both methods use the same standardized variables, so agreement is descriptive and does not establish independent validation.",
          lower_cut, upper_cut, ari),
  paste(profile_sentences, collapse = " ")))
worst <- which.min(recall)
never <- cond_levels[colSums(conf) == 0]
never_text <- if (length(never)) paste0(" No test observations were predicted as ",
                                       paste(never, collapse = " or "), "; precision is therefore undefined for those classes.") else ""
write_findings("lda", c(
  sprintf("The fixed stratified split uses %d training and %d test rows. LDA correctly classifies %d of %d test observations, giving accuracy %s. Always predicting the training-majority class (%s) gives %s on the same test rows, so the accuracy difference is %.1f percentage points. Macro recall is %s, which weights the five conditions equally.",
          nrow(train), nrow(test), sum(diag(conf)), sum(conf), pct(accuracy), majority_class,
          pct(majority_accuracy), 100 * (accuracy - majority_accuracy), pct(macro_recall)),
  paste0(sprintf("Class recall ranges from %s to %s; the lowest recall includes %s (%d correctly classified out of %d). Each heatmap cell shows a count and its percentage of the actual class, and the gold diagonal marks correct predictions. Off-diagonal mass therefore makes the model's errors visible even when class sizes differ.",
                pct(min(recall)), pct(max(recall)), cond_levels[worst], conf[worst, worst], rowSums(conf)[worst]),
         never_text),
  "These estimates concern held-out rows from the same campaign. Repeated measurements at shared settings can occur in both partitions; this is not validation on independent runs or unseen operating conditions. LDA's common within-class covariance and approximately Gaussian predictor assumptions are also approximations for bounded QBER and zero-heavy SKR."))
slope_base <- reg_slopes[reg_slopes$Condition == "Baseline", ]
interaction_row <- which(rownames(manova_stats) == "VOA_c:NoiseCondition")
write_findings("regression", c(
  sprintf("Multivariate multiple regression jointly models QBER and log-SKR using attenuation centered at %.2f dB, five noise-condition levels, and their interaction. Baseline is the reference condition. At the center, its fitted QBER is %.4f and fitted log-SKR is %.3f; its slopes per additional dB are %+.4f and %+.3f, respectively. Condition coefficients describe differences at the center, and interaction coefficients describe differences from the baseline slopes.",
          voa_center, reg_coef["(Intercept)", "QBER"], reg_coef["(Intercept)", "LogSKR"],
          slope_base$QBER_per_dB, slope_base$LogSKR_per_dB),
  sprintf("The within-sample $R^2$ values are %.3f for QBER and %.3f for log-SKR. For the attenuation-by-condition interaction, the sequential Pillai statistic is %.3f, with approximate $F=%.2f$ and reference $p$ %s. The displayed MANOVA tests are sequential in formula order; with an interaction present, main-effect coefficients must be interpreted at the reference condition or centered attenuation.",
          reg_summary$R_squared[1], reg_summary$R_squared[2],
          manova_stats[interaction_row, "Pillai"], manova_stats[interaction_row, "approx F"],
          tex_p(manova_stats[interaction_row, "Pr(>F)"])),
  sprintf("The residual-versus-fitted and normal Q--Q plots assess curvature, unequal spread, and non-normal tails. Among observed rows, %d fitted QBER values fall outside $[0,1]$ and %d fitted log-SKR values fall below zero, indicating whether the unconstrained linear model crosses physical boundaries. Standard errors and reference tests assume an error structure that repeated within-setting measurements may violate. The regression is therefore a descriptive association model for this campaign, not evidence of causation or confirmed population-level significance.",
          outside_qber, outside_logskr)))
# Explicit, computed explanations for the assignment's interpretation criteria.
append_finding <- function(name, paragraph)
  cat("\n\n", paragraph, "\n", file = file.path("results", paste0("findings_", name, ".tex")),
      append = TRUE, sep = "")
if (all(c("Baseline", "12 dBm") %in% as.character(dat$NoiseCondition[dat$VOA_dB == 28]))) {
  matched_mean <- vapply(c("Baseline", "12 dBm"), function(g)
    mean(dat$SecureKeyRate_bps[dat$NoiseCondition == g & dat$VOA_dB == 28]), numeric(1))
  append_finding("eda", sprintf(
    "At the shared attenuation of 28 dB, mean SKR is %.0f bits/s for baseline and %.0f bits/s for 12 dBm. The 12 dBm condition also has %.1f\\%% zero-throughput readings overall. This matched-setting comparison explains why its high pooled median is not evidence of better performance under stronger injection.",
    matched_mean[1], matched_mean[2], cond_desc$Zero_SKR_pct[5]))
}
if (pca$rotation["VOA_dB", 1] * pca$rotation["QBER", 1] > 0 &&
    pca$rotation["QBER", 1] * pca$rotation["LogSKR", 1] < 0)
  append_finding("pca", "PC1 can therefore be interpreted as an attenuation/error--throughput deterioration axis: it contrasts greater attenuation and error with lower key throughput. PC2 mainly captures the attenuation contrast left after that shared direction is removed. The near-one-dimensional pattern explains why the variables can be compressed, while the smaller contrasts remain relevant for distinguishing conditions.")
if (all(km_cross > 0))
  append_finding("clusters", "Every fitted cluster contains observations from all five experimental conditions. The clusters therefore describe near-collapse, reduced-throughput, and higher-throughput performance regimes rather than recovering the injection labels. Labels are ordered by increasing mean log-SKR, so their numeric ordering follows throughput in this solution.")
append_finding("lda", sprintf(
  "The largest absolute standardized coefficient in LD1 is %s (%+.3f), and in LD2 it is %s (%+.3f). These identify the strongest conditional contributions to the two plotted directions, not independent causal effects. Because attenuation, error, and throughput are strongly correlated, substantial coefficients can coexist with poor class separation.",
  escape_tex(sub("_z$", "", ld1_order[1])), lda_scaling[ld1_order[1], 1],
  escape_tex(sub("_z$", "", ld2_order[1])), lda_scaling[ld2_order[1], 2]))
error_details <- vapply(seq_along(cond_levels), function(i) {
  row <- as.numeric(conf[i, ])
  row[i] <- 0
  j <- which.max(row)
  if (sum(row) == 0)
    return(sprintf("%s has recall %s with no errors.", cond_levels[i], pct(recall[i])))
  sprintf("%s has recall %s; its most frequent error is prediction as %s (%d of %d actual observations, %s).",
          cond_levels[i], pct(recall[i]), cond_levels[j], conf[i, j], rowSums(conf)[i],
          pct(conf[i, j] / rowSums(conf)[i]))
}, character(1))
append_finding("lda", paste(error_details, collapse = " "))

summary_lines <- c(sprintf("Observations: %d", nrow(dat)),
                   sprintf("Missing values: %d", missing_total),
                   sprintf("PC1 variance explained: %.2f%%", 100 * pca_var[1]),
                   sprintf("PC2 variance explained: %.2f%%", 100 * pca_var[2]),
                   sprintf("PC1+PC2 cumulative variance: %.2f%%", 100 * pca_cum[2]),
                   sprintf("PCs needed for >=80%% variance: %d", n_pc_80),
                   sprintf("Candidate k by elbow geometry: %d", optimal_k),
                   sprintf("K-means silhouette: %.3f; Ward silhouette: %.3f; ARI: %.3f", km_sil, hc_sil, ari),
                   sprintf("LDA held-out accuracy: %.2f%%; majority baseline: %.2f%%; macro recall: %.2f%%",
                           100 * accuracy, 100 * majority_accuracy, 100 * macro_recall),
                   sprintf("Regression R-squared: QBER %.3f; LogSKR %.3f",
                           reg_summary$R_squared[1], reg_summary$R_squared[2]))
writeLines(summary_lines, "results/analysis_summary.txt")
capture.output(sessionInfo(), file = "results/sessionInfo.txt")
saveRDS(list(pca = pca, kmeans = km, hierarchical = hc, lda = lda_fit,
             regression = reg_fit, train_rows = train_idx, seeds = seeds),
        "results/fitted_models.rds")
message("Analysis complete. Generated CSV summaries, 15 figures, LaTeX tables and findings.")
message(paste(summary_lines, collapse = "\n"))
