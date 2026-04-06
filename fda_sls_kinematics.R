# Functional Data Analysis: Single Leg Squat Kinematics PRE vs POST
# B-spline smoothing with pointwise 95% CI for paired comparisons
# Analyzes ankle, knee, hip angles in sagittal and frontal planes

library(readxl)
library(fda)

# --- User Settings ---

input_dir  <- "data/sls_fda"
output_dir <- "output"

pre_file  <- file.path(input_dir, "pre_fda_sls.xlsx")
post_file <- file.path(input_dir, "post_fda_sls.xlsx")

# Subjects to exclude (coordinate system errors or missing data)
exclude_subs <- c("sub_09", "sub_13", "sub_16")

sheets <- c("ankle_x", "knee_x", "hip_x", "ankle_y", "knee_y", "hip_y")

angle_labels <- c(
  ankle_x = "Ankle - Sagittal Plane",
  knee_x  = "Knee - Sagittal Plane",
  hip_x   = "Hip - Sagittal Plane",
  ankle_y = "Ankle - Frontal Plane",
  knee_y  = "Knee - Frontal Plane",
  hip_y   = "Hip - Frontal Plane"
)

# FDA parameters
n_basis      <- 19
spline_order <- 4
lambda_val   <- 1e-4
alpha        <- 0.05

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- Helper: Read and Clean One Sheet ---

read_angle_sheet <- function(filepath, sheet_name, exclude) {
  dat <- read_excel(filepath, sheet = sheet_name)
  
  time_col <- dat[[1]]
  dat <- dat[, -1]
  
  sub_ids <- sub("^(sub_\\d+)_.*", "\\1", colnames(dat))
  colnames(dat) <- sub_ids
  
  keep <- !(sub_ids %in% exclude)
  dat <- dat[, keep, drop = FALSE]
  
  mat <- as.matrix(dat)
  storage.mode(mat) <- "double"
  
  if (any(is.na(mat))) {
    na_subs <- colnames(mat)[colSums(is.na(mat)) > 0]
    warning(paste0("NAs found in ", sheet_name, " for: ", paste(na_subs, collapse = ", ")))
  }
  
  list(matrix = mat, time = as.numeric(time_col), subjects = colnames(dat))
}

# --- Read All Data ---

cat("Reading data...\n")

pre_data  <- lapply(sheets, function(s) read_angle_sheet(pre_file, s, exclude_subs))
post_data <- lapply(sheets, function(s) read_angle_sheet(post_file, s, exclude_subs))
names(pre_data)  <- sheets
names(post_data) <- sheets

for (s in sheets) {
  if (!identical(pre_data[[s]]$subjects, post_data[[s]]$subjects)) {
    stop(paste("Subject mismatch in sheet:", s))
  }
}

n_sub    <- ncol(pre_data[[1]]$matrix)
n_time   <- nrow(pre_data[[1]]$matrix)
time_pts <- pre_data[[1]]$time

cat(sprintf("Subjects retained: %d\n", n_sub))
cat(sprintf("Time points: %d (0-100%% SLS cycle)\n\n", n_time))

# --- Create B-spline Basis ---

t_range <- range(time_pts)

basis_obj <- create.bspline.basis(rangeval = t_range, nbasis = n_basis, norder = spline_order)
fdPar_obj <- fdPar(fdobj = basis_obj, Lfdobj = 2, lambda = lambda_val)

cat(sprintf("B-spline basis: %d functions, order %d, lambda = %g\n\n", n_basis, spline_order, lambda_val))

# --- Smooth All Curves ---

cat("Smoothing curves...\n")

smooth_pre  <- list()
smooth_post <- list()

for (s in sheets) {
  smooth_pre[[s]]  <- smooth.basis(argvals = time_pts, y = pre_data[[s]]$matrix, fdParobj = fdPar_obj)
  smooth_post[[s]] <- smooth.basis(argvals = time_pts, y = post_data[[s]]$matrix, fdParobj = fdPar_obj)
}

# --- Evaluate and Compute Statistics ---

eval_pts <- seq(t_range[1], t_range[2], length.out = 101)
results <- list()

for (s in sheets) {
  pre_vals  <- eval.fd(eval_pts, smooth_pre[[s]]$fd)
  post_vals <- eval.fd(eval_pts, smooth_post[[s]]$fd)
  diff_vals <- post_vals - pre_vals
  
  diff_mean <- rowMeans(diff_vals)
  diff_sd   <- apply(diff_vals, 1, sd)
  diff_se   <- diff_sd / sqrt(n_sub)
  
  t_crit <- qt(1 - alpha / 2, df = n_sub - 1)
  
  ci_lower <- diff_mean - t_crit * diff_se
  ci_upper <- diff_mean + t_crit * diff_se
  sig <- (ci_lower > 0) | (ci_upper < 0)
  
  pre_mean  <- rowMeans(pre_vals)
  pre_sd    <- apply(pre_vals, 1, sd)
  post_mean <- rowMeans(post_vals)
  post_sd   <- apply(post_vals, 1, sd)
  
  results[[s]] <- list(
    pre_vals = pre_vals, post_vals = post_vals, diff_vals = diff_vals,
    pre_mean = pre_mean, pre_sd = pre_sd, post_mean = post_mean, post_sd = post_sd,
    diff_mean = diff_mean, diff_sd = diff_sd, diff_se = diff_se,
    ci_lower = ci_lower, ci_upper = ci_upper, sig = sig, t_crit = t_crit
  )
}

# --- Identify Significant Regions ---

find_sig_regions <- function(sig_vec, eval_pts) {
  if (!any(sig_vec)) return(NULL)
  
  rle_out <- rle(sig_vec)
  ends    <- cumsum(rle_out$lengths)
  starts  <- c(1, ends[-length(ends)] + 1)
  sig_idx <- which(rle_out$values == TRUE)
  
  data.frame(start_pct = eval_pts[starts[sig_idx]], end_pct = eval_pts[ends[sig_idx]])
}

# --- Results Summary ---

cat("\n========================================================================\n")
cat("RESULTS: Significant Regions (95% CI excludes zero)\n")
cat("========================================================================\n\n")

for (s in sheets) {
  r <- results[[s]]
  regions <- find_sig_regions(r$sig, eval_pts)
  
  cat(sprintf("--- %s ---\n", angle_labels[s]))
  
  if (is.null(regions) || nrow(regions) == 0) {
    cat("  No significant differences.\n\n")
    next
  }
  
  for (i in seq_len(nrow(regions))) {
    idx <- which(eval_pts >= regions$start_pct[i] & eval_pts <= regions$end_pct[i])
    peak_val <- r$diff_mean[idx][which.max(abs(r$diff_mean[idx]))]
    mean_val <- mean(r$diff_mean[idx])
    direction <- ifelse(mean_val > 0, "POST > PRE", "PRE > POST")
    
    cat(sprintf("  %3.0f%% – %3.0f%% | Mean: %+.2f° | Peak: %+.2f° | %s\n",
                regions$start_pct[i], regions$end_pct[i], mean_val, peak_val, direction))
  }
  cat("\n")
}

# --- Generate Figures ---

col_pre  <- "#2166AC"
col_post <- "#B2182B"
col_pre_rgb  <- col2rgb(col_pre) / 255
col_post_rgb <- col2rgb(col_post) / 255
col_diff <- "#333333"

# Sagittal Plane Figure
# pdf(file.path(output_dir, "sls_sagittal_plane.pdf"), width = 10, height = 12)

par(mfrow = c(3, 2), mar = c(4.5, 4.5, 2.5, 1), oma = c(0, 0, 2, 0))

for (s in c("ankle_x", "knee_x", "hip_x")) {
  r <- results[[s]]
  
  y_range <- range(c(r$pre_mean - r$pre_sd, r$pre_mean + r$pre_sd,
                     r$post_mean - r$post_sd, r$post_mean + r$post_sd))
  
  plot(eval_pts, r$pre_mean, type = "n", xlim = t_range, ylim = y_range,
       xlab = "SLS Cycle (%)", ylab = "Angle (°)", main = angle_labels[s])
  
  polygon(c(eval_pts, rev(eval_pts)), c(r$pre_mean + r$pre_sd, rev(r$pre_mean - r$pre_sd)),
          col = rgb(col_pre_rgb[1], col_pre_rgb[2], col_pre_rgb[3], 0.2), border = NA)
  polygon(c(eval_pts, rev(eval_pts)), c(r$post_mean + r$post_sd, rev(r$post_mean - r$post_sd)),
          col = rgb(col_post_rgb[1], col_post_rgb[2], col_post_rgb[3], 0.2), border = NA)
  
  lines(eval_pts, r$pre_mean, col = col_pre, lwd = 2.5)
  lines(eval_pts, r$post_mean, col = col_post, lwd = 2.5)
  legend("topright", legend = c("PRE", "POST"), col = c(col_pre, col_post), lwd = 2.5, bty = "n")
  
  # Difference plot
  y_range_diff <- range(c(r$ci_lower, r$ci_upper)) * 1.1
  
  plot(eval_pts, r$diff_mean, type = "n", xlim = t_range, ylim = y_range_diff,
       xlab = "SLS Cycle (%)", ylab = "Difference (°)", main = "POST – PRE ± 95% CI")
  abline(h = 0, lty = 2, col = "gray50")
  
  if (any(r$sig)) {
    regions <- find_sig_regions(r$sig, eval_pts)
    for (i in seq_len(nrow(regions))) {
      idx <- which(eval_pts >= regions$start_pct[i] & eval_pts <= regions$end_pct[i])
      polygon(c(eval_pts[idx], rev(eval_pts[idx])),
              c(rep(y_range_diff[1] * 1.05, length(idx)), rep(y_range_diff[2] * 1.05, length(idx))),
              col = rgb(1, 0.8, 0, 0.25), border = NA)
    }
  }
  
  polygon(c(eval_pts, rev(eval_pts)), c(r$ci_upper, rev(r$ci_lower)),
          col = rgb(0.5, 0.5, 0.5, 0.3), border = NA)
  lines(eval_pts, r$diff_mean, col = col_diff, lwd = 2.5)
}

mtext("Single Leg Squat - Sagittal Plane: PRE vs POST", outer = TRUE, cex = 1.2, font = 2)

# dev.off()

# Frontal Plane Figure
# pdf(file.path(output_dir, "sls_frontal_plane.pdf"), width = 10, height = 12)

par(mfrow = c(3, 2), mar = c(4.5, 4.5, 2.5, 1), oma = c(0, 0, 2, 0))

for (s in c("ankle_y", "knee_y", "hip_y")) {
  r <- results[[s]]
  
  y_range <- range(c(r$pre_mean - r$pre_sd, r$pre_mean + r$pre_sd,
                     r$post_mean - r$post_sd, r$post_mean + r$post_sd))
  
  plot(eval_pts, r$pre_mean, type = "n", xlim = t_range, ylim = y_range,
       xlab = "SLS Cycle (%)", ylab = "Angle (°)", main = angle_labels[s])
  
  polygon(c(eval_pts, rev(eval_pts)), c(r$pre_mean + r$pre_sd, rev(r$pre_mean - r$pre_sd)),
          col = rgb(col_pre_rgb[1], col_pre_rgb[2], col_pre_rgb[3], 0.2), border = NA)
  polygon(c(eval_pts, rev(eval_pts)), c(r$post_mean + r$post_sd, rev(r$post_mean - r$post_sd)),
          col = rgb(col_post_rgb[1], col_post_rgb[2], col_post_rgb[3], 0.2), border = NA)
  
  lines(eval_pts, r$pre_mean, col = col_pre, lwd = 2.5)
  lines(eval_pts, r$post_mean, col = col_post, lwd = 2.5)
  legend("topright", legend = c("PRE", "POST"), col = c(col_pre, col_post), lwd = 2.5, bty = "n")
  
  # Difference plot
  y_range_diff <- range(c(r$ci_lower, r$ci_upper)) * 1.1
  
  plot(eval_pts, r$diff_mean, type = "n", xlim = t_range, ylim = y_range_diff,
       xlab = "SLS Cycle (%)", ylab = "Difference (°)", main = "POST – PRE ± 95% CI")
  abline(h = 0, lty = 2, col = "gray50")
  
  if (any(r$sig)) {
    regions <- find_sig_regions(r$sig, eval_pts)
    for (i in seq_len(nrow(regions))) {
      idx <- which(eval_pts >= regions$start_pct[i] & eval_pts <= regions$end_pct[i])
      polygon(c(eval_pts[idx], rev(eval_pts[idx])),
              c(rep(y_range_diff[1] * 1.05, length(idx)), rep(y_range_diff[2] * 1.05, length(idx))),
              col = rgb(1, 0.8, 0, 0.25), border = NA)
    }
  }
  
  polygon(c(eval_pts, rev(eval_pts)), c(r$ci_upper, rev(r$ci_lower)),
          col = rgb(0.5, 0.5, 0.5, 0.3), border = NA)
  lines(eval_pts, r$diff_mean, col = col_diff, lwd = 2.5)
}

mtext("Single Leg Squat - Frontal Plane: PRE vs POST", outer = TRUE, cex = 1.2, font = 2)

# dev.off()

# --- Smoothing Diagnostics ---

cat("\n========================================================================\n")
cat("SMOOTHING DIAGNOSTICS\n")
cat("========================================================================\n\n")

gcv_vals <- smooth_pre[["knee_x"]]$gcv
cat(sprintf("Representative GCV (knee sagittal, PRE): Mean = %.6f\n", mean(gcv_vals, na.rm = TRUE)))
cat(sprintf("Basis: %d cubic B-splines | Lambda: %g | df = %d\n", n_basis, lambda_val, n_sub - 1))

cat("\n--- Analysis complete. ---\n")
