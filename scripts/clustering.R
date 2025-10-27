library(tidyverse)
library(cluster)
library(plotly)
library(htmlwidgets)

# Task 1 ------------------------------------------------------------------

# data generator
generate_hypercube_clusters <- function(n, k, side_length, noise_sd = 1.0) {
  # n = number of dimensions (and clusters)
  # k = points per cluster
  # side_length = L distance of cluster centers along each axis
  # noise_sd = cluster spread
  
  # create centers at positive corners of nD hypercube
  centers <- diag(rep(side_length, n))  # (L,0,...), (0,L,...), etc.
  
  # generate points around each cluster center
  clusters <- lapply(1:n, function(i) {
    matrix(rnorm(k * n, mean = centers[i, ], sd = noise_sd), ncol = n)
  })
  
  # combine into one data frame
  data <- do.call(rbind, clusters)
  data <- as.data.frame(data)
  data$cluster_true <- rep(1:n, each = k)
  return(data)
}

# set up simulation
set.seed(611)
dims <- c(6, 5, 4, 3, 2)     
side_lengths <- seq(10, 1, -1)  
k <- 100                       
noise_sd <- 1.0                
results <- data.frame()

# run simulation
for (n in dims) {
  for (L in side_lengths) {
    cat("Running simulation for dimension", n, ", side length", side_lengths)
    dat <- generate_hypercube_clusters(n, k, L, noise_sd)
    # Run Gap Statistic with k-means
    gap <- clusGap(
      dat[, 1:n],
      FUN = kmeans,
      K.max = 10,
      nstart = 20,
      iter.max = 50,
      B = 50  # number of Monte Carlo samples for gap statistic
    )
    # estimate optimal k using the standard gap statistic rule
    best_k <- with(gap, maxSE(Tab[,"gap"], Tab[,"SE.sim"]))
    # Store result
    results <- rbind(results,
                     data.frame(dimension = n,
                                side_length = L,
                                estimated_clusters = best_k))
    
    cat(sprintf("n=%d, L=%.1f, estimated k=%d\n", n, L, best_k))
  }
}

# visualization
p <- ggplot(results, aes(x = side_length, y = estimated_clusters)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(size = 2) +
  geom_hline(aes(yintercept = dimension),
             linetype = "dashed", color = "red", size = 0.8) +
  scale_x_reverse(breaks = seq(10, 1, -1)) +
  facet_wrap(~ dimension, ncol = 3, scales = "free_y") +
  labs(title = "Gap Statistic vs. Cluster Separation by Dimension",
       subtitle = "Dashed red line shows the true number of clusters (n)",
       x = "Hypercube side length (L)",
       y = "Estimated clusters (Gap Statistic)") +
  theme_minimal(base_size = 14)
print(p)
dir.create("output", showWarnings = FALSE)
ggsave("output/task1_gap_plot.png", p, width = 8, height = 5)


# When side_length is large, clusters are far apart, so the gap statistic will
# usually find the correct number of clusters. As side_length decreases, clusters
# begin to overlap, so the estimated number of clusters drops. 
# n=2: the algorithm consistently begins to reduce the estimate after L=2.
# n=3: the algorithm consistently begins to reduce the estimate after L=3.
# n=4: the algorithm consistently begins to reduce the estimate after L=2,
# but did not consistently detect the correct number of clusters to begin with.
# n=5: the algorithm consistently begins to reduce the estimate after a slight
# increase at L=3, but did not consistently detect the correct number of clusters to start.
# n=6: the algorithm consistently begins to reduce the estimate after L=4,
# but did not consistently detect the correct number of clusters to begin with.


# Task 2 ------------------------------------------------------------------

# generate data
generate_shell_clusters <- function(n_shells, k_per_shell, max_radius, noise_sd = 0.1) {
  # n_shells: number of concentric shells
  # k_per_shell: points per shell
  # max_radius: outermost radius (inner radius will be spaced evenly)
  # noise_sd: radial noise (thickness)
  # returns: data.frame with columns x,y,z and shell (true radius)
  
  if(n_shells < 1) stop("n_shells must be >= 1")
  # avoid zero inner radius: start at a small positive value, e.g. max_radius / (n_shells + 0.5)
  radii <- seq(max_radius / (n_shells + 0), max_radius, length.out = n_shells)
  
  pts_list <- lapply(radii, function(r) {
    # sample directions uniformly on sphere: sample normal coordinates and normalize
    mat <- matrix(rnorm(k_per_shell * 3), ncol = 3)
    norms <- sqrt(rowSums(mat^2))
    dirs <- mat / norms
    # radial jitter: radius + gaussian noise on radial coordinate
    radii_noise <- r + rnorm(k_per_shell, mean = 0, sd = noise_sd)
    coords <- dirs * radii_noise
    df <- as.data.frame(coords)
    names(df) <- c("x", "y", "z")
    df$shell <- paste0("r", round(r, 3))
    df$true_radius <- r
    df
  })
  do.call(rbind, pts_list)
}


# well-separated shells (max_radius = 8)
dat_wide <- generate_shell_clusters(n_shells = 4, k_per_shell = 100, max_radius = 8, noise_sd = 0.1)
dat_wide$group <- factor(dat_wide$shell)   # categorical label for color

p_wide <- plot_ly(dat_wide, x = ~x, y = ~y, z = ~z, color = ~group, colors = "Set1",
                  type = "scatter3d", mode = "markers",
                  marker = list(size = 3, opacity = 0.8)) %>%
  layout(title = "Concentric shells (well-separated): max_radius = 8",
         scene = list(xaxis = list(title = "x"),
                      yaxis = list(title = "y"),
                      zaxis = list(title = "z")))
p_wide
saveWidget(p_wide, "output/shells_wide_maxradius8.html", selfcontained = TRUE)

# close shells (max_radius = 2)
dat_close <- generate_shell_clusters(n_shells = 4, k_per_shell = 100, max_radius = 2, noise_sd = 0.1)
dat_close$group <- factor(dat_close$shell)

p_close <- plot_ly(dat_close, x = ~x, y = ~y, z = ~z, color = ~group, colors = "Set1",
                   type = "scatter3d", mode = "markers",
                   marker = list(size = 3, opacity = 0.8)) %>%
  layout(title = "Concentric shells (close): max_radius = 2",
         scene = list(xaxis = list(title = "x"),
                      yaxis = list(title = "y"),
                      zaxis = list(title = "z")))

p_close
saveWidget(p_close, "output/shells_close_maxradius2.html", selfcontained = TRUE)

# spectral clustering function
spectral_kmeans <- function(x, k, d_threshold = 1, nstart = 20, iter.max = 50) {
  # x: data matrix or data.frame (rows = points)
  # k: desired number of clusters
  # d_threshold: adjacency threshold (connect i and j if dist < d_threshold)
  
  xmat <- as.matrix(x)
  n <- nrow(xmat)
  
  # pairwise distances and adjacency
  dmat <- as.matrix(dist(xmat))
  A <- (dmat < d_threshold) * 1L
  diag(A) <- 0L
  
  # degree and Laplacian
  deg <- rowSums(A)
  D <- diag(deg)
  L <- D - A
  
  # symmetric normalized Laplacian L_sym = D^{-1/2} L D^{-1/2}
  # handle zero degrees: set invsqrt to 0 for degree 0 nodes
  inv_sqrt_deg <- ifelse(deg > 0, 1 / sqrt(deg), 0)
  D_inv_sqrt <- diag(inv_sqrt_deg)
  L_sym <- D_inv_sqrt %*% L %*% D_inv_sqrt
  
  # numerical symmetry enforcement
  L_sym <- (L_sym + t(L_sym)) / 2
  
  # eigen-decomposition (symmetric)
  eig <- eigen(L_sym, symmetric = TRUE)
  vals <- eig$values
  vecs <- eig$vectors
  
  # eigen() sorts decreasing by default; we want k smallest eigenvalues
  idx_small <- order(vals)[1:k]
  U <- vecs[, idx_small, drop = FALSE]
  
  # row-normalize rows of U (common step)
  row_norms <- sqrt(rowSums(U^2))
  row_norms[row_norms == 0] <- 1  # avoid divide by 0
  U_norm <- U / row_norms
  
  # cluster the rows with kmeans in this spectral embedding
  km <- kmeans(U_norm, centers = k, nstart = nstart, iter.max = iter.max)
  
  # return km directly
  return(km)
}

# simulation
n_shells <- 4
k_per_shell <- 100
noise_sd <- 0.1
d_threshold <- 1 

# radii values to test: from 10 down to 0
max_radii <- seq(10, 0, by = -1)

# container for results
results <- data.frame()

B_sim <- 30

for (R in max_radii) {
  message("Testing max_radius = ", R)
  
  # generate data
  dat <- generate_shell_clusters(n_shells, k_per_shell, max_radius = R, noise_sd = noise_sd)
  coords <- dat[, c("x", "y", "z")]
  
  # define wrapper for clusGap
  FUNspectral <- function(x, k) {
    spectral_kmeans(x, k, d_threshold = d_threshold, nstart = 20, iter.max = 50)
  }
  
  # run clusGap
  gap <- tryCatch(
    clusGap(coords, FUN = FUNspectral, K.max = 8, B = B_sim),
    error = function(e) {
      message("clusGap failed at R=", R, " error: ", e$message)
      return(NULL)
    }
  )
  
  if (!is.null(gap)) {
    # estimate best k using maxSE rule
    best_k <- maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"])
  } else {
    best_k <- NA
  }
  
  results <- rbind(results, data.frame(max_radius = R, estimated_k = best_k))
  message("max_radius=", R, " -> estimated_k = ", best_k)
}

# visualization
results$max_radius <- as.numeric(results$max_radius)
p <- ggplot(results, aes(x = max_radius, y = estimated_k)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = n_shells, linetype = "dashed", color = "red") +
  scale_x_reverse() +        # so 10 -> 0 left -> right
  scale_y_continuous(breaks = seq(0, max(na.omit(results$estimated_k), n_shells + 2), by = 1)) +
  labs(title = "Spectral Clustering (d_threshold = 1): Estimated #clusters vs max_radius",
       subtitle = paste0("n_shells = ", n_shells, ", k_per_shell = ", k_per_shell,
                         ", noise_sd = ", noise_sd, ", B = ", B_sim),
       x = "max_radius (outermost shell radius)",
       y = "Estimated number of clusters (Gap Statistic)") +
  theme_minimal(base_size = 13)
print(p)
print(results)
ggsave("output/task2_spectral_gap_plot.png", p, width = 8, height = 5, dpi = 150)

# Looking at the plot, the algorithm fails to detect the correct number of clusters. As max_radius decreases, 
# the shells move closer to each other. After max_radius=7.5, the failure point, estimated_k drops below 4 
# consistently which indicates spectral clustering no longer identifies 4 distinct components.
# A smaller threshold, 0.8 for example, would result in fewer edges, improving separation. But if too small, 
# within-shell points could be disconnected, causing fragmentation. A larger threshold would result in more 
# edges, so more inter-shell connections. Shells will merge earlier as the radius shrinks, so the algorithm 
# will fail at a larger max_radius.