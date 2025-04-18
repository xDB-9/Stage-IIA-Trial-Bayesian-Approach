# Load necessary libraries
library(ggplot2)
library(dplyr)
library(pracma)  # For numerical integration
library(gridExtra)  # For arranging plots side by side

# Parameters
n_max <- 30
interim_points <- c(10, 20, 30)
target_rate <- 0.4  # p1
null_rate <- 0.2    # p0
a0 <- 0.5          # Prior parameters
b0 <- 0.5
theta_u <- 0.9     # Efficacy threshold
theta_l <- 0.1     # Futility threshold
ev_threshold <- 0.5  # Evidence threshold for BEV

# Simulate data
set.seed(123)
responses <- rbinom(n_max, 1, 0.4)  # True response rate is 40%

# Calculate BEV at interim points
bev_values <- sapply(1:length(interim_points), function(i) {
  n <- interim_points[i]
  x <- sum(responses[1:n])
  posterior <- function(p) dbeta(p, a0 + x, b0 + n - x)
  p_seq <- seq(0, 1, length.out = 500)
  posterior_vals <- posterior(p_seq)
  info <- posterior_vals / dbeta(p_seq, a0, b0)  # Using prior as reference r(p)
  ei <- p_seq[info >= ev_threshold]
  if (length(ei) > 0) {
    h1_region <- ei[ei > null_rate]
    if (length(h1_region) > 0) {
      integral <- trapz(h1_region, posterior(h1_region))
      return(integral / trapz(p_seq, posterior(p_seq)))  # Normalize BEV
    }
  }
  return(0)  # If no evidence region, BEV is 0
})

# Function to create plot for a specific interim point
create_bev_plot <- function(n, x, bev) {
  plot_data <- data.frame(p = seq(0, 1, length.out = 500))
  plot_data$Prior <- dbeta(plot_data$p, a0, b0)
  
  # Add posterior curve for this n
  label <- paste0("Posterior (n=", n, ", x=", x, ")")
  plot_data[[label]] <- dbeta(plot_data$p, a0 + x, b0 + n - x)
  
  plot_data_long <- plot_data %>%
    pivot_longer(cols = -p, names_to = "distribution", values_to = "density") %>%
    filter(distribution %in% c("Prior", label))  # Ensure only two levels
  
  plot_data_long$distribution <- factor(
    plot_data_long$distribution,
    levels = c("Prior", label)
  )
  plot_data_long$distribution <- droplevels(plot_data_long$distribution)  # Drop unused levels
  
  # Shade EI for H1
  bev_shade <- function(n, x, threshold) {
    p_seq <- seq(0, 1, length.out = 500)
    post <- dbeta(p_seq, a0 + x, b0 + n - x)
    info <- post / dbeta(p_seq, a0, b0)
    ei <- p_seq[info >= threshold & p_seq > null_rate]
    if (length(ei) > 0) {
      data.frame(p = ei, density = post[match(ei, p_seq)])
    } else {
      data.frame(p = numeric(0), density = numeric(0))
    }
  }
  
  shade_data <- bev_shade(n, x, ev_threshold)
  
  ggplot(plot_data_long, aes(x = p, y = density, color = distribution, linetype = distribution)) +
    geom_line(size = 1.2) +
    geom_vline(xintercept = c(null_rate, target_rate), linetype = "dashed", color = "darkgray") +
    annotate("text", x = null_rate, y = 0, label = paste("p0 =", null_rate), hjust = -0.1, vjust = -1, color = "darkred") +
    annotate("text", x = target_rate, y = 0, label = paste("p1 =", target_rate), hjust = -0.1, vjust = -1, color = "darkgreen") +
    geom_area(data = shade_data, aes(x = p, y = density), fill = "lightblue", alpha = 0.4, color = NA, inherit.aes = FALSE) +
    geom_label(aes(x = 0.85, y = max(plot_data_long$density) * 0.8), 
               label = paste("BEV(H1) =", round(bev, 3)), vjust = 1,
               hjust = 1, fill = "white", color = "black", fontface = "bold", size = 5) +
    scale_color_manual(values = c("gray60", "navy")) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    labs(x = "Response Rate (p)",
         y = "Probability Density",
         title = paste("BEV at n =", n, ", x =", x),
         subtitle = "Shaded = EI for H1; Dotted lines = thresholds") +
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(face = "bold"))
}

# Generate plots for each interim point
p1 <- create_bev_plot(interim_points[1], sum(responses[1:interim_points[1]]), bev_values[1])
p2 <- create_bev_plot(interim_points[2], sum(responses[1:interim_points[2]]), bev_values[2])
p3 <- create_bev_plot(interim_points[3], sum(responses[1:interim_points[3]]), bev_values[3])

# Arrange plots side by side
grid.arrange(p1, p2, p3, ncol = 3)

# Print BEV values for confirmation
print("BEV values at interim points:")
print(round(bev_values, 3))
