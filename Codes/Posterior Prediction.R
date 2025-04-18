##################################
#Line representation
##################################
library(ggplot2)
library(dplyr)

# Parameters
n_max <- 30
interim_points <- c(10, 20, 30)
target_rate <- 0.4  # p1
null_rate <- 0.2    # p0

# Simulate data
set.seed(123)
responses <- rbinom(n_max, 1, 0.4)  # True response rate is 40%

# Create patient data
patient_data <- data.frame(
  patient_id = 1:n_max,
  response = responses,
  interim = ifelse(1:n_max %in% interim_points, 1, 0)
)

# Calculate running statistics
running_stats <- patient_data %>%
  mutate(
    cum_response = cumsum(response),
    response_rate = cum_response / patient_id,
    lower_ci = qbeta(0.025, cum_response + 0.5, patient_id - cum_response + 0.5),
    upper_ci = qbeta(0.975, cum_response + 0.5, patient_id - cum_response + 0.5)
  )

# Calculate dynamic posterior probabilities at interim points
pp_values <- sapply(interim_points, function(n) {
  x <- sum(responses[1:n])
  1 - pbeta(target_rate, x + 0.5, n - x + 0.5)  # P(p > target_rate | data)
})

# Create data for annotation of posterior probabilities
pp_data <- data.frame(
  x = interim_points,
  y = running_stats$response_rate[interim_points],
  pp = round(pp_values, 2)  # Rounded for display
)

# Plot
ggplot(running_stats) +
  # Patient dots
  geom_point(aes(x = patient_id, y = 0, color = factor(response)), 
             size = 8, shape = 16) +
  
  # Running response rate
  geom_line(aes(x = patient_id, y = response_rate), 
            linewidth = 1, color = "blue") +
  
  # Confidence interval ribbon
  geom_ribbon(aes(x = patient_id, ymin = lower_ci, ymax = upper_ci),
              alpha = 0.2, fill = "blue") +
  
  # Interim analysis vertical lines
  geom_vline(xintercept = interim_points, 
             linetype = "dashed", color = "darkgray") +
  
  # Posterior probability annotations
  geom_text(data = pp_data,
            aes(x = x, y = y + 0.1, label = paste("PP =", pp)),
            color = "purple", fontface = "bold") +
  
  # Target and null reference lines
  geom_hline(yintercept = target_rate, linetype = "dotted", color = "darkgreen") +
  annotate("text", x = 1, y = target_rate, label = "Target p1 = 0.4", 
           hjust = 0, vjust = -0.5, color = "darkgreen") +
  
  geom_hline(yintercept = null_rate, linetype = "dotted", color = "darkred") +
  annotate("text", x = 1, y = null_rate, label = "Null p0 = 0.2", 
           hjust = 0, vjust = -0.5, color = "darkred") +
  
  # Appearance and theming
  scale_color_manual(values = c("tomato", "forestgreen"),
                     labels = c("No Response", "Response"),
                     name = "Outcome") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(-0.1, 1)) +
  labs(x = "Patient Number", 
       y = "Response Rate",
       title = "Patient Enrollment Timeline with Interim Analyses",
       subtitle = "Running response rate with 95% credible interval and Bayesian PP") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")



#######################################
#Distribution plot
#######################################

set.seed(123)
n_max <- 30
interim_points <- c(10, 20, 30)
true_rate <- 0.4
responses <- rbinom(n_max, 1, true_rate)

# Prior parameters
a0 <- 0.5
b0 <- 0.5
p0 <- 0.2
p1 <- 0.4
theta_u <- 0.9
theta_l <- 0.1

# Cumulative responses at interim
cum_responses <- cumsum(responses)
response_obs <- cum_responses[interim_points]

# Posterior probability at each interim
pp_values <- sapply(1:length(interim_points), function(i) {
  n <- interim_points[i]
  x <- response_obs[i]
  1 - pbeta(p1, a0 + x, b0 + n - x)
})

# Determine stopping point
stop_index <- which(pp_values > theta_u | pp_values < theta_l)[1]
if (is.na(stop_index)) stop_index <- length(interim_points)
stop_n <- interim_points[stop_index]
stop_x <- response_obs[stop_index]
stop_pp <- round(pp_values[stop_index], 3)

# Create posterior distribution data
p_seq <- seq(0, 1, length.out = 500)
plot_data <- data.frame(p = p_seq)
plot_data$Prior <- dbeta(p_seq, a0, b0)

# Add posterior curves
for (i in 1:stop_index) {
  n <- interim_points[i]
  x <- response_obs[i]
  label <- paste0("Posterior (n=", n, ", x=", x, ")")
  plot_data[[label]] <- dbeta(p_seq, a0 + x, b0 + n - x)
}

# Convert to long format
plot_data_long <- plot_data %>%
  pivot_longer(cols = -p, names_to = "distribution", values_to = "density")

# Set factor levels for ordering
plot_data_long$distribution <- factor(
  plot_data_long$distribution,
  levels = c("Prior", paste0("Posterior (n=", interim_points[1:stop_index], ", x=", response_obs[1:stop_index], ")"))
)

# Posterior means for annotation
posterior_lines <- data.frame()
for (i in 1:stop_index) {
  n <- interim_points[i]
  x <- response_obs[i]
  mean_p <- (a0 + x) / (a0 + b0 + n)
  label <- paste0("n=", n)
  distribution <- paste0("Posterior (n=", n, ", x=", x, ")")
  posterior_lines <- rbind(posterior_lines, data.frame(
    p = mean_p, label = label, distribution = distribution
  ))
}

# Max density for scaling text position
max_y <- max(plot_data_long$density)
annot_label <- paste0("Stopped at n = ", stop_n, "\nPP(p > ", p1, ") = ", stop_pp)

# Plot
ggplot(plot_data_long, aes(x = p, y = density, color = distribution, linetype = distribution)) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = c(p0, p1), linetype = "dashed", color = "darkgray") +
  annotate("text", x = p0, y = 0, label = paste("p0 =", p0), hjust = -0.1, vjust = -1, color = "darkred") +
  annotate("text", x = p1, y = 0, label = paste("p1 =", p1), hjust = -0.1, vjust = -1, color = "darkgreen") +
  annotate("rect", xmin = p1, xmax = 1, ymin = 0, ymax = max_y, alpha = 0.1, fill = "green") +
  
  # Add posterior mean lines and readable labels
  geom_vline(data = posterior_lines, aes(xintercept = p, color = distribution), linetype = "dotted", linewidth = 1) +
  geom_text(data = posterior_lines, aes(x = p, y = max_y * 0.85, label = label, color = distribution),
            angle = 90, vjust = -0.3, hjust = 0, size = 4, fontface = "bold") +
  
  annotate("text", x = 0.7, y = max_y * 0.95, label = annot_label, hjust = 0, color = "black", fontface = "bold") +
  
  scale_color_manual(values = c("gray60", "skyblue", "royalblue", "navy")) +
  scale_linetype_manual(values = c("dashed", "solid", "solid", "solid")) +
  labs(x = "Response Rate (p)",
       y = "Probability Density",
       title = "Posterior Distributions with Interim Analyses",
       subtitle = "Shaded = efficacy zone; Dotted lines = posterior means with labels") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))