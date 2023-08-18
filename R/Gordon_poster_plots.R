model_filter <- c("DM recruitment by PP", "DM Static", "Non-linear 1", "Exponential")


# 270025
# 320878
# 355773
# 113443
# 223936
# 328824

n_plot |>
  ggplot() +
  aes(x = timestep_idx, y = observed) +
  geom_point() +
  geom_line() +
  facet_wrap(~agrp_prp_id, scales = "free") +
  theme_bw()

title_size <- 18
axis_size <- 18
strip_size <- 14


g1 <- n_plot |>
  filter(agrp_prp_id %in% c(270025),
         model %in% model_filter) |>
  mutate(observed = log(observed + 1)) |>
  ggplot() +
  aes(x = timestep_idx, y = median, ymin = lower95, ymax = upper95) +
  geom_ribbon(fill = "lightblue") +
  geom_point(aes(y = observed, shape = "Pigs removed")) +
  geom_line() +
  labs(title = "Figure 1: Latent pig abundance",
       y = "Log(abundnace + 1)",
       x = "Timestep",
       shape = "",
       fill = "Model") +
  facet_grid( ~ model) +
  theme_bw() +
  theme(legend.position = "none",
        title = element_text(size = title_size),
        axis.title = element_text(size = axis_size),
        axis.text = element_text(size = axis_size-5),
        strip.text = element_text(size = strip_size))
g1
ggsave("latent_abundance.jpeg",
       g1,
       device = "jpeg",
       dpi = "print")

g2 <- all_po |>
  filter(agrp_prp_id %in% c(270025),
         model %in% model_filter) |>
  mutate(observed = log(observed + 1)) |>
  group_by(model) |>
  mutate(x = 1:n()) |>
  ggplot() +
  aes(x = x, y = median, ymin = lower95, ymax = upper95) +
  geom_linerange(color = "lightblue") +
  geom_point(aes(y = observed, pch = "Removed"), color = "black") +
  labs(title = "Figure 2: Predicted take",
       y = "Log(take + 1)",
       x = "Removal event",
       shape = "",
       color = "Model") +
  # scale_color_discrete(values = "black") +
  facet_grid( ~ model) +
  theme_bw() +
  theme(legend.position = "none",
        title = element_text(size = title_size),
        axis.title = element_text(size = axis_size),
        axis.text = element_text(size = axis_size-5),
        strip.text = element_text(size = strip_size))

g2
ggsave("predicted_take.jpeg",
       g2,
       device = "jpeg",
       dpi = "print")

g3 <- all_po |>
  filter(agrp_prp_id %in% c(270025),
    model %in% model_filter) |>
  mutate(rmse = log(rmse)) |>
  group_by(model) |>
  mutate(x = 1:n()) |>
  ggplot() +
  aes(x = x, y = rmse) +
  geom_point() +
  labs(title = "Figure 3: Model performance",
       y = "Log(RMSE)",
       x = "Removal event",
       shape = "",
       color = "Model") +
  facet_grid(~ model) +
  theme_bw() +
  theme(legend.position = "none",
        title = element_text(size = title_size),
        axis.title = element_text(size = axis_size),
        axis.text = element_text(size = axis_size-5),
        strip.text = element_text(size = strip_size))
g3
ggsave("rmse.jpeg",
       g3,
       device = "jpeg",
       dpi = "print")

g4 <- all_params |>
  filter(model %in% model_filter) |>
  mutate(parameter = if_else(grepl("ruggedness", parameter), "ruggedness", parameter),
         parameter = if_else(grepl("road density", parameter), "road density", parameter),
         parameter = if_else(grepl("canopy cover", parameter), "canopy cover", parameter)) |>
  gg_obs_slopes() +
  labs(color = "Model",
       title = "Figure 7: Effect of terrain on take probability") +
  theme(#legend.position = "none",
        title = element_text(size = title_size),
        axis.title = element_text(size = axis_size),
        axis.text = element_text(size = axis_size-5),
        strip.text = element_text(size = strip_size))

g4
ggsave("obs_model.jpeg",
       g4,
       device = "jpeg",
       dpi = "print")


g5 <- all_params |>
  filter(model %in% model_filter) |>
  filter(parameter == "Population growth") |>
  ggplot() +
  aes(x = model) +
  geom_point(aes(y = median), size = 6, position = position_dodge(width=0.5)) +
  geom_linerange(aes(ymin = lower95, ymax = upper95), linewidth = 1, position = position_dodge(width=0.5)) +
  geom_linerange(aes(ymin = lower75, ymax = upper75), linewidth = 3, position = position_dodge(width=0.5)) +
  theme_bw() +
  labs(title = "Figure 4: Instantaneous growth rate",
       y = "Posterior estimate",
       x = "Model") +
  theme(legend.position = "none",
        title = element_text(size = title_size-5),
        axis.title = element_text(size = axis_size),
        axis.text = element_text(size = axis_size-5),
        strip.text = element_text(size = strip_size))
g5
ggsave("lambda.jpeg",
       g5,
       width = 5,
       height = 5,
       units = "in",
       device = "jpeg",
       dpi = "print")


g6 <- all_params |>
  filter(model %in% model_filter) |>
  filter(grepl("phi", node)) |>
  ggplot() +
  aes(x = model) +
  geom_point(aes(y = median), size = 6, position = position_dodge(width=0.5)) +
  geom_linerange(aes(ymin = lower95, ymax = upper95), linewidth = 1, position = position_dodge(width=0.5)) +
  geom_linerange(aes(ymin = lower75, ymax = upper75), linewidth = 3, position = position_dodge(width=0.5)) +
  theme_bw() +
  labs(title = "Figure 5: Apparent survival rate",
       y = "",
       x = "Model") +
  theme(legend.position = "none",
        title = element_text(size = title_size-5),
        axis.title = element_text(size = axis_size),
        axis.text = element_text(size = axis_size-5),
        strip.text = element_text(size = strip_size))

g6
ggsave("phi.jpeg",
       g6,
       width = 5,
       height = 5,
       units = "in",
       device = "jpeg",
       dpi = "print")

every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}

g6 <- all_params |>
  filter(model %in% "DM recruitment by PP") |>
  filter(grepl("log_zeta[", node, fixed = TRUE)) |>
  # arrange(parameter)
  mutate(parameter = as.numeric(parameter)) |>
  ggplot() +
  aes(x = parameter) +
  geom_point(aes(y = median), size = 2, position = position_dodge(width=0.5)) +
  geom_linerange(aes(ymin = lower95, ymax = upper95), linewidth = 0.5, position = position_dodge(width=0.5)) +
  # geom_linerange(aes(ymin = lower75, ymax = upper75), linewidth = 1, position = position_dodge(width=0.5)) +
  theme_bw() +
  # scale_x_discrete(breaks = every_nth(n = 10)) +
  labs(title = "Figure 6: Per capita recruitment",
       y = "Posterior estimate",
       x = "Primary period") +
  theme(legend.position = "none",
        title = element_text(size = title_size),
        axis.title = element_text(size = axis_size),
        axis.text = element_text(size = axis_size-5),
        strip.text = element_text(size = strip_size))

g6
ggsave("zeta_pp.jpeg",
       g6,
       device = "jpeg",
       dpi = "print")

g7 <- all_params |>
  filter(model %in% model_filter) |>
  filter(grepl("capita recruitment", parameter),
         node != "tau_zeta") |>
  # arrange(parameter)
  mutate(parameter = as.numeric(parameter)) |>
  ggplot() +
  aes(x = parameter) +
  geom_point(aes(y = median), size = 2, position = position_dodge(width=0.5)) +
  geom_linerange(aes(ymin = lower95, ymax = upper95), linewidth = 0.5, position = position_dodge(width=0.5)) +
  # geom_linerange(aes(ymin = lower75, ymax = upper75), linewidth = 1, position = position_dodge(width=0.5)) +
  theme_bw() +
  # scale_x_discrete(breaks = every_nth(n = 10)) +
  labs(title = "Figure 6: Per capita recruitment",
       y = "Posterior estimate",
       x = "Primary period") +
  theme(legend.position = "none",
        title = element_text(size = title_size),
        axis.title = element_text(size = axis_size),
        axis.text = element_text(size = axis_size-5),
        strip.text = element_text(size = strip_size))

g7
ggsave("zeta.jpeg",
       g6,
       device = "jpeg",
       dpi = "print")
