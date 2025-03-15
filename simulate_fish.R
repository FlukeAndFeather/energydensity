library(tidyverse)

set.seed(123)

# Sample size
n <- 120

# Species and their size ranges
species <- c("anchovy", "clubhook", "lanternfish", "market")
size_range <- list(anchovy = c(5, 40),
                   clubhook = c(4, 50),
                   lanternfish = c(2, 20),
                   market = c(5, 80))

# Coefficients for simulation

# ed ("true" energy density)

# This specifies ed is Gamma-distributed with a log link to the linear
# predictor. Contrast this with Gaussian regression on log-transformed ed. X is
# the design matrix (intercept, species, mass, mass * species; anchovy is the
# reference level)

# ed ~ Gamma(shape, rate)
# shape = inverse_phi
# rate = inverse_phi / mu
# mu = exp(X * betas)
# betas[1] ~ Cauchy(0, 10)
# betas[2:] ~ Cauchy(0, 2.5)
# inverse_phi ~ Exponential(1)

ed_betas <- c(intercept = 2,
              clubhook = -1.8,
              lanternfish = -1.5,
              market = -1.2,
              mass = -0.008,
              mass_clubhook = 0.024,
              mass_lanternfish = 0.146,
              mass_market = 0.0165)
ed_invphi = 10

expected_ed <- function(s, m) {
  exp(model.matrix(~ s * m) %*% ed_betas)
}
simulate_ed <- function(s, m) {
  mu <- expected_ed(s, m)
  rgamma(n = length(s), shape = ed_invphi, rate = ed_invphi / mu)
}

# Simulate data
fish <- tibble(
  species = sample(species, size = n, replace = TRUE),
  mass_g = runif(n,
                 min = map_dbl(species, \(s) size_range[[s]][1]),
                 max = map_dbl(species, \(s) size_range[[s]][2])),
  ed_kjg = simulate_ed(species, mass_g)
)

fish_expected <- expand_grid(
  species = species,
  mass_frac = seq(0, 1, length.out = 100)
) %>%
  mutate(
    mass_g = map2_dbl(
      species,
      mass_frac,
      \(s, f) (size_range[[s]][2] - size_range[[s]][1]) * f + size_range[[s]][1]
    ),
    ed_kjg = expected_ed(species, mass_g)
  )
ggplot(fish, aes(mass_g, ed_kjg, color = species)) +
  geom_line(data = fish_expected, linewidth = 1.5) +
  geom_point(shape = 21, size = 2.5) +
  labs(x = "Mass (g)",
       y = "ED (kJ g^-1)") +
  theme_classic() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.99, 0.99),
        legend.justification = c(1, 1))
