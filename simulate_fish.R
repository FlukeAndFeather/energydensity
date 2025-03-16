library(brms)
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
# the design matrix (intercept, species, region, mass, mass * species, mass *
# region; anchovy and north are the reference levels for species and region.)

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
              south = 0.1,
              mass_clubhook = 0.024,
              mass_lanternfish = 0.146,
              mass_market = 0.0165,
              mass_south = 0.004)
ed_invphi = 10

expected_ed <- function(s, r, m) {
  exp(model.matrix(~ s * m + r * m) %*% ed_betas)
}
simulate_ed <- function(s, r, m) {
  mu <- expected_ed(s, r, m)
  rgamma(n = length(s), shape = ed_invphi, rate = ed_invphi / mu)
}

# Simulate data
fish <- tibble(
  species = sample(species, size = n, replace = TRUE),
  region = sample(c("north", "south"), size = n, replace = TRUE),
  mass_g = runif(n,
                 min = map_dbl(species, \(s) size_range[[s]][1]),
                 max = map_dbl(species, \(s) size_range[[s]][2])),
  ed_kjg = simulate_ed(species, region, mass_g)
)

fish_expected <- expand_grid(
  species = species,
  region = c("north", "south"),
  mass_frac = seq(0, 1, length.out = 100)
) %>%
  mutate(
    mass_g = map2_dbl(
      species,
      mass_frac,
      \(s, f) (size_range[[s]][2] - size_range[[s]][1]) * f + size_range[[s]][1]
    ),
    ed_kjg = expected_ed(species, region, mass_g)
  )

ggplot(fish, aes(mass_g, ed_kjg, color = species)) +
  geom_line(aes(linetype = region), fish_expected, linewidth = 1.5) +
  geom_point(aes(shape = region), size = 2.5) +
  labs(x = "Mass (g)",
       y = "ED (kJ g^-1)") +
  scale_shape_manual(values = c(19, 21)) +
  scale_linetype_manual(values = c("dotted", "dashed")) +
  theme_classic() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.99, 0.99),
        legend.justification = c(1, 1))

ed_mod <- brm(ed_kjg ~ mass_g * species + mass_g * region,
              fish,
              family = Gamma(link = "log"),
              prior = c(
                prior(cauchy(0, 10), class = "Intercept"),
                prior(cauchy(0, 5), class = "b"),
                prior(exponential(1), class = "shape")
              ),
              cores = 4)

sim2brm <- tribble(
  ~sim_param,         ~brm_param,
  "intercept",        "b_Intercept",
  "clubhook",         "b_speciesclubhook",
  "lanternfish",      "b_specieslanternfish",
  "market",           "b_speciesmarket",
  "mass",             "b_mass_g",
  "south",            "b_regionsouth",
  "mass_clubhook",    "b_mass_g:speciesclubhook",
  "mass_lanternfish", "b_mass_g:specieslanternfish",
  "mass_market",      "b_mass_g:speciesmarket",
  "mass_south",       "b_mass_g:regionsouth"
)
hdis <- bayestestR::hdi(ed_mod, ci = 0.95)
maps <- bayestestR::map_estimate(ed_mod)
tibble(
  param = names(ed_betas),
  value = ed_betas
) %>%
  left_join(sim2brm, by = c(param = "sim_param")) %>%
  mutate(hdi_lo = map_dbl(brm_param, \(p) hdis$CI_low[hdis$Parameter == p]),
         hdi_hi = map_dbl(brm_param, \(p) hdis$CI_high[hdis$Parameter == p]),
         map = map_dbl(brm_param, \(p) maps$MAP_Estimate[maps$Parameter == p])) %>%
  ggplot(aes(y = param)) +
  geom_pointrange(aes(x = map, xmin = hdi_lo, xmax = hdi_hi),
                  color = "cornflowerblue",
                  linewidth = 1.5) +
  geom_point(aes(x = value), color = "firebrick") +
  theme_bw()


