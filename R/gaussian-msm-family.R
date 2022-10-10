#' Custom family for the Gaussian Marginal Structural Model
#'
#'

.gaussian_stanvars <- function(multilevel, ...) {

  ## Force multilevel option to FALSE, needs further testing
  multilevel <- FALSE

  # Code for the population-level gaussian pseudo-likelihood
  "/* Weighted Log PDF of the Gaussian Pseudo-Likelihood for Population Level Effects
   * Args:
   *   y: the response vector of length N
   *   mu: the linear predictor
   *   sigma: noise parameter
   *   w_tilde: the realized inverse probability weights
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real normal_ipw_lpdf(vector y, vector mu, real sigma, vector w_tilde) {
     real weighted_term;
     weighted_term = 0.00;
     weighted_term = weighted_term + w_tilde * normal_lpdf(y | mu, sigma);
     return weighted_term;
    }
  " -> stan_pseudo_gaussian_density

  # Code for the latent inverse probability of treatment weights
  "
  // Statistics from the Design Stage Model
  vector<lower = 0>[N] wts_lambda; // Location of the Population-Level Weights
  vector<lower = 0>[N] wts_delta; // Scale of the Population-Level Weights

  // Prior on the scale of the weights
  real<lower = 0> delta_prior_shape1;
  real<lower = 0> delta_prior_shape2;
  " -> latent_ipweights_data

  # Code for the latent scale of the weights
  latent_ipweights_par <- "vector<lower=0>[N] wts_z; // Standardized Latent IP Weights"

  # Code for the regularized observation-level ip weights
  "
  // Compute the IPT Weights
  vector[N] w_tilde; // IPT Weights
  w_tilde = wts_lambda + wts_delta * wts_z[1];
  "  -> latent_ipweights_tpar

  # Code for the beta prior on the scale of the weights
  "// Sampling the Weights
  wts_z ~ beta(delta_prior_shape1, delta_prior_shape2);
  " -> latent_ipweights_prior

  # If multilevel is true, use the double weighted approach from
  # Savitsky and Williams (2022) that re-weights the group level effects
  if (isTRUE(multilevel)) {

    "/* Weighted Log PDF of the Gaussian Pseudo-Likelihood for Group Level Effects
   * Args:
   *   y: the response vector of length N
   *   mu: the linear predictor
   *   sigma: noise parameter
   *   w_tilde: the realized inverse probability weights
   *   N: the number of observations
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real normal_ipw_re_lpdf(vector y, vector mu, real sigma, vector w_tilde, int N) {
     real weighted_term;
     weighted_term = 0.00;
     for (n in 1:N) {
       weighted_term = weighted_term + w_tilde[n] * normal_lpdf(y[n] | mu[n], sigma);
       }
     return weighted_term;
    }" -> stan_pseudo_gaussian_density_group

    # Code for the latent inverse probability of treatment weights
    "
    // Statistics from the Design Stage Model for the group level weights
    vector<lower = 0>[J_1] wts_upsilon; // Location of the Group-Level Weights
    vector<lower = 0>[J_1] wts_xi; // Scale of the Group-Level Weights

    // Prior on the scale of the weights
    real<lower = 0> sd_prior_shape1;
    real<lower = 0> sd_prior_shape2;
    " -> latent_group_ipweights_data

    # Code for the latent scale of the weights
    latent_group_ipweights_par <- "vector<lower=0>[N_1] wts_u; // Standardized Latent IP Weights"

    # Code for the regularized observation-level ip weights
    "
    // Compute the IP Weights for the Group-Level Effects
    vector[N_1] w_hat; // IPT Weights
    w_hat = wts_upsilon + wts_xi * wts_u[1];
    "  -> latent_group_ipweights_tpar

    # Code for the beta prior on the scale of the weights
    "// Sampling the Group-Level Latent Weights
    wts_u ~ beta(xi_prior_shape1, xi_prior_shape2);
    " -> latent_group_ipweights_prior

    # Code to actually resample the random effects
    # @todo figure out how to do this in the context of brms' non-centered
    # parameterization
  }
}

