#' Internal Functions for Building the Stanvars Objects for the Custom Families
#'
#' @param family Name of the custom family for the Bayesian MSM. Currently only
#' supports "normal_ipw", "poisson_ipw", "zero_inflated_poisson_ipw".
#'
#' @param weight_args A named list containing, at a minimum, two vectors called
#' `lambda` and `delta` containing the location and scale of the stabilized
#' inverse probability weights respectively.
#'
#' @param weight_priors A named list containing priors for the regularized
#' weights to be for the two shape parameters of the beta distribution.
#'
#' @param multilevel Logical indicating whether the formula passed includes a
#' random effects structure. This is currently ignored but planned for future
#' development. See https://github.com/ajnafa/IPWBayes/issues/1 for the current
#' status
#'
#' @param ... Additional arguments reserved for future development
#'
#' @importFrom brms stanvar
#'
.make_stanvars <- function(family, weight_args, weight_priors, multilevel, ...) {

  # stanvars object for the data block
  data_vars <- .prep_weights_data(weight_args, multilevel, weight_priors)

  # stanvars object for the parameters and model block
  model_pars <- .prep_model_parameters(multilevel)

  # stanvars object for the functions block
  stan_functions <- .prep_stan_functions(family, multilevel)

  # Combine these into one
  msm_stanvars <- c(data_vars, model_pars, stan_functions)

  return(msm_stanvars)
}




# A function for building the stanvars object for the parameters blocks
.prep_model_parameters <- function(multilevel) {

  # Stanvar for the prior on the scale of the weights
  ipweights_par <- stanvar(
    scode = "vector<lower=0>[N] wts_z; // Standardized Latent IP Weights",
    block = "parameters"
  )

  # Stanvar for the regularized observation-level ip weights
  ipweights_tpar <- stanvar(
    scode = "
    // Compute the IPT Weights
    vector[N] w_tilde; // IPT Weights
    w_tilde = wts_lambda + wts_delta * wts_z[1];
    ",
    block = "tparameters"
  )

  # Stanvar for sampling the beta prior on the scale of the weights
  ipweights_prior_sampling <- stanvar(
    scode = "
    // Sampling the Weights
    wts_z ~ beta(delta_prior_alpha, delta_prior_beta);
    ",
    block = "model",
    position = "start"
  )

  # Add everything together
  ipweights_parameters <- c(ipweights_par, ipweights_tpar, ipweights_prior_sampling)

  # Return the combined stanvars
  return(ipweights_parameters)
}


# A function for building the stanvars object for the data block
.prep_weights_data <- function(weight_args, multilevel, weight_priors) {

  # Check that weights are passed as a named list
  stopifnot("Weights must be a named list containing vectors for at least the lambda and delta parameters" = is.list(weight_args))

  # Extract the data for lambda and delta
  lambda <- weight_args$lambda
  delta <- weight_args$delta

  # Extract the data for the weight priors
  delta_prior_alpha <- weight_priors[1]
  delta_prior_beta <- weight_priors[2]

  # Check that at least lambda and delta are present
  stopifnot(exprs = {
    !is.null(lambda) && !is.null(delta)
    all.equal(length(lambda), length(delta))
  })

  # Location of the inverse probability of treatment weights
  ipweights_lambda <- stanvar(
    x = lambda,
    scode = "// Statistics from the Design Stage Model
    vector<lower = 0>[N] wts_lambda; // Location of the Population-Level Weights",
    block = "data"
  )

  # Scale of the inverse probability of treatment weights
  ipweights_delta <- stanvar(
    x = delta,
    scode = "vector<lower = 0>[N] wts_delta; // Scale of the Population-Level Weights",
    block = "data"
  )

  # Shape prameter for the beta prior on the weights
  delta_prior_alpha <- stanvar(
    x = delta_prior_alpha,
    scode = "real<lower = 0> delta_prior_alpha;",
    block = "data"
  )

  # Shape prameter for the beta prior on the weights
  delta_prior_beta <- stanvar(
    x = delta_prior_beta,
    scode = "real<lower = 0> delta_prior_beta;",
    block = "data"
  )

  # Check whether to use group level weights
  if (multilevel) {

    # Print a message saying multilevel support is pending
    print("Support for weighted random effects is still under development")

  }

  # Add everything together
  ipw_weights_data <- c(ipweights_delta, ipweights_lambda, delta_prior_alpha, delta_prior_beta)

  # Return the combined stanvars object
  return(ipw_weights_data)
}


# A function for building the stanvars object for the functions block
.prep_stan_functions <- function(family, multilevel) {

  # Check that family argument is one of those supported
  stopifnot("Currently supported families are poisson_ipw, normal_ipw, and zero_inflated_poisson_ipw" = family %in% c("poisson_ipw", "normal_ipw", "zero_inflated_poisson_ipw"))


  if (family == "poisson_ipw") {

    "/* Weighted Log PMF of the Poisson Pseudo-Likelihood
    * Args:
    *   y: the response vector of length
    *   mu: the linear predictor
    *   w_tilde: the realized inverse probability weights
      * Returns:
      *   a scalar to be added to the log posterior
    */
    real poisson_ipw_lpmf(int y, vector mu, vector w_tilde){
      real weighted_term;
      weighted_term = 0.00;
      weighted_term = weighted_term + w_tilde * poisson_log_lpmf(y | mu);
      return weighted_term;
    }" -> stan_pseudo_likelihood
  }


  if (family == "normal_ipw") {

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
      }" -> stan_pseudo_likelihood

  }


  if (family == "zero_inflated_poisson_ipw") {

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
      }" -> stan_pseudo_likelihood
  }

  # stanvar object for the pseudo-likelihood
  stan_pseudo_likelihood <- stanvar(
    scode = stan_pseudo_likelihood,
    block = "functions"
  )

  # Check whether to use group level weights
  if (multilevel) {

    # Print a message saying multilevel support is pending
    print("Support for weighted random effects is still under development")

    # "/* Weighted Log PDF of the Gaussian Pseudo-Likelihood for Group Level Effects
    # * Args:
    #   *   y: the response vector of length N
    # *   mu: the linear predictor
    # *   sigma: noise parameter
    # *   w_tilde: the realized inverse probability weights
    # *   N: the number of observations
    # * Returns:
    #   *   a scalar to be added to the log posterior
    # */
    #   real normal_ipw_re_lpdf(vector y, vector mu, real sigma, vector w_tilde, int N) {
    #     real weighted_term;
    #     weighted_term = 0.00;
    #     for (n in 1:N) {
    #       weighted_term = weighted_term + w_tilde[n] * normal_lpdf(y[n] | mu[n], sigma);
    #     }
    #     return weighted_term;
    #   }" -> stan_weighted_random_effects
  }

  # Return the stanvar for the pseudo-likelihood
  return(stan_pseudo_likelihood)
}
