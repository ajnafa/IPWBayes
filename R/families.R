#' Custom Bayesian Marginal Structural Model Families
#'
#' The \pkg{IPWBayes} package currently provides three model families
#' for estimating marginal structural models and is primarily a wrapper
#' around the \pkg{brms} package.
#'
#' @importFrom brms custom_family

# Guassian Pseudo-Likelihood
normal_ipw <- custom_family(
  name = "normal_ipw",
  dpars = c("mu", "sigma"),
  links = "identity",
  type = "real",
  lb = c(NA, 0),
  vars = c("w_tilde", "N"),
  loop = FALSE
)

# Poisson Pseudo-Likelihood
poisson_ipw <- custom_family(
  name = "poisson_ipw",
  dpars = "mu",
  links = "identity",
  type = "int",
  lb = 0,
  vars = c("w_tilde", "N"),
  loop = FALSE
)
