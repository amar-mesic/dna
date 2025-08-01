#' @title Sample log normal model(s) with parameters according to priors
#'
#' @param number_of_contributors Integer
#' @param sampling_parameters List. Needs to contain:
#' \itemize{
#'  \item min_template. Numeric of length one.
#'  \item max_template. Numeric of length one.
#'  \item degradation_shape. Numeric of length one.
#'  \item degradation_scale. Numeric of length one.
#' }
#' @param model_settings List. See \link{log_normal_model}.
#' @details In simulation studies involving many mixed DNA profiles, one often needs to generate various samples with different model parameters. This function samples a log normal model with parameters according to prior distributions. The template parameter for each contributor is sampled uniformly between \code{min_template} and \code{max_template}. The degradation parameter for each contributor is sampled from a gamma distribution with parameters \code{degradation_shape} and \code{degradation_scale}.
#' @return When \code{length(number_of_contributors)==1}, a single \link{log_normal_model} of class \code{pg_model}. Otherwise, a list of these.
#' @examples
#' gf <- gf_configuration()
#'
#' sampling_parameters <- list(min_template = 50., max_template = 1000.,
#'                             degradation_shape = 2.5, degradation_scale = 1e-3)
#'
#' model_no_stutter <- sample_log_normal_model(number_of_contributors = 1,
#'                                             sampling_parameters = sampling_parameters,
#'                                             model_settings = gf$log_normal_settings)
#' @export
sample_log_normal_model <- function(number_of_contributors, sampling_parameters, model_settings){


  params <- sample_log_normal_parameters_amar(number_of_contributors, sampling_parameters, model_settings)

  if (is.null(model_settings$stutter_model)){
    model <- log_normal_model(template = params$template, degradation = params$degradation,
                              LSAE = params$LSAE, c2 = params$c2,
                              model_settings = model_settings)
  }
  else{
    k2 <- sample_log_normal_stutter_variance(model_settings$stutter_variability)

    model <- log_normal_model(template = params$template, degradation = params$degradation,
                              LSAE = params$LSAE, c2 = params$c2,
                              k2 = params$k2,
                              model_settings = model_settings)
  }

  model
}
