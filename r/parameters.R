#An R6 object to store all the parameters for an instance of the model.

Parameters <- R6::R6Class("Parameters",
	public = list(
		growth_rates = NULL,
		alpha_values = NULL,
		error_variance = NULL,
		
		initialize = function(growth_rates, alpha_values, error_variance = NULL) {
			self$growth_rates <- growth_rates
			self$alpha_values <- alpha_values
			self$error_variance <- error_variance
		}
	)
)

parameters_weighted_mean <- function(parameters_list, weights) {
	growth_rates_list <- lapply(parameters_list, function(p){p$growth_rates})
	alpha_values_list <- lapply(parameters_list, function(p){p$alpha_values})
	error_variance_list <- lapply(parameters_list, function(p){p$error_variance})
	weighted_growth_rates_list <- lapply(1:length(weights), function(i){growth_rates_list[[i]] * weights[i]})
	weighted_alpha_values_list <- lapply(1:length(weights), function(i){alpha_values_list[[i]] * weights[i]})
	weighted_error_variance_list <- weights * error_variance_list
	growth_rates <- Reduce("+", weighted_growth_rates_list) / sum(weights)
	alpha_values <- Reduce("+", weighted_alpha_values_list) / sum(weights)
	error_variance <- sum(weighted_error_variance_list) / sum(weights)
	return(Parameters$new(growth_rates, alpha_values, error_variance))
}

parameters_mean <- function(parameters_list) {
	parameters_weighted_mean(parameters_list, rep(1, length(parameters_list)))
}
