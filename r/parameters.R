#An R6 object to store all the parameters for an instance of the model.

Parameters <- R6::R6Class("Parameters",
	public = list(
		growth_rates = NULL,
		alpha_values = NULL,
		
		initialize = function(growth_rates, alpha_values) {
			self$growth_rates <- growth_rates
			self$alpha_values <- alpha_values
		}
	)
)
