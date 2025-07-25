#An R6 object to store all the parameters for an instance of the model.

Parameters <- R6::R6Class("Parameters",
	public = list(
		growth_rates = NULL,
		alpha_values = NULL,
		error_variance = NULL,
		
		initialize = function(l) {
			#Takes a named list of all the parameters (or one-row data),
			#with the names as the column headings output by the C++ program
			l <- as.list(l) #In case we're given a one-row data frame.
			self$growth_rates <- unlist(l[grep("^r_", names(l))])
			if(length(self$growth_rates) == 0) {
				self$growth_rates <- exp(unlist(l[grep("^log_r_", names(l))]))
			}
			self$alpha_values <- matrix(
				unlist(l[grep("^alpha_", names(l))]),
				nrow = length(self$growth_rates),
				byrow = TRUE
			)
			self$error_variance <- l[["error_variance"]] #This may be NULL.
		}
	)
)
