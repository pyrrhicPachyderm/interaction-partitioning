#An R6 object to read, store and handle the raw data that was fed into the C++ program.
#Requires the three data files that were provided to the C++ object:
#the focal vector, the response vector, and the design matrix.

InputData <- R6::R6Class("InputData",
	public = list(
		is_per_capita = NULL,
		focal_vector = NULL,
		response_vector = NULL,
		design_matrix = NULL,
		num_species = NULL,
		num_obs = NULL,
		
		initialize = function(data_type, focal_vector_file_name, response_vector_file_name, design_matrix_file_name) {
			if(data_type == "indv") self$is_per_capita = TRUE
			else if(data_type == "pop") self$is_per_capita = FALSE
			else stop("Unrecognised input data type.")
			
			self$focal_vector <- read.table(focal_vector_file_name)[[1]] + 1 #Add one to make it 1-indexed, as R prefers.
			self$response_vector <- read.table(response_vector_file_name)[[1]]
			self$design_matrix <- as.matrix(read.table(design_matrix_file_name))
			
			self$num_species <- ncol(self$design_matrix)
			self$num_obs <- nrow(self$design_matrix)
		},
		
		get_fitted_values = function(parameters) {
			focal_growth_rate <- parameters$growth_rates[self$focal_vector]
			if(!self$is_per_capita) {
				focal_density <- mapply(
					function(obs,focal){self$design_matrix[obs,focal]},
					1:num_obs, self$focal_vector
				)
			} else {
				focal_density <- rep(1, num_obs)
			}
			alpha_values_on_focal <- parameters$alpha_values[self$focal_vector,]
			
			intrinsic_growth <- focal_growth_rate
			total_competition <- rowSums(alpha_values_on_focal * self$design_matrix)
			fitted_values <- intrinsic_growth * focal_density * (1 - total_competition)
			
			return(fitted_values)
		},
		
		get_residuals = function(parameters) {
			self$get_fitted_values(parameters) - self$response_vector
		},
		
		get_partial_residuals = function(parameters, row_index, col_index) {
			#The partial residuals for a particular alpha value: alpha_{row_index,col_index}.
			residuals <- self$get_residuals(parameters)
			partial_residuals_x <- self$get_partial_residuals_x(row_index, col_index)
			partial_residuals <- residuals[self$focal_vector == row_index] + parameters$alpha_values[row_index, col_index] * partial_residuals_x
			return(partial_residuals)
		},
		
		get_partial_residuals_x = function(row_index, col_index) {
			#The numbers to go on the x axis of a partial residuals plot.
			self$design_matrix[self$focal_vector == row_index, col_index]
		}
	)
)
