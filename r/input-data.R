#An R6 object to read, store and handle the raw data that was fed into the C++ program.
#Requires the three data files that were provided to the C++ object:
#the focal vector, the response vector, and the design matrix.

InputData <- R6::R6Class("InputData",
	public = list(
		focal_vector = NULL,
		response_vector = NULL,
		design_matrix = NULL,
		is_per_capita = NULL,
		num_species = NULL,
		num_obs = NULL,
		
		initialize = function(focal_vector_file_name, response_vector_file_name, design_matrix_file_name, is_per_capita=FALSE) {
			self$focal_vector <- read.table(focal_vector_file_name)[[1]]
			self$response_vector <- read.table(response_vector_file_name)[[1]]
			self$design_matrix <- as.matrix(read.table(design_matrix_file_name))
			self$is_per_capita <- is_per_capita
			
			self$num_species <- ncol(self$design_matrix)
			self$num_obs <- nrow(self$design_matrix)
		},
		
		get_fitted_values = function(parameters) {
			focal_growth_rate <- parameters$growth_rates[self$focal_vector]
			focal_density <- mapply(
				function(obs,focal){self$design_matrix[obs,focal]},
				1:num_obs, self$focal_vector
			)
			alpha_values_on_focal <- parameters$alpha_values[self$focal_vector,]
			
			intrinsic_growth <- focal_growth_rate
			if(!self$is_per_capita) intrinsic_growth <- intrinsic_growth * focal_density
			total_competition <- rowSums(alpha_values_on_focal * self$design_matrix)
			fitted_values <- intrinsic_growth * (1 - total_competition)
			
			return(fitted_values)
		}
	)
)
