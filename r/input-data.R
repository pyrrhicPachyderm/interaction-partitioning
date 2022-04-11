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
		}
	)
)
