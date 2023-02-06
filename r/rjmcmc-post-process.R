#A set of helper functions for post-processing the data from reversible jump MCMC.
#Depends upon post-process.R

library(magrittr)

get_chain_lengths <- function(chain_ids) {
	as.vector(table(chain_ids))
}

RJMCMCData <- R6::R6Class("RJMCMCData",
	inherit = Data,
	
	public = list(
		chain_lengths = NULL, #The length of each MCMC chain.
		
		initialize = function(data_file_name, species_names) {
			super$initialize(data_file_name, species_names)
			
			if(!is.null(self$chain_id)) self$chain_lengths <- get_chain_lengths(self$chain_id)
		}
	),
	
	active = list(
		weighted_growth_coclassification_matrix = function() {
			return(private$annotate_matrix(private$get_average_coclassification_matrix(self$growth_groupings)))
		},
		weighted_row_coclassification_matrix = function() {
			return(private$annotate_matrix(private$get_average_coclassification_matrix(self$row_groupings)))
		},
		weighted_col_coclassification_matrix = function() {
			return(private$annotate_matrix(private$get_average_coclassification_matrix(self$col_groupings)))
		}
	)
)
