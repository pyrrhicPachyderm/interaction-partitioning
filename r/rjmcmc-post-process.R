#A set of helper functions for post-processing the data from reversible jump MCMC.
#Depends upon post-process.R

library(magrittr)

RJMCMCData <- R6::R6Class("RJMCMCData",
	inherit = Data,
	
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
