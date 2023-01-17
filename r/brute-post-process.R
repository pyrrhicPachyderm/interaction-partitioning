#A set of helper functions for post-processing the data from brute forcing.
#Depends upon post-process.R

library(magrittr)

prob_perm_test <- function(selected_probs, all_probs, num_iter = 1000) {
	#Tests whether the sum of selected_probs is more or less (as appropriate) than the sum of a random sample of all_probs of the same length.
	#Returns the p value.
	#TODO: Check that this is a valid statistical test.
	selected_sum <- sum(selected_probs)
	resamples <- replicate(num_iter-1,
		sum(sample(all_probs, length(selected_probs)))
	)
	resamples <- c(selected_sum, resamples)
	
	#Try the test both ways around, return whichever is appropriate.
	p_less <- mean(resamples <= selected_sum)
	p_more <- mean(resamples >= selected_sum)
	return(min(p_less, p_more))
}

BruteData <- R6::R6Class("BruteData",
	inherit = Data,
	
	public = list(
		weight_column = NULL, #The name of the column to use for weights.
		
		initialize = function(data_file_name, species_names, weight_column = c("aic_weight", "aicc_weight")) {
			super$initialize(data_file_name, species_names)
			self$weight_column <- match.arg(weight_column)
		}
	),
	
	active = list(
		average_parameters = function() {
			Parameters$new(lapply(
				self$parameters,
				weighted.mean,
				w = self$statistics[[self$weight_column]]
			))
		},
		
		unweighted_average_coclassification_matrix = function() {
			return(private$annotate_matrix(private$get_average_coclassification_matrix(self$row_groupings)))
		},
		weighted_row_coclassification_matrix = function() {
			return(private$annotate_matrix(private$get_weighted_coclassification_matrix(self$row_groupings, self$statistics[[self$weight_column]])))
		},
		weighted_col_coclassification_matrix = function() {
			return(private$annotate_matrix(private$get_weighted_coclassification_matrix(self$col_groupings, self$statistics[[self$weight_column]])))
		},
		
		equivalently_grouped_weight = function() {
			return(sum(self$equivalently_grouped_statistics[[self$weight_column]]))
		},
		equivalently_grouped_expected_weight = function() {
			return(nrow(self$equivalently_grouped_statistics) / nrow(self$statistics))
		},
		equivalently_grouped_weight_p = function() {
			#p value of a permutation test
			return(prob_perm_test(self$equivalently_grouped_statistics[[self$weight_column]], self$statistics[[self$weight_column]]))
		}
	)
)
