#A set of helper functions for post-processing the data from brute forcing.
#Depends upon post-process.R

library(magrittr)

prob_perm_test <- function(selected_probs, all_probs, num_iter = 1000) {
	#Tests whether the sum of selected_probs is less than the sum of a random sample of all_probs of the same length.
	#Returns the p value.
	#TODO: Check that this is a valid statistical test.
	selected_sum <- sum(selected_probs)
	resamples <- replicate(num_iter-1,
		sum(sample(all_probs, length(selected_probs)))
	)
	resamples <- c(selected_sum, resamples)
	p <- mean(resamples <= selected_sum)
	return(p)
}

BruteData <- R6::R6Class("BruteData",
	inherit = Data,
	
	active = list(
		weighted_row_coclassification_matrix = function() {
			return(private$annotate_matrix(private$get_weighted_coclassification_matrix(self$row_groupings, self$statistics$aic_weight)))
		},
		weighted_col_coclassification_matrix = function() {
			return(private$annotate_matrix(private$get_weighted_coclassification_matrix(self$col_groupings, self$statistics$aic_weight)))
		},
		
		equivalently_grouped_weight = function() {
			return(sum(self$equivalently_grouped_statistics$aic_weight))
		},
		equivalently_grouped_expected_weight = function() {
			return(nrow(self$equivalently_grouped_statistics) / nrow(self$statistics))
		},
		equivalently_grouped_weight_p = function() {
			#p value of a permutation test
			return(prob_perm_test(self$equivalently_grouped_statistics$aic_weight, self$statistics$aic_weight))
		}
	)
)
