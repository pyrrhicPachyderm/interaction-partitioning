#A set of helper functions for post-processing the data from reversible jump MCMC.
#Depends upon post-process.R

library(magrittr)

get_chain_lengths <- function(chain_ids) {
	as.vector(table(chain_ids))
}

get_rhat <- function(parameter, chain_id) {
	#Calculates the Gelman-Rubin statistic for a single parameter.
	#Per Bayesian Data Analysis, Third Edition, Gelman et al., pgs. 284-285.
	
	#First, split into a list of chains.
	chains <- split(parameter, chain_id)
	#Split each chain in two.
	chains <- unlist(lapply(chains, function(v){
		list(v[1:(length(v)/2)], v[(length(v)/2+1):length(v)])
	}), recursive = FALSE)
	
	#Verify all chains are the same length (necessary for the calculation of n).
	stopifnot(length(unique(sapply(chains, length))) == 1)
	
	#Perform the calculations.
	n <- length(chains[[1]])
	b <- n * var(sapply(chains, mean))
	w <- mean(sapply(chains, var))
	v <- (n-1) / n * w + 1 / n * b
	rhat <- sqrt(v / w)
	
	return(rhat)
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
