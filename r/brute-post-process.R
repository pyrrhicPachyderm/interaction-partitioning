#A set of helper functions for post-processing the data from brute forcing.
#Works in an object-oriented fashion.
#Constructing the central object requires a data file produced by the C++ program, with column names preserved.
#Additionally, it requires a vector of species names.
#These must be provided in the order that the C++ program indexed the species.

information_criterion_weights <- function(information_criterion) {
	#Calculates Akaike weights from AIC values (or Schwarz weights from BIC values, etc.).
	information_criterion <- information_criterion - min(information_criterion)
	weights <- exp(-0.5 * information_criterion)
	weights <- weights / sum(weights)
	return(weights)
}

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
	public = list(
		num_species = NULL,
		species_names = NULL,
		row_groupings = NULL, #A data frame of row groupings, named by the species names.
		col_groupings = NULL, #As above, for col groupings.
		residuals = NULL, #A data frame of residuals.
		statistics = NULL, #A data frame of all output values besides the above.
		
		initialize = function(data_file_name, species_names) {
			data_table <- read.table(data_file_name, header=TRUE)
			
			self$num_species <- length(species_names)
			self$species_names <- species_names
			
			private$validate_data(data_table, self$num_species)
			
			#Extract the groupings and other data.
			#Add one to the groupings to make them 1-indexed, as R prefers.
			self$row_groupings <- data_table[,grep("row_group", names(data_table))] + 1
			self$col_groupings <- data_table[,grep("col_group", names(data_table))] + 1
			names(self$row_groupings) <- species_names
			names(self$col_groupings) <- species_names
			self$residuals <- data_table[,grep("residual", names(data_table))]
			self$statistics <- data_table[,grep("row_group|col_group|residual", names(data_table), invert=TRUE)]
			
			#Add additional columns of statistics.
			if(!is.null(self$statistics$aic)) private$add_aic_weights()
		}
	),
	
	private = list(
		validate_data = function(data_table, num_species) {
			#Checks whether a data table is for the correct number of species.
			if(
				length(grep("row_group", names(data_table))) != num_species |
				length(grep("col_group", names(data_table))) != num_species
			) {
				stop("Data table has the wrong number of species.")
			}
		},
		
		add_aic_weights = function() {
			self$statistics$aic_weight <- information_criterion_weights(self$statistics$aic)
		},
		
		match_grouping = function(desired_grouping, grouping_table) {
			#Returns the indices at which a row of grouping_table is equal to the vector desired_grouping.
			
			elementwise_match <- lapply(1:ncol(grouping_table), function(i) {
				return(grouping_table[[i]] == desired_grouping[i])
			})
			
			rowwise_match <- Reduce("&", elementwise_match)
			
			return(which(rowwise_match))
		},
		
		get_statistics_row = function(row_grouping, col_grouping) {
			#row_grouping and col_grouping should be vectors representing the groupings.
			#Extracts one row of statistics, matching those groupings.
			row_matches <- private$match_grouping(row_grouping, self$row_groupings)
			col_matches <- private$match_grouping(col_grouping, self$col_groupings)
			return(self$statistics[intersect(row_matches, col_matches),])
		},
		
		match_groupings = function(grouping_table_1, grouping_table_2) {
			#Returns the indices at which the rows of the two grouping tables are equivalent.
			rowwise_match <- sapply(1:nrow(grouping_table_1), function(i) {
				return(all(as.vector(as.matrix(grouping_table_1[i,])) == as.vector(as.matrix(grouping_table_2[i,]))))
			})
			return(rowwise_match)
		},
		
		annotate_matrix = function(mat) {
			#Annotates a num_species by num_species matrix with the species names.
			rownames(mat) <- self$species_names
			colnames(mat) <- self$species_names
			return(mat)
		},
		
		get_coclassification_matrix = function(grouping_table, index) {
			#Returns a binary matrix of whether two species are grouped together in a given row a grouping table.
			grouping <- as.vector(as.matrix(grouping_table[index,]))
			mat <- as.matrix(sapply(grouping, function(group) {
				return(grouping == group)
			}))
			return(mat)
		},
		
		get_weighted_coclassification_matrix = function(grouping_table) {
			#Returns a coclassification matrix, as above, weighted by AIC weights.
			weights <- self$statistics$aic_weight
			weighted_matrices <- lapply(1:nrow(grouping_table), function(i) {
				return(private$get_coclassification_matrix(grouping_table, i) * weights[i])
			})
			return(Reduce("+", weighted_matrices))
		}
	),
	
	active = list(
		equivalently_grouped_statistics = function() {
			#Where the row grouping is the same as the column grouping.
			return(self$statistics[private$match_groupings(self$row_groupings,self$col_groupings),])
		},
		fully_grouped_statistics = function() {
			grouping <- rep(1, self$num_species)
			return(private$get_statistics_row(grouping, grouping))
		},
		fully_separated_statistics = function() {
			grouping <- 1:self$num_species
			return(private$get_statistics_row(grouping, grouping))
		},
		
		weighted_row_coclassification_matrix = function() {
			return(private$annotate_matrix(private$get_weighted_coclassification_matrix(self$row_groupings)))
		},
		weighted_col_coclassification_matrix = function() {
			return(private$annotate_matrix(private$get_weighted_coclassification_matrix(self$col_groupings)))
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
