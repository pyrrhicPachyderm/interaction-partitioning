#A set of helper functions for post-processing the data from the C++ programs.
#Works in an object-oriented fashion.
#Constructing the central object requires a data file produced by the C++ program, with column names preserved.
#Additionally, it requires a vector of species names.
#These must be provided in the order that the C++ program indexed the species.

library(magrittr)

information_criterion_weights <- function(information_criterion) {
	#Calculates Akaike weights from AIC values (or Schwarz weights from BIC values, etc.).
	information_criterion <- information_criterion - min(information_criterion)
	weights <- exp(-0.5 * information_criterion)
	weights <- weights / sum(weights)
	return(weights)
}

as_coclassification_matrix <- function(grouping) {
	#Converts a grouping to a binary coclassification matrix.
	mat <- sapply(grouping, function(group) {
		return(grouping == group)
	})
	return(mat)
}

Data <- R6::R6Class("Data",
	public = list(
		num_row_species = NULL,
		num_col_species = NULL,
		species_names = NULL,
		row_groupings = NULL, #A data frame of row groupings, named by the species names.
		col_groupings = NULL, #As above, for col groupings.
		parameters = NULL, #A data frame of parameter values.
		statistics = NULL, #A data frame of all output values besides the above.
		
		initialize = function(data_file_name, species_names) {
			data_table <- read.table(data_file_name, header=TRUE)
			
			self$species_names <- species_names
			
			#Extract the groupings and other data.
			#Add one to the groupings to make them 1-indexed, as R prefers.
			self$row_groupings <- data_table[,grep("row_group", names(data_table))] + 1
			self$col_groupings <- data_table[,grep("col_group", names(data_table))] + 1
			self$num_row_species <- ncol(self$row_groupings)
			self$num_col_species <- ncol(self$col_groupings)
			names(self$row_groupings) <- species_names[1:self$num_row_species]
			names(self$col_groupings) <- species_names[1:self$num_col_species]
			self$parameters <- data_table[,grep("parameters", names(data_table))]
			self$statistics <- data_table[,grep("_[0-9]*$|parameters", names(data_table), invert=TRUE)]
			
			#Add additional columns of statistics.
			if(!is.null(self$statistics$aic)) private$add_aic_weights()
			if(!is.null(self$statistics$aicc)) private$add_aicc_weights()
		},
		
		get_row_grouping = function(index) {
			as.vector(as.matrix(self$row_groupings[index,]))
		},
		get_col_grouping = function(index) {
			as.vector(as.matrix(self$col_groupings[index,]))
		},
		get_parameters = function(index) {
			Parameters$new(self$parameters[index,])
		}
	),
	
	private = list(
		add_aic_weights = function() {
			self$statistics$aic_weight <- information_criterion_weights(self$statistics$aic)
		},
		
		add_aicc_weights = function() {
			self$statistics$aicc_weight <- information_criterion_weights(self$statistics$aicc)
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
			#Extracts those rows of statistics matching those groupings.
			row_matches <- private$match_grouping(row_grouping, self$row_groupings)
			col_matches <- private$match_grouping(col_grouping, self$col_groupings)
			return(self$statistics[intersect(row_matches, col_matches),])
		},
		
		match_groupings = function(grouping_table_1, grouping_table_2) {
			#Returns the indices at which the rows of the two grouping tables are equivalent.
			elementwise_match <- lapply(1:ncol(grouping_table_1), function(i) {
				return(grouping_table_1[[i]] == grouping_table_2[[i]])
			})
			
			rowwise_match <- Reduce("&", elementwise_match)
			
			return(which(rowwise_match))
		},
		
		annotate_matrix = function(mat) {
			#Annotates a matrix with the species names.
			rownames(mat) <- self$species_names[1:nrow(mat)]
			colnames(mat) <- self$species_names[1:ncol(mat)]
			return(mat)
		},
		
		get_weighted_coclassification_matrix = function(grouping_table, weights) {
			#Returns a coclassification matrix, averaged across the whole table.
			#Uses a weighted average by a given set of weights.
			sapply(grouping_table, function(group1) {
				sapply(grouping_table, function(group2) {
					sum((group1 == group2) * weights) / sum(weights)
				})
			})
		},
		
		get_average_coclassification_matrix = function(grouping_table) {
			#As above, but with all weights equal.
			weights <- rep(1, nrow(grouping_table))
			return(private$get_weighted_coclassification_matrix(grouping_table, weights))
		}
	),
	
	active = list(
		equivalently_grouped_statistics = function() {
			#Where the row grouping is the same as the column grouping.
			return(self$statistics[private$match_groupings(self$row_groupings,self$col_groupings),])
		},
		fully_grouped_statistics = function() {
			row_grouping <- rep(1, self$num_row_species)
			col_grouping <- rep(1, self$num_col_species)
			return(private$get_statistics_row(row_grouping, col_grouping))
		},
		fully_separated_statistics = function() {
			row_grouping <- 1:self$num_row_species
			col_grouping <- 1:self$num_col_species
			return(private$get_statistics_row(row_grouping, col_grouping))
		},
		
		fully_grouped_coclassification_matrix = function(num_species) {
			return(matrix(1, num_species, num_species))
		},
		fully_separated_coclassification_matrix = function(num_species) {
			return(diag(num_species))
		},
		
		fully_grouped_row_coclassification_matrix = function(num_species) {return(fully_grouped_coclassification_matrix(self$num_row_species))},
		fully_grouped_col_coclassification_matrix = function(num_species) {return(fully_grouped_coclassification_matrix(self$num_col_species))},
		fully_separated_row_coclassification_matrix = function(num_species) {return(fully_separated_coclassification_matrix(self$num_row_species))},
		fully_separated_col_coclassification_matrix = function(num_species) {return(fully_separated_coclassification_matrix(self$num_col_species))}
	)
)
