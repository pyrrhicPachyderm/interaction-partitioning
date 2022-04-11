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

Data <- R6::R6Class("Data",
	public = list(
		num_species = NULL,
		species_names = NULL,
		row_groupings = NULL, #A data frame of row groupings, named by the species names.
		col_groupings = NULL, #As above, for col groupings.
		parameters = NULL, #A data frame of parameter values.
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
			self$parameters <- data_table[,grep("parameters", names(data_table))]
			private$objectify_parameters()
			self$statistics <- data_table[,grep("_[0-9]*$|parameters", names(data_table), invert=TRUE)]
			
			#Add additional columns of statistics.
			if(!is.null(self$statistics$aic)) private$add_aic_weights()
			if(!is.null(self$statistics$aicc)) private$add_aicc_weights()
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
		
		objectify_parameters = function() {
			#Turns self$parameters from a data frame of raw parameter values to a vector of R6 Parameters objects.
			
			growth_rates_df <- self$parameters[,grep("parameters_r_", names(self$parameters))]
			alpha_values_df <- self$parameters[,grep("parameters_alpha_", names(self$parameters))]
			
			self$parameters <- sapply(1:nrow(self$parameters), function(i) {
				growth_rates <- as.vector(as.matrix(growth_rates_df[i,]))
				row_grouping <- as.vector(as.matrix(self$row_groupings[i,]))
				col_grouping <- as.vector(as.matrix(self$col_groupings[i,]))
				group_alpha_values <- matrix(
					as.vector(as.matrix(alpha_values_df[i,])),
					nrow = self$num_species,
					byrow = TRUE
				)
				alpha_values <- group_alpha_values[row_grouping,col_grouping]
				
				return(Parameters$new(growth_rates, alpha_values))
			})
		},
		
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
			#Returns a binary matrix of whether two species are grouped together in a given row of a grouping table.
			grouping <- as.vector(as.matrix(grouping_table[index,]))
			mat <- as.matrix(sapply(grouping, function(group) {
				return(grouping == group)
			}))
			return(mat)
		},
		
		get_weighted_coclassification_matrix = function(grouping_table, weights) {
			#Returns a coclassification matrix, as above, averaged across the whole table.
			#Uses a weighted average by a given set of weights.
			weighted_matrices <- lapply(1:nrow(grouping_table), function(i) {
				return(private$get_coclassification_matrix(grouping_table, i) * weights[i])
			})
			return(Reduce("+", weighted_matrices))
		},
		
		get_average_coclassification_matrix = function(grouping_table) {
			#As above, but with all weights equal.
			weights <- rep(1 / nrow(grouping_table), nrow(grouping_table))
			return(private$get_weighted_coclassification_matrix(grouping_table, weights))
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
		
		fully_grouped_coclassification_matrix = function() {
			return(matrix(1, self$num_species, self$num_species))
		},
		fully_separated_coclassification_matrix = function() {
			return(diag(self$num_species))
		}
	)
)
