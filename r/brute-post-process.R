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

BruteData <- R6::R6Class("BruteData",
	public = list(
		num_species = NULL,
		species_names = NULL,
		row_groupings = NULL, #A data frame of row groupings, named by the species names.
		col_groupings = NULL, #As above, for col groupings.
		statistics = NULL, #A data frame of all output values besides groupings.
		
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
			self$statistics <- data_table[,grep("row_group|col_group", names(data_table), invert=TRUE)]
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
		}
	),
	
	active = list(
		fully_grouped_statistics = function() {
			grouping <- rep(1, self$num_species)
			return(private$get_statistics_row(grouping, grouping))
		},
		fully_separated_statistics = function() {
			grouping <- 1:self$num_species
			return(private$get_statistics_row(grouping, grouping))
		}
	)
)
