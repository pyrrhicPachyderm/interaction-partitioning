#A set of helper functions for post-processing the data from brute forcing.
#Works in an object-oriented fashion.
#Constructing the central object requires a data file produced by the C++ program, with column names preserved.
#Additionally, it requires a vector of species names.
#These must be provided in the order that the C++ program indexed the species.

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
			self$row_groupings <- data_table[,grep("row_group", names(data_table))]
			self$col_groupings <- data_table[,grep("col_group", names(data_table))]
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
		}
	)
)
