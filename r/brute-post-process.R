#A set of helper functions for post-processing the data from brute forcing.
#Does not read in any data; data must all be passed in by arguments.
#Data is assumed to be passed as produced by the C++ program, however, with column names preserved.
#Additionally, some of the functions rely on a vector of species names.
#These must be provided in the order that the C++ program indexed the species.

validate_data <- function(data, species_names) {
	#Ensures that this is data for the correct number of species.
	if(
		length(grep("row_group", names(data))) != length(species_names) |
		length(grep("col_group", names(data))) != length(species_names)
	) {
		stop("Data frame has the wrong number of species.")
	}
}
