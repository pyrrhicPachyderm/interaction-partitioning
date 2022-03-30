#A function for printing a grouped alpha matrix as a table using knitr.
#Expects the species names, the row grouping, the column grouping, and the alpha matrix.
#The groupings are given as a vector, [number of species] in length, with the number of the group each species belongs to.
#The matrix is defined in terms of the groups, and hence is [number of row groups] by [number of column groups].

#First, we need to determine an ordering of the species such that all species that are grouped together are adjacent in all groupings.
#In the absence of an elegant method, I'll do this by an exhaustive search of orderings.
#An ordering is given as a vector of integers: a permutation of the original order.
#Using R's subsetting, species can be reordered according to an ordering using species_names[ordering].
#First, a helper function, to check if a grouping is in consecutive blocks.
is_consecutive_blocks <- function(grouping) {
	current_group <- grouping[1]
	completed_groups <- c()
	for(g in grouping[-1]) {
		if(g != current_group) {
			if(g %in% completed_groups) return(FALSE)
			completed_groups <- c(completed_groups, current_group)
			current_group <- g
		}
	}
	return(TRUE)
}
