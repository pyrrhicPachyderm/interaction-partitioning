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
#Second, a wrapper around is_consecutive_blocks that takes the base ordering and a grouping.
forms_consecutive_blocks <- function(grouping, ordering) {
	is_consecutive_blocks(grouping[ordering])
}
#Third, another helper function to get a list of all orderings.
#In C++, to make use of std::next_permutation().
Rcpp::cppFunction("List get_all_orderings(int n) {
	List all_orderings = List::create();
	std::vector<int> ordering(n);
	for(int i = 0; i < n; i++) ordering[i] = i+1;
	do {
		all_orderings.push_back(NumericVector::import(ordering.begin(), ordering.end()));
	} while(std::next_permutation(ordering.begin(), ordering.end()));
	return all_orderings;
}")
#Now, the function to exhaustively search and find an appropriate ordering.
#The function takes a list of groupings, and returns an ordering.
#It is hence sufficiently general to be constrained by any number of groupings, not only two.
get_ordering <- function(groupings) {
	all_orderings <- get_all_orderings(length(groupings[[1]]))
	for(ordering in all_orderings) {
		if(all(sapply(groupings, forms_consecutive_blocks, ordering))) {
			return(ordering)
		}
	}
	stop("Cannot find an ordering of species that makes all groups consecutive")
}
