#A function for printing a grouped alpha matrix as a table using knitr.
#Expects the species names, the row grouping, the column grouping, and the alpha matrix.
#The groupings are given as a vector, [number of species] in length, with the number of the group each species belongs to.
#The matrix is the full expanded matrix, not defined in terms of the groups.
#Hence the matrix is [number of species] by [number of species].

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

#Finally, the function to actually print the grouped matrix.
grouped_alpha_matrix <- function(species_names, row_grouping, col_grouping, mat) {
	num_species <- length(species_names)
	
	#Reorder the species to place grouped species adjacent.
	ordering <- get_ordering(list(row_grouping, col_grouping))
	species_names <- species_names[ordering]
	row_grouping <- row_grouping[ordering]
	col_grouping <- col_grouping[ordering]
	
	#Configuration options.
	use_dashed_lines <- TRUE
	dashed_line_spec <- "1pt/1pt" #Dash/gap.
	header_angle <- 30
	
	#Define some helper functions/variables.
	alignment <- "c"
	
	#The beginning and end of the table.
	begin <- paste(c("\\begin{tabular}{", rep(alignment, num_species), "l}"), collapse="")
	end <- "\\end{tabular}"
	
	#The row with the species names as column headers.
	get_species_column_label <- function(index) {
		#TODO: Sort out the width of dividers.
		#\rotatebox makes things wider, so we enclose it in a fixed width \makebox.
		sprintf("\\makebox[1em][l]{\\rotatebox{%d}{\\emph{%s}}}", header_angle, species_names[index])
	}
	species_column_labels <- sapply(1:num_species, get_species_column_label)
	label_line <- paste(c(species_column_labels,"\\\\"), collapse="&")
	
	get_species_row_label <- function(index) {
		sprintf("\\emph{%s}", species_names[index])
	}
	
	#A helper function to get the content of a particular group; simply the text to print in the table.
	get_group_content <- function(row_index, col_index) {
		#TODO: Format with the correct number of decimal places.
		sprintf("\\num{%f}", mat[row_index, col_index])
	}
	
	#A helper function to get the content of a particular group.
	#Wrapped in \multicolumn and \multirow to make it the right size for that group.
	get_expanded_group_content <- function(row_index, col_index) {
		num_rows <- sum(row_grouping == row_grouping[row_index])
		num_cols <- sum(col_grouping == col_grouping[col_index])
		#TODO: Add vertical rules and deal with the width properly.
		sprintf("\\multicolumn{%d}{%s}{\\multirow{%d}*{%s}}", num_cols, alignment, num_rows, get_group_content(row_index, col_index))
	}
	
	#A function to get the contents of any cell of the table.
	#Including cells that have no content, as their content is handled elsewhere using \multicolumn or \multirow.
	#When using \multicolumn, we must skip alignment tab characters, `&`.
	#As such, each cell will include the `&` from the end of the cell.
	#The species names are going on the end, so they simply shalln't have an `&`.
	get_cell <- function(row_index, col_index) {
		#First, check if this is the top row of the current group.
		if(row_index == 1 || row_grouping[row_index] != row_grouping[row_index-1]) {
			#If this is the first column of the current group, we define the content here.
			#Otherwise, we leave it blank.
			#If leaving it blank, we don't even put `&` in; \multicolumn handles this.
			if(col_index == 1 || col_grouping[col_index] != col_grouping[col_index-1]) {
				return(sprintf("%s&", get_expanded_group_content(row_index, col_index)))
			} else {
				return("")
			}
		} else {
			#If this is a subsequent row, the content has already been added above, so this is empty.
			#We just need `&` to tab past everything.
			return("&")
		}
	}
	
	get_row <- function(row_index) {
		cells <- sapply(1:num_species, function(col_index){get_cell(row_index,col_index)})
		paste(c(cells, get_species_row_label(row_index), "\\\\"), collapse="")
	}
	
	row_lines <- paste(sapply(1:num_species, get_row), collapse="")
	
	cat(begin)
	cat(label_line)
	cat(row_lines)
	cat(end)
}
