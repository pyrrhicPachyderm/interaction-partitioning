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

#Anther helper function, for use in zero padding numbers.
#Returns the number of digits in the integer part of a number.
get_integer_part_width <- function(x) {
	max(1, 1 + log10(abs(x)))
}

#A helper function that converts numbers to letters, as numbers can't be used in TeX macro names.
#Takes a non-negative number; some adding of 1 is necessary to convert to 1-indexing.
#0: a, ..., 25: z, 26: aa, 27: ab, ..., 701: zz, 702: aaa, etc.
number_to_letters <- function(n) {
	num_letters <- 1
	while(length(letters)^num_letters <= n) {
		n <- n - length(letters)^num_letters
		num_letters <- num_letters + 1
	}
	result <- ""
	for(i in 1:num_letters) {
		n_remainder <- n %% (length(letters)^i)
		letter_index <- floor(n_remainder / (length(letters)^(i-1))) #Zero-indexed.
		result <- paste0(letters[letter_index + 1], result) #Plus 1 to convert to 1-indexed.
	}
	return(result)
}

#Finally, the function to actually print the grouped matrix.
#We want what is essentially a static variable, to provide unique names to TeX macro.
#In particular, the macro that will serve as a length to store the column width.
#Per https://stackoverflow.com/a/59126599, we use a local environment.
grouped_alpha_matrix <- local({call_num <- 0; function(species_names, row_grouping, col_grouping, mat) {
	num_species <- length(species_names)
	
	#Reorder the species to place grouped species adjacent.
	ordering <- get_ordering(list(row_grouping, col_grouping))
	species_names <- species_names[ordering]
	row_grouping <- row_grouping[ordering]
	col_grouping <- col_grouping[ordering]
	
	#Configuration options.
	decimal_places <- 3
	use_dashed_lines <- TRUE
	dashed_line_spec <- "1pt/1pt" #Dash/gap.
	header_angle <- 30
	
	#Get a unique alphabetic id for this call to the function, using the "static variable" call_num.
	#Use it to define a length macro.
	call_id <- number_to_letters(call_num)
	length_macro <- sprintf("\\colwidth%s", call_id)
	
	#Define some helper functions/variables.
	alignment <- function(num_cols) {
		sprintf(">{\\centering\\arraybackslash}p{\\dimexpr %d %s\\relax}", num_cols, length_macro)
	}
	
	#Determine the width of the widest cell entry.
	#Following the decimal point will always be decimal_places characters.
	#So we only need to know how many digits precede the decimal places.
	#And whether there're any negative signs.
	integer_part_pad_width <- max(get_integer_part_width(mat))
	do_negative_pad <- any(mat < 0)
	
	#Determine the width of the widest column as a TeX length.
	widest_number <- paste0(c(rep("0", integer_part_pad_width), ".", rep("0", decimal_places)), collapse="")
	if(do_negative_pad) widest_number <- paste0("{-}", widest_number)
	length_definition <- sprintf("\\newlength{%s}\n\\settowidth{%s}{$%s$}\n", length_macro, length_macro, widest_number)
	
	#The beginning and end of the table.
	begin <- paste(c("\\begin{tabular}{", rep(alignment(1), num_species), "l}"), collapse="")
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
	#The use of multirow we can't use dcolumn or siunitx to align on decimal points, so we have to do it
	#by fixing the number of decimal places, and using phantoms to pad to a fixed width.
	get_group_content <- function(row_index, col_index) {
		number <- mat[row_index, col_index]
		#If there is a minus sign, it should be unary, not binary.
		#This is also necessary for alignment, as phantom minus signs will be treated as unary.
		#So the minus sign must be wrapped in braces, which means printing it separately.
		number_string <- sprintf(paste0("%.",decimal_places,"f"), abs(number)) #Format with the correct number of decimal places.
		if(number < 0) {
			number_string <- paste0("{-}", number_string) #Prepend a unary minus if necessary.
		}
		hphantom <- paste(rep("0", integer_part_pad_width - get_integer_part_width(number)), collapse="") #Pad with phantom zeroes.
		if(do_negative_pad && number >= 0) {
			hphantom <- paste0(hphantom, "-") #Pad with an phantom negative sign if necessary.
		}
		result <- sprintf(paste0("$\\hphantom{%s}%s$"), hphantom, number_string)
		return(result)
	}
	
	#A helper function to get the content of a particular group.
	#Wrapped in \multicolumn and \multirow to make it the right size for that group.
	get_expanded_group_content <- function(row_index, col_index) {
		num_rows <- sum(row_grouping == row_grouping[row_index])
		num_cols <- sum(col_grouping == col_grouping[col_index])
		#TODO: Add vertical rules.
		sprintf("\\multicolumn{%d}{%s}{\\multirow{%d}*{%s}}", num_cols, alignment(num_cols), num_rows, get_group_content(row_index, col_index))
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
	
	cat(length_definition)
	cat(begin)
	cat(label_line)
	cat(row_lines)
	cat(end)
}})
