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
#Next, we need another helper function to get a list of all orderings.
#Ordinarily, I would do this in C++, to make use of std::next_permutation().
#However, you can't conveniently save an Rcpp function to an rda file.
#So, third, I implement next_permutation in R.
next_permutation <- function(perm) {
	#Using the linear-time implementation described here: https://wordaligned.org/articles/next-permutation
	#First, check that the vector is long enough to have permutations.
	if(length(perm) < 2) return(perm)
	#Second, find the longest monotonically decreasing tail.
	head_of_tail <- NA
	for(i in (length(perm)-1):1) {
		if(perm[i+1] > perm[i]) {
			head_of_tail <- i+1
			break
		}
	}
	#If head_of_tail is NA, the entire permutation is monotonically decreasing, just reverse it.
	if(is.na(head_of_tail)) {
		return(perm[length(perm):1])
	}
	#Else, swap the end of the head with the the last element of the tail that is greater than it.
	temp <- perm[head_of_tail-1]
	for(i in length(perm):head_of_tail) {
		if(perm[i] > temp) {
			perm[head_of_tail-1] <- perm[i]
			perm[i] <- temp
			break
		}
	}
	#Reverse the tail.
	perm[head_of_tail:length(perm)] <- perm[length(perm):head_of_tail]
	return(perm)
}
#Fourth, another helper function to get a list of all orderings.
get_all_orderings <- function(n) {
	all_orderings <- vector("list", factorial(n))
	all_orderings[[1]] <- 1:n
	for(i in 2:length(all_orderings)) {
		all_orderings[[i]] <- next_permutation(all_orderings[[i-1]])
	}
	return(all_orderings)
}
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

#Finally, the function to actually print the grouped matrix.
grouped_alpha_matrix <- function(species_names, row_grouping, col_grouping, mat, is_ordered=TRUE, array_stretch=1.5) {
	num_species <- length(species_names)
	
	if(is_ordered) {
		#Assert that the ordering is correct.
		if(!forms_consecutive_blocks(row_grouping, 1:num_species) || !forms_consecutive_blocks(col_grouping, 1:num_species)) {
			stop("Given ordering does not make all groups consecutive")
		}
	} else {
		#Reorder the species to place grouped species adjacent.
		ordering <- get_ordering(list(row_grouping, col_grouping))
		species_names <- species_names[ordering]
		row_grouping <- row_grouping[ordering]
		col_grouping <- col_grouping[ordering]
		mat <- mat[ordering,ordering]
	}
	
	#Configuration options.
	decimal_places <- 3
	use_dashed_lines_inner <- TRUE
	use_dashed_lines_outer <- FALSE
	dashed_line_spec <- "3pt/3pt" #Dash/gap.
	header_angle <- 30
	
	#Define some helper functions/variables.
	alignment <- function(num_cols) {
		sprintf(">{\\centering\\arraybackslash}p{\\dimexpr %d %s\\relax}", num_cols, length_macro)
	}
	length_macro <- "\\colwidth"
	divider <- function(use_dashed_lines) {
		ifelse(use_dashed_lines, sprintf(";{%s}",dashed_line_spec), "|")
	}
	cline <- function(use_dashed_lines) {
		ifelse(use_dashed_lines,
			function(col1, col2){sprintf("\\cdashline{%d-%d}[%s]", col1, col2, dashed_line_spec)},
			function(col1, col2){sprintf("\\cline{%d-%d}", col1, col2)}
		)
	}
	outer_divider <- divider(use_dashed_lines_outer)
	inner_divider <- divider(use_dashed_lines_inner)
	outer_cline <- cline(use_dashed_lines_outer)
	inner_cline <- cline(use_dashed_lines_inner)
	
	#Determine the width of the widest cell entry.
	#Following the decimal point will always be decimal_places characters.
	#So we only need to know how many digits precede the decimal places.
	#And whether there're any negative signs.
	integer_part_pad_width <- max(get_integer_part_width(mat))
	do_negative_pad <- any(mat < 0)
	
	#Determine the width of the widest column as a TeX length.
	widest_number <- paste0(c(rep("0", integer_part_pad_width), ".", rep("0", decimal_places)), collapse="")
	if(do_negative_pad) widest_number <- paste0("{-}", widest_number)
	length_definition <- sprintf("\\ifdefined%s\\relax\\else\\newlength{%s}\\fi\n\\settowidth{%s}{$%s$}\n", length_macro, length_macro, length_macro, widest_number)
	
	#The beginning and end of the table.
	begin <- paste(c(
		"\\begingroup",
		"\\renewcommand*\\arraystretch{", array_stretch, "}",
		"\\begin{tabular}{", rep(alignment(1), num_species), "l}"
	), collapse="")
	end <- "\\end{tabular}\\endgroup"
	
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
	#Also adds vertical rules.
	#Adding the vertical dividers is interfered with by multirow, so we need some empty \multicolumns.
	#is_empty_content produces such.
	get_expanded_group_content <- function(row_index, col_index, is_first_group, is_last_group, is_empty_content) {
		num_rows <- sum(row_grouping == row_grouping[row_index])
		num_cols <- sum(col_grouping == col_grouping[col_index])
		pre_divider <- ifelse(is_first_group, outer_divider, inner_divider)
		post_divider <- ifelse(is_last_group, outer_divider, "") #No post divider if next column has a pre divider.
		content <- ifelse(is_empty_content, "",
			sprintf("\\multirow{%d}*{%s}", num_rows, get_group_content(row_index, col_index))
		)
		sprintf("\\multicolumn{%d}{%s%s%s}{%s}", num_cols, pre_divider, alignment(num_cols), post_divider, content)
	}
	
	#A function to get the contents of any cell of the table.
	#Including cells that have no content, as their content is handled elsewhere using \multicolumn or \multirow.
	#When using \multicolumn, we must skip alignment tab characters, `&`.
	#As such, each cell will include the `&` from the end of the cell.
	#The species names are going on the end, so they simply shalln't have an `&`.
	get_cell <- function(row_index, col_index) {
		#First, check if this is the first column of the current group.
		if(col_index == 1 || col_grouping[col_index] != col_grouping[col_index-1]) {
			#If this is the first column, we need a multicolumn.
			#We also check if this is the first row of the current group.
			#If so, the multicolumn has actual content.
			is_first_group <- col_index == 1
			is_last_group <- col_grouping[col_index] == col_grouping[num_species]
			has_content <- row_index == 1 || row_grouping[row_index] != row_grouping[row_index-1]
			return(sprintf("%s&", get_expanded_group_content(row_index, col_index, is_first_group, is_last_group, !has_content)))
		} else {
			#If this is not the first column, we leave it blank.
			#We don't even put `&` in; \multicolumn handles this.
			return("")
		}
	}
	
	get_row <- function(row_index) {
		#First, check if this is the top row of the current group.
		if(row_index == 1 || row_grouping[row_index] != row_grouping[row_index-1]) {
			pre_cline <- ifelse(row_index == 1, outer_cline(1, num_species), inner_cline(1, num_species))
		} else {
			pre_cline <- ""
		}
		#A post cline is only necessary for the very last row.
		post_cline <- ifelse(row_index == num_species, outer_cline(1, num_species), "")
		
		cells <- sapply(1:num_species, function(col_index){get_cell(row_index,col_index)})
		paste(c(pre_cline, cells, get_species_row_label(row_index), "\\\\", post_cline), collapse="")
	}
	
	row_lines <- paste(sapply(1:num_species, get_row), collapse="")
	
	return(paste0(
		length_definition,
		begin,
		label_line,
		row_lines,
		end
	))
}
