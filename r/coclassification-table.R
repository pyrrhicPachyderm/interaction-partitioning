#A function for printing a coclassification matrix as a coloured table using knitr.
#Expects a coclassification matrix (with numeric entries in [0,1]), with row and column names.
#Uses a pgfplots colour map and xcolor colour name for the diagonal.

library(magrittr)

strip_booktabs_rules <- function(tab) {
	#Strips out all the rules and such added by booktabs.
	#This produces a table without rules.
	#This actually removes dependency on booktabs.
	#But it's easier to remove booktab's rules than *all* the standard rules.
	tab %>%
		sub("\\\\toprule", "", .) %>%
		sub("\\\\midrule", "", .) %>%
		sub("\\\\addlinespace", "", .) %>%
		sub("\\\\bottomrule", "", .)
}

weighted_coclassification_kable <- function(mat, colourmap="hot", diagonal_colour="lightgray", digits=3) {
	#colourmap is the name of a pgfplots colourmap.
	
	array_stretch <- 1.5
	header_angle <- 30
	
	#\dimexpr can't work with floating point multipliers, so we need rationals.
	as_rational_string <- function(x) {
		as.character(gmp::as.bigq(x))
	}
	
	#\baselineskip seemingly doesn't work inside a tikzpicture, so we save it elsewhere.
	length_definition <- paste(
		"\\ifdefined\\coclassificationtextheight\\else\\newlength{\\coclassificationtextheight}\\fi",
		"\\setlength{\\coclassificationtextheight}{\\baselineskip}",
	sep="\n")
	
	#The height of the whole table.
	table_height <- paste0("{\\dimexpr\\coclassificationtextheight * ",as_rational_string(array_stretch)," * ",nrow(mat),"}")
	#The difference between the baseline of the bottom row and the bottom of the table.
	baseline_drop <- paste0("{\\dimexpr -\\coclassificationtextheight * ",as_rational_string(array_stretch-1)," / 2}")
	
	colourbar_code <- paste(
		"\\begin{tikzpicture}[baseline]",
			"\\begin{axis}[",
				paste0("colormap name=",colourmap,","),
				"colorbar left,",
				"point meta min=0,",
				"point meta max=1,",
				"colorbar style={",
					paste0("height=",table_height,","),
					"at={(parent axis.south east)},",
					"anchor=south west,",
					paste0("yshift=",baseline_drop,","),
				"},",
				"%Disable the axes themselves; we only want the colourbar.",
				"hide axis, scale only axis, width=0pt, height=0pt,",
			"]\\end{axis}",
		"\\end{tikzpicture}",
	sep="\n")
	
	format_num <- function(num) {
		sprintf(paste0("%.",digits,"f"), num) %>%
			sub("0.", ".", .) #Strip the leading 0.
	}
	
	#LaTeX doesn't like us using \pgfplotscolormapaccess mid-table.
	#So we will use it beforehand, define colours, then acces them.
	#This requires us to generate a unique colour name for each cell.
	#We will make it verbose to reduce the risk of a clash.
	get_colour_name <- function(row, col) {
		paste("weighted coclassification table colour", row, col)
	}
	
	#Now, the colour definitions.
	colour_definitions <- sapply(1:ncol(mat), function(col) {
		sapply(1:nrow(mat), function(row) {
			paste0(
				"\\pgfplotscolormapaccess[0:1]{",mat[row,col],"}{",colourmap,"}\n",
				"\\definecolor{",get_colour_name(row,col),"}{rgb}{\\pgfmathresult}"
			)
		})
	})
	colour_definitions <- paste(as.vector(colour_definitions), collapse="\n")
	
	#To include the \cellcolor commands, the entries of the table must be strings.
	#So we convert to a string matrix, inserting cell colourings.
	string_mat <- sapply(1:ncol(mat), function(col) {
		sapply(1:nrow(mat), function(row) {
			if(row == col) { #No numbers on the diagonals.
				return(paste0("\\cellcolor{",diagonal_colour,"}"))
			}
			value <- format_num(mat[row,col])
			paste0("\\cellcolor{",get_colour_name(row,col),"}",value)
		})
	})
	
	#\rotatebox makes things wider, so we enclose in a fixed width \makebox.
	colnames(string_mat) <- paste0("\\makebox[1em][l]{\\rotatebox{",header_angle,"}{\\emph{",colnames(mat),"}}}")
	#We don't set row names for mat, because we want the species labels on the right instead of the left.
	#We cbind instead.
	string_mat <- cbind(string_mat, matrix(
		paste0("\\emph{",rownames(mat),"}"),
	ncol=1))
	
	matrix_tab <- knitr::kable(string_mat, escape=FALSE, booktabs=TRUE,
		align=c(rep("c",ncol(mat)),"l")
	) %>%
		strip_booktabs_rules() %>%
		sub("\\begin{tabular}", "\\begin{tabular}[b]", ., fixed=TRUE) #Baseline at the bottom, for colourbar alignment.
	
	#Need \bgroup \egroup to contain the redefinition of \arraystretch.
	matrix_tab <- paste(
		paste0("\\bgroup\n\\renewcommand{\\arraystretch}{",array_stretch,"}"),
		matrix_tab,
		"\\egroup",
	sep="\n")
	
	outer_tab <- paste(
		"\\begin{tabular}{rl}",
		colourbar_code,
		"&",
		matrix_tab,
		"\\end{tabular}",
	sep="\n")
	
	cat(length_definition)
	cat(colour_definitions)
	cat(outer_tab)
}