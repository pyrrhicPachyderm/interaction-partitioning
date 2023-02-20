#A set of helper functions for using R string in rnw files.

#Converts a string to a binomial name, with \emph{}.
#Works correctly on a vector of binomial names as well.
#If a string can't be parsed as a binomial name, it is returned unchanged.
#This is intentional; it allows "other" or the like to be included in lists of species names, with italicisation.
#Can handle an abbreviated genus; this is still treated as a valid binomial name.
#A genus name followed by "sp", "sp.", "spp", or "spp." has the genus name italicised, but not the remainder.
as_binomial_name <- function(names) {
	#NB: Backslashes need to be quadrupled: escaped once as they're written, and once more as they go through sub.
	#When matching a word (a sequence of one or more non-space characters), we also require that it doesn't contain backslashes.
	#Hence "[^\ ]". This ensures we don't match words with LaTeX commands.
	italics_command <- "\\\\emph"
	#First, handle the "sp", "sp.", "spp" or "spp." cases.
	names <- sub("^([^\\\\ ]+) (sp|sp\\.|spp|spp\\.)$", paste0(italics_command, "{\\1} \\2"), names)
	#Handle all other cases.
	names <- sub("^([^\\\\ ]+ [^\\\\ ]+)$", paste0(italics_command, "{\\1}"), names)
	return(names)
}

#Converts a vector of strings to a single string.
#Concatenates the strings, separating them with a comma and space,
#except the last two, which are separated with "and".
#If end_punct is set, it is appended to the end, unless the final character is already equal to end_punct.
#This is useful if the list might end with "sp.", say.
#Terminating "}" characters will be ignored when deciding whether to append end_punct.
as_english_list <- function(strings, end_punct = "", oxford_comma = TRUE) {
	if(length(strings) == 1) { #Handle the degenerate edge case.
		result <- strings
	} else {
		result <- paste0(
			paste(strings[1:(length(strings)-1)], collapse = ", "),
			ifelse(oxford_comma, ",", ""),
			" and ",
			strings[length(strings)]
		)
	}
	#We can't directly use grep for end_punct, as it'll probably be ".", which is a special character for grep.
	#Instead, we sub() any terminating "}" from the end, then cut a substring of the correct length, and check for equality.
	trimmed_result <- sub("}*$", "", result)
	if(substr(trimmed_result, nchar(trimmed_result) - nchar(end_punct) + 1, nchar(trimmed_result)) != end_punct) {
		result <- paste0(result, end_punct)
	}
	return(result)
}
