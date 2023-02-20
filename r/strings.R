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
