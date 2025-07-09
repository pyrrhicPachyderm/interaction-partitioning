permanova_table <- function(dist_matrix, predictor, nperm = 4999) {
	return(vegan::adonis2(dist_matrix ~ predictor, as.data.frame(predictor), permutations = nperm))
}

permanova_f <- function(dist_matrix, predictor, nperm = 4999) {
	return(permanova_table(dist_matrix, predictor, nperm)[1,"F"])
}

permanova_df <- function(dist_matrix, predictor, nperm = 4999) {
	return(as.vector(permanova_table(dist_matrix, predictor, nperm)[1:2,"Df"]))
}

permanova_p <- function(dist_matrix, predictor, nperm = 4999) {
	return(permanova_table(dist_matrix, predictor, nperm)[1,"Pr(>F)"])
}
