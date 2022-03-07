#A Procrustes test based upon two distance matrices.
#Assuming the matrices are n by n, the distances are first converted to n points in n-1 dimensions.
#This requires that the matrices be metric and Euclidean, so an additional function to check for that is included.
#Finally, a Procrustes test is performed.

#TODO: This performs a Procrustes test *with* scaling.
#The coclassification matrices have meaningful distances, bounded between 0 and 1.
#So we may want an unscaled Procrustes test.

distance_matrix_is_euclidean <- function(m) {
	return(ade4::is.euclid(as.dist(m)))
}

dist_to_coords <- function(m) {
	return(cmdscale(as.dist(m), nrow(m)-1))
}

procrustes_test <- function(m1, m2) {
	return(vegan::protest(dist_to_coords(m1), dist_to_coords(m2)))
}

procrustes_ss <- function(m1, m2) {
	return(procrustes_test(m1, m2)$ss)
}

procrustes_p <- function(m1, m2) {
	return(procrustes_test(m1, m2)$signif)
}
