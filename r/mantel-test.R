#ape's implementation of the Mantel test doesn't give the correlation; it gives an unnormalised test statistic.
#vegan's implementation doesn't give a p value; it gives a significance value with no control over the tail of the test.
#So here I've rolled my own implementation.

permute_matrix <- function(m) {
	s <- sample(1:nrow(m))
	return(m[s,s])
}

mantel_cor <- function(m1, m2) {
	#Grab the upper triangles of both, and get the correlation between them.
	return(cor(m1[upper.tri(m1)], m2[upper.tri(m2)]))
}

mantel_p <- function(m1, m2, nperm = 999, alternative = c("two.sided", "less", "greater")) {
	alternative <- match.arg(alternative)
	obs_cor <- mantel_cor(m1, m2)
	perm_cors <- replicate(nperm, mantel_cor(m1, permute_matrix(m2)))
	
	num_less <- sum(perm_cors <= obs_cor)
	num_greater <- sum(perm_cors >= obs_cor)
	num_counted <- switch(alternative,
		two.sided = 2*min(num_less, num_greater),
		less = num_less,
		greater = num_greater
	)
	p <- (num_counted+1)/(nperm+1)
	return(min(1, p))
}
