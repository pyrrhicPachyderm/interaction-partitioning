#A set of helper functions for building distance matrices, for comparison with the coclassification matrices.

build_dist_matrix <- function(species_names, data, dist_func, force_symmetric = FALSE) {
	#Builds and annotates a distance matrix between species.
	#Requires a vector of species names, and a vector of data about the species.
	#dist_func must take two elements of the data vector, and return a distance.
	#All elements, even the diagonal, will be filled with dist_func, so it should return 0 if given the same data twice.
	#As usual for a distance function, it should be symmetric; if it isn't, this can be forced.
	#Forcing symmetry is only useful when it should be symmetric, but isn't (due to machine error or stochasticity),
	#as no guarantees are made about how symmetry is forced.
	
	m <- sapply(data, function(col_data){
		sapply(data, function(row_data){
			dist_func(col_data, row_data)
		})
	})
	if(force_symmetric) {
		m[lower.tri(m)] <- t(m)[lower.tri(m)]
	}
	rownames(m) <- species_names
	colnames(m) <- species_names
	return(m)
}
