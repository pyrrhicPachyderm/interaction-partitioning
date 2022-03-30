#A function for printing a grouped alpha matrix as a table using knitr.
#Expects the species names, the row grouping, the column grouping, and the alpha matrix.
#The groupings are given as a vector, [number of species] in length, with the number of the group each species belongs to.
#The matrix is defined in terms of the groups, and hence is [number of row groups] by [number of column groups].
