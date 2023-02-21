#An R6 object to read, store and handle the priors file fed that was into the C++ program.
#Depends on strings.R.

Priors <- R6::R6Class("Priors",
	public = list(
		priors_table = NULL,
		
		initialize = function(priors_file_name) {
			self$priors_table <- read.table(priors_file_name, fill = TRUE, comment.char = "#")
		},
		
		#Intended for use in math mode.
		get_as_latex = function(row) {
			name <- self$priors_table[row,1]
			parameters <- self$priors_table[row,-1]
			parameters <- parameters[!is.na(parameters)]
			parameter_strings <- paste0("\\num[round-mode=figures,round-precision=3,round-pad=false]{", parameters, "}")
			string <- paste0("\\text{", name, "}\\left(", paste(parameter_strings, collapse = ","), "\\right)")
			return(string)
		}
	)
)
