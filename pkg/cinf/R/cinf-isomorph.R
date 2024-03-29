# Searching isomorphisms

# Finds all substructure isomorphisms
find_substr_isomorph <- function (substr_lab, substr_ct, str_lab, str_ct) {
  
  isomorph <- function(stp) {
    for (i in 1:str_size)
      if (!used[i] && substr_lab[stp]==str_lab[i]) {
        cc[stp] <<- i
        if (stp > 1)
          for (j in 1:(stp-1))
            if ((substr_ct[stp,j]!=0) && (substr_ct[stp,j]!=str_ct[i,cc[j]])) { 
              to_exit <- TRUE
              next
            }
        if (to_exit) {
          to_exit <- FALSE
          next
        }
        if (stp == substr_size) {
          num_matches <<- num_matches + 1
          isom_list[[num_matches]] <<- cc
        }
        else {
          used[i] <<- TRUE
          isomorph(stp+1)
          used[i] <<- FALSE
        }
      }
  }

  substr_size <- length(substr_lab)
  str_size <- length(str_lab)
  used <- logical(str_size)
  cc <- integer(substr_size)
  to_exit <- FALSE
  num_matches <- 0
  isom_list <- list()
  isomorph(1)
  isom_list
} 

