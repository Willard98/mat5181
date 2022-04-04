disp.str <- function(str) {
   for (i in 1:length(str)) {
      cat(str[i]);cat("\n")
   }
}

disp.str.g <- function(str) {
#   if (dev.list()[[1]] !> 1) 
   for (i in 1:length(str)) {
      cat(str[i]);cat("\n")
   }
}

