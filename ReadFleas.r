   # Set the file to be read
(d.file <- paste(data.dir, "flea.dat", sep=""))
(d.column <- paste(data.dir, "flea.col", sep=""))

headers <- scan(d.column, "" )
n.var <- length(headers)
d.flea <- matrix(scan(d.file), ncol=n.var, byrow=T)
flea.species <- c(rep("C",21),rep("Hp",22),rep("Hk",31))
species <- c(rep(1,21),rep(2,22),rep(3,31))

dimnames(d.flea) <- list(1:length(flea.species),headers)

df.flea <- data.frame(d.flea)
cat("New variables\n")
cat("  headers      - variable names\n")
cat("  n.var        - number of variables\n")
cat("  d.flea       - data matrix\n")
cat("  flea.species - C,Hp,Hk\n")
cat("  species      - 1,2,3\n")
cat("  df.flea      - data frame\n")

