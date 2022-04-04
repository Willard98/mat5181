

# Setting a path ----------------------------------------------------------

# Let's say you want to make a path to a particular file
# and you aren't sure how R will read it because it changes based on operating system.
# We can use the file.path() function to add the right separators for your operating
# system by using the .Platform$file.sep function
data.dir <- file.path("F:", "Storage", "School", "4601_F2020", "L1", "Data", .Platform$file.sep)
code.dir <- file.path("F:", "Storage", "School", "4601_F2020", "L1", "Code", .Platform$file.sep)

# Please set the data directory and code directory to
# the folder you plan to store your code and data for this course!
# You might also set the working directory there using setwd("enterPathHere")


# Setting your working directory ------------------------------------------

# To figure out where your directory is right now, type in console
getwd()
# To set it to something else:
setwd(code.dir)

# Be careful to consider the side effects of changing your working directory:
# Relative file references in your code (for datasets, source files, etc)
# will become invalid when you change working directories.
# The way paths are constructed depends on your platform
# If you're using Windows, you use a / in R instead of the \ that Windows uses
# so you need to change any path you paste or type into R.

# Type into your Console to find out what to use on your computer:
.Platform$file.sep


# Load the Fleas Dataset --------------------------------------------------

# Be careful here. If you built your data.dir as above,
# it already has the final / built into it and does not need
# a "sep = .Platform$file.sep" argument
# Instead, specify that there is no separator using sep = ""

# Specify the path to the flea.dat file - the data
(d.file <- paste(data.dir, "flea.dat", sep = ""))
# Specify the path to the flea.col - the column headers,
# which are in a separate file here
(d.column <- paste(data.dir, "flea.col", sep = ""))

# Load the column headers.
# There are plenty of "read.xyz" functions
# read.csv()
# read.table()
# scan()
# as well as many others from various packages

# Let's try scan()
# Notice that it fails! And it gives us a weird error.
# What does this mean? It's expecting "real" or numeric values
# but it got character value "tars1" instead
(headers <- scan(d.column))

# Instead, specify "character" or an empty string
(headers <- scan(d.column, what = "character"))
(headers <- scan(d.column, what = ""))

# An alternative is to read the column headers as a table
(headers <- read.table(d.column))
# Then convert them to a vector
# Notice that you have to specify that you are converting the first column of the table
# using headers$V1 or headers[, 1]
(headers <- as.character(headers$V1))

# The length of the vector gives the number of variables
(n.var <- length(headers))

# If we use read.table() here, it works well with this particular data
# without any other specifications like number of rows or columns
# or byrow = TRUE
# This is NOT always the case! Often, you must specify that it is
# tab-separated data, newline separated data, comma-separated data, and so on.
# To do so in read.table, use sep = "", or sep = "," etc.

# Suppose you use scan()
# Notice the problem?
d.flea.s <- scan(d.file)
typeof(d.flea.s)
class(d.flea.s)
length(d.flea.s)

# It reads in the entire dataset in a single vector!
# Let's try "wrapping" the scan() function
# in the matrix() function. The first argument
# in the matrix() function is "data = something".
# So let's pass the scan() function on our dataset as that argument.
# Notice the problem again?
d.flea.s <- matrix(scan(d.file))
# It did the same thing, but now it's a matrix with one column!
d.flea.s[1:5,]

# Let's try again, but with the scan() function we specify character as before
# and with the matrix() function we specify the number of columns that we got earlier.
# Then we tell it to fill in the data by row.
d.flea.str <- matrix(scan(d.file, "character"), ncol = n.var, byrow = TRUE)
d.flea.str[1:5,]
typeof(d.flea.str)
class(d.flea.str)
dim(d.flea.str)
length(d.flea.str)

# Now we've got our data in a tabular matrix format, but the numbers are in quotation marks
# so they're actually characters (like we asked to happen in the scan() function)
# Why might we do this? Let's imagine that some of the data has been entered as "Missing"
# or as some other word. What's going to happen if we read in the entire dataset as numeric values?
# It's likely to mess up!

# Now, take a look at the first value in the table.
# How can we make this into numeric data?
# If we use the as.numeric() function, that value becomes a number,
# but notice that it doesn't stay that way. We haven't told R to save over the value.
d.flea.str[1, 1]
as.numeric(d.flea.str[1,1])
d.flea.str[1, 1]

# Let's try taking the d.flea.str character-type data and replacing d.flea.s with it.
# We'll do it by using the as.numeric() function to make it into numeric data
# and putting it in the "data" parameter of the matrix() function.
d.flea.s <- matrix(as.numeric(d.flea.str), ncol = n.var)
d.flea.s[1:5,]
str(d.flea.s)
class(d.flea.s)
typeof(d.flea.s[, 1])

# What if we want to skip all of this work?
# Let's use the read.table() function instead.
d.flea <- read.table(d.file)
d.flea[1:5,]
str(d.flea)
class(d.flea)

# The read.table function automatically made our data into a nice table!
# In fact, we now have what is called a data.frame, which is a useful object for a lot of things.
# Notice that a data.frame is also a list. The columns of a dataframe are stored in lists.
# If we want to extract a column, we can use the number of the column, or we can use the name
# of the column in quotation marks. Remember to use double square brackets with lists: [[]]
# If we want to look at specific rows of that column, add another set of single square brackets.
# We can also use single square brackets with matrix/data.frame notation.
# And we can also use $ notation with the variable name.
typeof(d.flea)
names(d.flea)
View(d.flea[["V1"]])
View(d.flea[, "V1"])
View(d.flea[, 1])
d.flea[["V1"]][1:5]
d.flea[1:5, "V1"]
d.flea$V1[1:5]
head(d.flea)
tail(d.flea)

# To look at a specific group of elements, we need to get fancier.
# Recall that we can concatenate a group of things together using c().
# The dim() function gives us the dimensions of the table: rows x columns
dim(d.flea)
# But how do we use those? Let's say we want the number of rows in a table.
dim(d.flea)[1:2]
# Let's assign that to a variable called numRows. Be careful what you use for a variable name!
# If you assign a variable name to the same name as a function it will replace the function
# until you restart R, so you will lose easy access to that function until then.
# We could also use the nrow() and ncol() functions.
(numRows <- dim(d.flea)[1])
(numRows <- nrow(d.flea))
(numCols <- dim(d.flea)[2])
(numCols <- ncol(d.flea))
# Say we want to look at rows that are multiples of 5 and columns that are multiples of 2.
# We'll make two sequences.
(fives <- seq(from = 5, to = numRows, by = 5))
(twos <- seq(from = 2, to = numCols, by = 2))
# Since R is a functional programming language, we could skip a step.
(fives <- seq(from = 5, to = nrow(d.flea), by = 5))
(twos <- seq(from = 2, to = ncol(d.flea), by = 2))

# Let's take a look at only those elements!
d.flea.s[fives, twos]
# We might also skip another step.
d.flea.s[seq(from = 5, to = nrow(d.flea), by = 5), seq(from = 2, to = ncol(d.flea), by = 2)]
# Compare to the original
d.flea.s[5, 2]
d.flea.s[c(5,10,15,20), c(2,4,6)]
d.flea[c(5,10,15,20), c(2,4,6)]
View(d.flea)


# Rename Column Headers/Names ---------------------------------------------

# Now let's get rid of those pesky default column headers/names.
dim(d.flea.s)
# Recall that d.flea.s is a matrix
class(d.flea.s)
names(d.flea.s)
# We have to use the dimnames() function to replace row and column names in a matrix.
# Notice that currently the dimnames are null.
dimnames(d.flea.s)
# We use a list to replace it. The first is the rownames, and we will number the rows.
# The second element will be the list of column headers.
print(dim(d.flea.s)[1])
print(1:dim(d.flea.s)[1])
print(headers)
dimnames(d.flea.s) <- list(1:dim(d.flea.s)[1], headers)
# Success!
d.flea.s[c(5,10,15,20), c(2,4,6)]

# How can we do that with the columns of a data.frame?
class(d.flea)
# We have two options. The colnames() function, or the names() function.
colnames(d.flea)
names(d.flea)
# Either way, place your data.frame within the function and assign our headers to it
# The item you're assigning must be less than or equal to the number of columns
# or you will get an error. Typically, you want to put a character vector here.
colnames(d.flea) <- headers
names(d.flea) <- headers
names(d.flea) <- c("tars1", "tars2", "head", "aede1", "aede2", "aede3")
# What if you just want to replace one name? Say "head" should be "tail" instead.
(names(d.flea))
names(d.flea)[3] <- "tail"
names(d.flea)
d.flea[c(5,10,15,20), ]
# Let's fix that again.
names(d.flea)[3] <- "head"


# Labelled Data -----------------------------------------------------------

# This is a well-known dataset, so we know how many of each type of flea beetle is in this data
# Let's create some labels. First, make a sequence of 74 values.
# There are three types of beetle: "Concinna" (C), "Heptapot" (Hp), and "Heikert" (Hk) beetles.
flea.species <- c(rep("C",21),rep("Hp",22),rep("Hk",31))
View(flea.species)
# We will also make a sequence of numbers to go with them.
species <- c(rep(1,21),rep(2,22),rep(3,31))
View(species)

# We have row names already in a file. Point to the correct path.
(d.row <- paste(data.dir, "flea.row", sep = ""))
# Scan them in. Notice that they come in with quotation marks.
(row.headers <- scan(d.row, ""))
# Let's tell it not to use quotation marks by using the noquote() function.
# Notice that row.headers is still a character object, though.
(row.headers <- noquote(scan(d.row, "")))
typeof(row.headers)

# Make sure we're using a data.frame
is.data.frame(d.flea)
class(d.flea)
# Just to be sure, let's create a separate object called df.flea that is a data.frame
df.flea <- data.frame(d.flea.s)
# What happens if we try to add row.headers as our row names, though?
# We get an error. Rownames must be unique. They cannot be repeated values like we have here.
rownames(df.flea) <- row.headers


# Running code from a file ------------------------------------------------

# A few extra functions have been created for us to use in this next section.
# Ask R to load them. We do not need to change them to make these functions work,
# but they will help automatically show some things later.
source(paste(code.dir, "DispStr.r", sep = ""))
source(paste(code.dir, "pairs_ext.r", sep = ""))
source(paste(code.dir, "MakeStereo.r", sep = ""))



# Esquisse  ---------------------------------------------------------------

# install.packages("esquisse")
esquisse::esquisser()


# Rattle ------------------------------------------------------------------
# install.packages("rattle")
rattle::rattle()
library(rattle)
rattle()
d.flea$species <- species
rattle()


# R Commander -------------------------------------------------------------

install.packages("Rcmdr")
library(Rcmdr)


# Visualization of Pairs of Variables -------------------------------------

# We can look at a scatterplot matrix of our data values
# The pairs function makes a scatterplot of each pair of variables
# and displays a table of the results. If you are using RStudio,
# the Plots frame will display it.
# Notice that the plots are symmetric, so the top-right triangle is a reflection
# of the lower-left triangle.
# The graphics.off() function closes a plot. If you are using R without R studio
# this is useful for closing the current plots you have open.
# In RStudio, it will clear the plots in your Plots frame.
graphics.off()
pairs(d.flea)

# This can be made more informative by displaying the histograms of the variables
# and the correlations between them.
# Note that the size of the correlation determines the size of the number.
# Remember that correlations range between -1 and 1 and describe a linear relationship
# between two variables.
# Here we will use the pairs_ext.r file, which has added panel.cor and panel.hist
# to R for us to use.
graphics.off()
pairs(d.flea, upper.panel = panel.cor, diag.panel = panel.hist,
      main = "Pairs Plot of Flea Species with Diagonal Histogram and Correlations")

# If we want to debug the functions here, we could use the debug() function
# which will open the pairs_ext.r file and let us skim through the code.
# This is especially useful if we have tried to make a change to the input of pairs
# alongside using panel.cor and panel.hist and it has caused an error.
# However, this might be a little complicated if you are new to R or to programming,
# so feel free to skip this for now and come back when you are more comfortable!

# debug(panel.cor)
# graphics.off()
# pairs(d.flea, upper.panel = panel.cor, diag.panel = panel.hist)
# undebug(panel.cor)


# Let's try to make our plot look a little nicer.
# The pch parameter lets us change the size or type of points in the scatterplot.
# Possible numbers go from 1 to 25.
# The cex parameter changes the size of your points. Compare using cex = 1.5 to cex = 0.7.
# The xaxt and yaxt parameters = "n" means the axis labels are removed.
graphics.off()
pairs(d.flea, col = species + 1, pch = 19, cex = 1.5,
      xaxt = "n", yaxt = "n",
      main = "Pairs Plot of Flea Species with Colour and Big Points")
pairs(d.flea, col = species + 1, pch = 19, cex = 0.7,
      xaxt = "n", yaxt = "n",
      main = "Pairs Plot of Flea Species with Colour and Small Points")
# We could also make it so that the shape is different based on species instead of the colour
pairs(d.flea, col = "cornflowerblue", pch = species + 5, cex = 0.75,
      xaxt = "n", yaxt = "n",
      main = "Pairs Plot of Flea Species with Point-Shape based on Colour")
# Or both shape and colour could vary based on species
pairs(d.flea, col = species + 1, pch = species + 5, cex = 0.75,
      xaxt = "n", yaxt = "n",
      main = "Pairs Plot of Flea Species with Colour")
# Or we could scale the points based on species, since we've used 1, 2, and 3 for species labels
pairs(d.flea, col = species + 1, pch = 19, cex = 0.5 * species,
      xaxt = "n", yaxt = "n",
      main = "Pairs Plot of Flea Species with Colour")

# We can limit which pairs to include.
# The y ~ x1 + x2 + ... formula notation can be used as the first parameter
# Leave the "y" value  blank in this case, and specify data to be the dataset
names(d.flea)
pairs(~ aede1 + aede2 + aede3, data = d.flea,
      col = species + 1, pch = 19,
      xaxt = "n", yaxt = "n",
      main = "Pairs Plot of Flea Species Subset using Formula")
# Alternatively, use square bracket notation to select columns in the x parameter.
pairs(x = d.flea[, 1:4],
      col = species + 1, pch = 19,
      xaxt = "n", yaxt = "n",
      main = "Pairs Plot of Flea Species Subset using Square Brackets")
# Or use the c() function to select specific columns
pairs(x = d.flea[, c(1,3,5)],
      col = species + 1, pch = 19,
      xaxt = "n", yaxt = "n",
      main = "Pairs Plot of Flea Species Subset using concatenation function")
# We can also change the diagonal labels
pairs(d.flea[, c(1,3,5)],
      col = species + 1, pch = 19,
      xaxt = "n", yaxt = "n",
      labels = c("Beetle Tars 1", "Beetle Head", "Beetle Aede 2"),
      main = "Pairs Plot of Flea Species Subset using Custom Labels")



# Using a Library ---------------------------------------------------------

# The ggplot2 package is considered by many to be one of the best visualization
# packages available in R, alongside packages like plotly, lattice, and a few others.
# The Base R plotting functions are extremely powerful, so they are all you need.
# But ggplot2 offers an alternative "language" to plot in, and many people
# have built libraries to augment ggplot2.
# ggplot2 can take some time to learn, but it is a very useful package.
# # The main author of ggplot2 has made a free tutorial book available online at
# https://ggplot2-book.org/
# Install the ggplot2 library and the GGally library, then load them.
# You only have to run install.packages() once, unless you update R to a new version
# or otherwise alter the program substantially
# You do need to load it each time you re-open R!
# If your R session crashes, you must also reload your packages.
install.packages("ggplot2")
install.packages("GGally")
library("ggplot2")
library("GGally")

# The ggpairs plot is similar to the pairs plot, but uses ggplot instead.
# It also builds in the correlation along the top half
# and density plots along the diagonal, much like we saw above!
ggpairs(d.flea)

# The essence of ggplot is the use of geometric layers and aesthetics, one on top of another.
ggpairs(d.flea, mapping = ggplot2::aes(colour = as.factor(species)),
        title = "GGPairs Plot", lower = list(continuous = wrap("smooth", alpha = 0.3, size = 0.1)),
        )
# What does all of this mean?
# The mapping = aes() parameter is calling to ggplot2 instead of GGally.
# It wants to set the "aesthetic" of the entire plot to use species to decide the colour.
# The default ggplot function is typically like this:
# ggplot(datasetName, aes(x, y, colour))
# We could extend it easily using + notation. Say we call our plot gp.
# gp <- ggplot(datasetName, aes(x, y, colour)) + geom_point() + theme_bw()
# The "geom_point() function makes the basic ggplot a scatterplot
# The theme_bw() function uses a pre-made black&white theme for the plot.
# Within geom_point() we could add another aes() function
# that modifies the aesthetics of the scatterplot points.
# We could add different items to change almost any part of our plot.
# There are hundreds of options for legends, axis labels, points, lines, and so on.
# Here is a more complicated version of the same plot.
ggpairs(d.flea, title = "GGPairs Plot",
               mapping = ggplot2::aes(colour = as.factor(species)),
               lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.1),
                            discrete = "blank", combo = "blank"),
               diag = list(discrete="barDiag",
                           continuous = wrap("densityDiag", alpha = 0.5 )),
               upper = list(combo = wrap("box_no_facet", alpha = 0.5),
                            continuous = wrap("cor", size = 4))) +
   theme(panel.grid.major = element_blank())    # remove gridlines

# The entire ggpairs function has many parameters
# In the "upper" and "lower" parameters we can choose boxplots ("box_no_facet") for
# combined continuous/discrete data or discrete data only can use "facetbar".
# Correlations can be used for continuous data.
# We can leave it blank by choosing "na"
# The bottom panel can use "points" for continuous,
# a facet histogram ("facethist") for combined data
# or a facet bar ("facetbar") for discete data.
# Again, we can leave it blank by choosing "na"
# Along the diagonal we have a choice of "densityDiag" for density plots for continuous data.
# We have "barDiag" for bar plots for discrete data, and "naDiag" for leaving it blank.
# If we place "alpha" between 0 and 1 in the aes() function, it makes the values semi-transparent.
# There is a lot more you can learn to do with these plots.
ggpairs(
   data = d.flea,
   mapping = ggplot2::aes(colour = as.factor(species), alpha = 0.4),
   columns = 1:6,
   title = "Smaller GGPairs Plot",
   lower = list("points"),
   xlab = "X Axis Title",
   ylab = "Y Axis Title",
   axisLabels = "show"
)


# Coplots in Base R -------------------------------------------------------

# Conditioning plots allow us to create plots that display one variable against another
# conditional on a third (and fourth...). For example, plot y against x, but
# make it conditional on two more variables, a and b.
# The variables must be numeric or factor variables.
# Here we have y = aede3, x = tars1, and a = aede1.
# This is a way to visualize higher dimensional data by fixing the values
# of some of the variables and visualizing the others in two dimensions.
graphics.off()
df.flea <- data.frame(d.flea[, 1:6])
plot(y = d.flea$aede3, x = d.flea$tars1)
coplot(aede3 ~ tars1 | aede1, data = df.flea)

# Overlap sets how much the conditioning variables overlap each other.
# When the conditioning variable is an ordered categorical variable with many levels,
# we can group the values of the conditioning variable into bins and allow
# this overlap. The intervals are displayed in a sensible order. The default is 6 intervals.
# The overlap can be seen in the overlap bars of the plot of the conditional variable.
# The typical plot parameters like pch (type of point), cex (size of points), and col (colour)
# can be used as usual.
# Let's try changing the overlap and watch how this changes the plotted points and the overlap bars.
# We set up a vector of possible overlap values and loop through them to re-run the plot with each one.
graphics.off()
ol <- c(0.1, 0.2, 0.3, 0.4, 0.5)
for(i in 1:5){
coplot(aede3 ~ tars1 | aede1, data = df.flea,
       overlap = ol[i],  col = species + 1, pch = 16, cex = 0.8)
}
# We can set the overlap in both the x and y directions.
# We can also set multiple conditioning variables. Notice that we use * to indicate these.
# The type = "p" parameter may be changed to type = "smooth" to use a loess (local polynomial) smoother
# when the displayed data is more linear. Other options like
# col.line = "YourColour" display the smoothing line in a chosen colour,
# and lwd = YourNumber set the line thickness.
# Note that we can use both points and a smoothing line simultaneously using
# c("p", "smooth") though it is not always sensible to do so.
graphics.off()
coplot(aede3 ~ tars1 | aede1 * head, data = df.flea,
       overlap = c(0.5, 0.25), col = species + 1, pch = 16, cex = 0.8,
       type = "p")


# xyplots -----------------------------------------------------------------

# Coplots in the lattice package

library(lattice)

# We will create an extra function to use alongside our next plot type,
# which is the version of the coplot from the lattice package.
# Be sure to install and load the lattice package!
# This function allows us to set fixed interval lengths for the plot
equal.space <- function(data, count) {
      # range(data) gives the max and min of the variable data.
      # diff takes the difference between the two values so
      # diffs gives the width of each interval.
   diffs <- diff(range(data))/count
      # min(data)+diffs*(0:(count-1)) gives the starting values
      #    for the intervals.
      # min(data)+diffs*(1:count) gives the ending values
      #    for the intervals.
      # cbind treats two(or more) vectors as column vectors
      #    and binds them as columns of a matrix.
   intervals <- cbind(min(data)+diffs*(0:(count-1)),
                      min(data)+diffs*(1:count))
      # shingle takes the interval structure and the data
      #    and breaks the data into the appropriate groups.
   return (shingle(data, intervals))
}

# equal.count() converts data into a shingle
# A shingle is a type of data structure that generalizes factors into
# pseudo-continuous variables
# The first plot uses equal cases in each grouping
C1 <- equal.count(df.flea$aede1, number = 6, overlap = 0.1)
xyplot(aede3 ~ tars1 | C1, data = df.flea, pch = 19)
# The second plot uses equal spacing in each grouping
C2 <- equal.space(df.flea$aede1, 6)
xyplot(aede3 ~ tars1 | C2, data = df.flea, pch = 19)

# A smiilar plot can be built in ggplot, but it is more difficult
# to deal with the overlaps. More complicated options and additions
# could be used to adjust this plot until it closely resembled the
# default xyplot. This means we have more control over customization
# but comes at the expense of needing to know how to built it in the
# first place.

ggplot(df.flea, aes(x = tars1, y = aede3)) +
   geom_point(size = 0.5) +
   facet_wrap(~ cut_number(C1, 6))


# Ellipse Outline Function ------------------------------------------------

# Import the ellipseOutline.r file
source(paste(d.R <- code.dir, "ellipseOutline.r", sep=""))


# Synthetic Data ----------------------------------------------------------

# We want to generate a synthetic dataset.
# Create an empty object
ec.t1 <- {}

# Using a for loop with t running from -20 up to 20 we fill that object
# by running the ellipse.outline() function using the value of t on each loop
for (t in -20:20) {
   ec.t1 <- rbind(ec.t1, cbind(ellipse.outline(20, 20, 10, 5, t, 0, (200-t^2)/10), t))
}
# This results in a matrix with 25,994 rows and 4 columns
str(ec.t1)
class(ec.t1)
typeof(ec.t1)
dim(ec.t1)
# Convert it to a data.frame
ec.t1 <- data.frame(ec.t1[sample(dim(ec.t1)[1], dim(ec.t1)[1]),])
# Now run a pairs() scatterplot matrix
# This might crash R if your system is a bit on the older side
pairs(ec.t1, upper.panel = panel.cor, diag.panel = panel.hist)

# An alternative is again the GGally::ggpairs() function.
# Note that if you are using RStudio, make sure your Plots
# panel is made large enough enough that it can display the plot.
# Otherwise it may return an error instead of displaying the plot.
ggpairs(
   data = ec.t1,
   mapping = ggplot2::aes(alpha = 0.4),
   columns = 1:4,
   title = "Plot of ec.t1",
   lower = list("points"),
   xlab = "X Axis Title",
   ylab = "Y Axis Title",
   axisLabels = "show"
)


# Extract Equal Spacings from Synthetic Dataset ---------------------------

X <- equal.space(ec.t1$x, 25)
Y <- equal.space(ec.t1$y, 25)
Z <- equal.space(ec.t1$z, 25)
T1 <- equal.space(ec.t1$t, 25)

# R tries to maximize the use of the plotting window, but
# this can be less helpful if the data has a shape because it could
# obscure that shape or squish/stretch it
# Use the aspect parameter to force R to use equal axis scales
# In base R, the plot() function allows you to set the axis scales:
# Use the xlim = c(xmin, xmax) and ylim = c(ymin, ymax) parameters.
# In ggplot, use + xlim(min, max) and + ylim(min, max) to do the same thing.
xyplot(z ~ x | Y, data = ec.t1, pch=".", main ="z ~ x | Y",
                  aspect = diff(range(ec.t1$z))/diff(range(ec.t1$x)))

# The x11() function opens a separate window to view visualizations in.
# This is an older graphical display. It works well with the
# basic R programming environment, but can sometimes cause trouble
# with RStudio, so try to troubleshoot and if you can't get it working
# contact your TA for help.
# If you use a Mac, you must install XQuartz first. Go to
# https://www.xquartz.org/
# and download the free software. Once enabled, be sure to set permissions
# to allow R to access XQuartz.
# Note that x11() has additional arguments you may you. Type ?x11 in your
# console to look at them.
x11()
# Once the x11 window is open, run the next line of code and click on the window
# to view the plot. If you close the x11 window the plot will not appear.
# You may also simply run the plot without the x11 window open
# to view it in RStudio's Plots pane.

# The next few plots show the synthetic data on other axes, from the "top down",
# for example. Each one uses a different conditioning variable.
xyplot(y ~ x | Z, data = ec.t1, pch=".", main ="y ~ x | Z",
                  aspect = diff(range(ec.t1$y))/diff(range(ec.t1$x)))
x11()
xyplot(z ~ y | X, data = ec.t1, pch=".", main ="z ~ y | X",
                  aspect = diff(range(ec.t1$z))/diff(range(ec.t1$y)))

x11()
xyplot(z ~ x | T1, data = ec.t1, pch=".", main ="z ~ x | T",
                  aspect = diff(range(ec.t1$z))/diff(range(ec.t1$x)))

Z5 <- equal.space(ec.t1$z, 5)
T5 <- equal.space(ec.t1$t, 5)
x11()
xyplot(z ~ x | T5*Z5, data = ec.t1, main ="z ~ x | T5*Z5", pch=".",
                   aspect = diff(range(ec.t1$z))/diff(range(ec.t1$x)))

# We can create a grid of plots showing various cases and conditions
r <- 1
c <- 1
for (i in -20:15) {               # Loop through i from -20 to 15
   ind <- ec.t1$t==i              # Get the cases for which the t value == i
   X <- ec.t1$x[ind]              # And the corresponding x,y,z values
   Y <- ec.t1$y[ind]
   Z <- ec.t1$z[ind]
      # In the following - ( ?cloud)
      # print - displays the
      # cloud - a function that creates a cloud of points,
      #         with xlim, ylim, zlim (the range of values on the axes)
      #         set to the maximum range (x) to give proper scaling.
      # subpanel - the function use to plot the points.
      # groups - allows classes to be identified.
      # screen - sets the viewpoint.
      # split  - c(col, row, cols, rows)
      # more   -
   print(cloud(Z ~ X*Y, xlim = range(ec.t1$x),
                ylim = range(ec.t1$x),zlim = range(ec.t1$x),
                subpanel = panel.superpose, groups = rep(1,dim(ec.t1)[1]),
                screen = list(z = 10, x = -80, y = 0), data = ec.t1),
                split = c(c, r, 6, 6), more = TRUE)
   c <- c+1
   if (c%%6 == 1) {   # Remainder mod 6
      c <- 1
      r <- r+1
   }
}



# High Dimensional Plotting -------------------------------------------------------------

# There are more advanced ways of doing high dimensional
# plotting than making a static grid.


# tourr  -----------------------------------------------------------

# Since GGobi is no longer maintained, you may download the tourr package instead.
# This is also by Cook and Wickham, and was published in 2019.
# tourr provides an update for some of the GGobi methods.
# https://www.jstatsoft.org/article/view/v040i02/v40i02.pdf
# provides a lengthy discussion of the features of tourr.
# https://cran.r-project.org/web/packages/tourr/tourr.pdf
# provides the Cran package outline with descriptions of the functions available.
# We will use tourr here but leave the GGobi code (below) in case you want to try it.

# Run the ReadFleas.r file.
# You may need to check the path to make sure the code of
# ReadFleas.r does not include sep = "/" if you have followed along with data.dir
# and code.dir instructions used earlier (or as below)
# Notice that your packages are NO LONGER installed! This is because you changed
# versions of R. You will have to install packages you need for each version.
# This can lead to some difficulties if the default Cran version is not compatible
# with the version of R you are using, but there are ways to work around that problem.
data.dir <- file.path("F:", "Storage", "School", "4601_F2020", "L1", "Data", .Platform$file.sep)
code.dir <- file.path("F:", "Storage", "School", "4601_F2020", "L1", "Code", .Platform$file.sep)
species <- c(rep(1,21),rep(2,22),rep(3,31))
source(paste(code.dir, "ReadFleas.r", sep=""))

# Install the tourr package and the ash package (needed for a method in tourr)
# install.packages("ash")
# install.packages("tourr")
library(tourr)

# You may crash occasionally. That might mean you need to reload everything
# from the start of this section onward (including packages!).
# When you're done viewing an animated plot (in x11 or otherwise),
# hit Esc in the Console (if using RStudio) to stop the animation.
# You may use some of the usual plotting parameters.
# Use a lower max_frames and fps to avoid crashing your system.

# animate_scatmat() opens with a scatterplot matrix of projections and animates them
# to rotate around in the number of dimensions specified in the d parameter.
# Some methods can be very intensive for your computer to run, so beware!
x11()
animate_scatmat(df.flea, grand_tour(d = 3), col = species + 1, cex = 1.3,
                max_frames = 20, fps = 1)

# There are the following tour methods. Try them out!

# 1. grand_tour()
graphics.off()
x11()
animate(df.flea,
        tour_path = grand_tour(),
        display = display_xy(),
        fps = 10)

graphics.off()
x11()
animate_dist(df.flea)

graphics.off()
x11()
animate_xy(df.flea)

graphics.off()
x11()
animate_pcp(df.flea)

graphics.off()
x11()
animate_pcp(df.flea, grand_tour(4))

# 2. local_tour()
graphics.off()
x11()
animate_xy(df.flea[, 1:3], tour_path = local_tour(start = basis_random(3, 2), angle = pi/4))

graphics.off()
x11()
animate_xy(flea[, 1:3], local_tour(basis_init(3, 2)))

graphics.off()
x11()
animate_xy(flea[, 1:3], local_tour(basis_init(3, 2), 0.2))

graphics.off()
x11()
animate_xy(flea[, 1:3], local_tour(basis_random(3, 2), 0.2))


# 3. little_tour()
graphics.off()
x11()
animate(df.flea,
        tour_path = little_tour(),
        display = display_xy(),
        fps = 10)

graphics.off()
x11()
animate_xy(df.flea, little_tour())

graphics.off()
x11()
animate_pcp(df.flea, little_tour(3))

graphics.off()
x11()
animate_scatmat(df.flea, little_tour(3))

graphics.off()
x11()
animate_pcp(df.flea, little_tour(4))

# 4. guided_tour()
graphics.off()
x11()
animate_xy(df.flea[, 1:3], guided_tour(holes()), sphere = TRUE)

graphics.off()
x11()
animate_xy(df.flea[, 1:6], guided_tour(holes()), sphere = TRUE)

graphics.off()
x11()
animate_dist(df.flea[, 1:6], guided_tour(holes(), 1), sphere = TRUE)

graphics.off()
clrs <- c("#486030", "#c03018", "#f0a800")
col <- clrs[as.numeric(flea$species)]
x11()
animate_xy(df.flea[, 1:6], guided_tour(lda_pp(flea[,7])), sphere = TRUE, col=col)

# 5. dependence_tour()
graphics.off()
x11()
animate_xy(flea[, 1:3], dependence_tour(c(1, 2, 2)))

graphics.off()
x11()
animate_xy(flea[, 1:4], dependence_tour(c(1, 2, 1, 2)))

graphics.off()
x11()
animate_pcp(flea[, 1:6], dependence_tour(c(1, 2, 3, 2, 1, 3)))

# 6. frozen_tour()
frozen <- matrix(NA, nrow = 4, ncol = 2)
frozen[3, ] <- .5
graphics.off()
x11()
animate_xy(df.flea[, 1:4], frozen_tour(2, frozen))

frozen <- matrix(NA, nrow = 5, ncol = 2)
frozen[3, ] <- .5
frozen[4, ] <- c(-.2, .2)
graphics.off()
x11()
animate_xy(df.flea[, 1:5], frozen_tour(2, frozen))

# 7. frozen_guided_tour()
frozen <- matrix(NA, nrow = 4, ncol = 2)
frozen[3, ] <- .5
graphics.off()
x11()
animate_xy(flea[, 1:4], frozen_guided_tour(frozen, holes()))



# Parallel Coordinates Display --------------------------------------------

# Parallel coordinates are good for looking at
# high dimensional data.
# With the tour made active, select a [Transient]
# brush and a color.
# Notice how all the displays respond to the brushing.
# This enables us to identify corresponding values.
str(ozone)
graphics.off()
x11()
animate_image(ozone)
display(g[1], "Parallel Coordinates Display")

graphics.off()
x11()
animate_pcp(flea[, 1:6], dependence_tour(c(1, 2, 3, 2, 1, 3)))



# Stereo ------------------------------------------------------------------

source(paste(code.dir, "MakeStereo.r", sep=""))
   #Stereo
   # Use the 3 variables that seemed useful
make.Stereo(d.flea[,c(1,5,6)], species, Main="Flea beetles", asp="F",
           Xlab="tars1", Ylab="aede2", Zlab="aede3")

make.Stereo(d.flea[,c(6, 5, 1)], species, Main="Flea beetles", asp="F",
           Zlab="tars1", Ylab="aede2", Xlab="aede3")

# From tourr:
animate_stereo(flea[, 1:6])



# RGL ---------------------------------------------------------------------

# Make sure you install the rgl package
library(rgl)

# Now let's take a look at some interactive 3D plots.
plot3d(df.flea$tars1, df.flea$aede3, df.flea$aede1, col = species + 1)

plot3d(d.flea[,1], d.flea[,5], d.flea[,6], xlab="tars1",ylab="aede2",zlab="aede3",
       col=species+1, size=0.5, type="s")

for (j in seq(0, 90, 10)) {
   for(i in 0:360) {
      rgl.viewpoint(i, j);
   }
}



# GGobi -------------------------------------------------------------------

# GGobi and tourr -------------------------------------------------------------------

# GGobi has not been maintained since around 2012, it still can work
# under the right circumstances. There is a standalone version available at
# http://ggobi.org/ for Windows, Mac, and Linux. Both the standalone and package
# versions tend to be a little unstable due to their age.
# # GGobi provides some helpful ways to visualize higher-dimensional data
# Professor Di Cook at Monash University in Australia created it alongside her
# student Hadley Wickham (who later created ggplot and many other very popular R packages)

# Since GGobi is so out of date now, it no longer works with the current version of R
# (R version 4.0.2 for Fall 2020). If you would like to try it out, you may run a
# previous version of R. It is recommended that you try something around version 2.x or earlier.
# You can use your preferred search engine, such as AskJeeves, Bing,
# or Altavista to find a previous version of R, or the Windows versions are available at
# https://cran.r-project.org/bin/windows/base/old/
# and the Mac versions are available at
# https://cran.r-project.org/bin/macosx/old/
# Choose one that is earlier than R 4.0 and install it (3.5 is not early enough!).
# It will not replace your current version. They may be run alongside each other.
# If you are using RStudio, restart it. Once you have done so, go to:
# Tool -> Global Options -> General -> Under "R Version" click "Change.
# A list of your installed R versions should appear.
# Select the earlier version. You must restart RStudio again.

data.dir <- file.path("F:", "Storage", "School", "4601_F2020", "L1", "Data", .Platform$file.sep)
code.dir <- file.path("F:", "Storage", "School", "4601_F2020", "L1", "Code", .Platform$file.sep)
source(paste(code.dir, "ReadFleas.r", sep=""))

install.packages("rggobi")
library(rggobi)


g <- ggobi(d.flea)

# This opens with a scatterplot of the projections on
# one plane - i.e. a part of the scatterplot matrix.
# We can get a scatterplot matrix in Ggobi
display(g[1], "Scatterplot Matrix")
# or - [Display][New scatterplot matrix]
# from the Ggobi console.
# (This may not display all pairs by default
#  - you will need to select the other variables.
# A view that gives a better look at the data is selected by
#  - [ViewMode][2D Tour]

display(g[1], "2D Tour")
# This gives a 2D "tour" of the 6 dimensional data.
# The portion of each variable in the view is shown by the
# representation of the axis in the bottom corner
# (and on the console).
# As the tour runs 3 clusters will appear.
# When they do, you can click [Pause] and apply brushing
# - [Interaction][Brush]

# As we move the brush, the data points change colour as we
# pass over them but return to normal when we move away.
# To change this behaviour, select
# [Transient] and change to [Persistent]


(old.col <- glyph_colour(g[1]))

# It turns out that we know the species of the
# flea beetles so we can compare the clustering that
# we observed with the true classification.
(noquote(rbind(flea.species, old.col)))

# Set the points to a single colour and style.
glyph_colour(g[1]) <-rep(2, 74)
glyph_type(g[1]) <-rep(4, 74)
# Another use of brushing is in linked plots,
# so we can open other displays.
display(g[1], "Parallel Coordinates Display")
# Parallel coordinates are good for looking at
# high dimensional data.
# With the tour made active, select a [Transient]
# brush and a color.
# Notice how all the displays respond to the brushing.
# This enables us to identify corresponding values.

# In order to proceed with other aspects of linking
# we will use the colours based on the species.
glyph_colour(g[1]) <- c(rep(6,21),rep(4,22),rep(9,31))

cols <- rep(6, 74)
cols[which(d.flea[,6] < 95)] <- 9
cols[which(d.flea[,1] < 160)] <- 4
glyph_colour(g[1]) <- cols

close(g)

source(paste(code.dir, "MakeStereo.r", sep="/"))

#Stereo
# Use the 3 variables that seemed useful
make.Stereo(d.flea[,c(1,5,6)], species, Main="Flea beetles", asp="F",
            Xlab="tars1", Ylab="aede2", Zlab="aede3")

make.Stereo(d.flea[,c(6, 5, 1)], species, Main="Flea beetles", asp="F",
            Zlab="tars1", Ylab="aede2", Xlab="aede3")

# rgl
library(rgl)

plot3d(d.flea[,1], d.flea[,5], d.flea[,6], xlab="tars1",ylab="aede2",zlab="aede3",
       col=species+1, size=0.5, type="s")

for (j in seq(0, 90, 10)) {
   for(i in 0:360) {
      rgl.viewpoint(i, j);
   }
}


# Randu
# The Randu random number data set.
# First look at the pairs data.

d.file <- paste(data.dir, "randu.dat", sep = "/")
d.randu <- read.table(d.file)
pairs(d.randu, upper.panel=panel.cor, diag.panel=panel.hist)

library(rggobi)
g <- ggobi(d.randu)
# Use a [View][2D tour]
close(g)

#====== Prim7 ============
prim <- read.table(paste(data.dir, "prim7.dat",sep="/"))
g <- ggobi(prim)

new.col <- rep(1, 500)
col.2 <- c(2,3,4,14,15,16,17,18,21,23,30,34,37,41,43,46,49,50,53,54,55,
           57,58,63,65,66,69,70,72,73,74,75,77,78,79,85,86,88,90,91,92,
           93,94,95,99,100,102,104,105,106,107,109,110,113,114,116,120,
           121,124,125,126,127,129,130,133,139,140,141,143,145,147,150,
           152,153,157,158,159,160,161,164,166,169,172,175,176,177,178,
           180,185,194,195,198,200,203,204,209,210,211,212,218,219,220,
           222,223,226,228,229,233,234,236,238,240,242,244,245,246,248,
           249,252,253,257,259,263,264,265,266,267,269,270,273,277,278,
           280,281,282,283,284,286,292,294,296,297,300,305,310,311,314,
           315,317,323,331,332,333,334,335,341,342,343,346,351,356,359,
           360,361,362,365,370,372,374,375,377,378,379,380,383,386,388,
           389,390,391,393,397,398,400,402,403,405,407,408,413,414,415,
           417,418,419,420,425,427,428,429,430,432,433,434,436,437,438,
           440,444,445,447,448,452,453,454,455,456,463,465,467,470,471,
           473,476,477,478,480,481,482,484,485,487,488,489,490,491,494,
           497)
col.3 <- c(11,20,27,33,47,51,60,61,62,98,115,118,119,132,155,186,191,
           193,202,205,207,208,213,225,230,231,232,235,239,243,250,251,
           268,272,295,312,316,338,339,345,349,354,358,364,366,376,381,
           395,401,421,422,446,460,496)
col.5 <- c(5,8,13,19,26,32,39,48,56,71,81,96,111,136,137,144,149,156,
           162,165,188,199,201,216,255,262,274,279,289,291,301,320,322,
           326,327,329,344,348,353,363,367,369,384,399,404,406,411,423,
           441,442,443,469,474,479,483,495,499,500)
col.8 <- c(7,29,31,36,89,101,117,131,138,154,173,187,190,192,196,197,
           206,247,254,256,258,287,290,298,299,309,324,325,385,387,464)
col.9 <- c(1,12,22,24,25,44,45,52,64,83,103,108,122,123,134,135,146,151,
           167,168,170,174,179,181,184,221,224,237,261,271,285,293,304,
           306,307,308,319,328,337,352,355,357,368,396,410,424,426,435,
           439,449,451,458,461,462,466,472,475,493)


new.col[col.2] <-  2
glyph_colour(g[1]) <-  new.col

new.col[col.3] <-  3
glyph_colour(g[1]) <-  new.col

new.col[col.5] <-  5
glyph_colour(g[1]) <-  new.col

new.col[col.8] <-  8
glyph_colour(g[1]) <-  new.col

new.col[col.9] <-  9
glyph_colour(g[1]) <-  new.col

prim.lin <- read.table(paste(data.dir, "prim7.lines",sep="/"))
edges(g[1]) <- prim.lin

close(g)

g <- ggobi(paste(data.dir, "cube6.xml",sep="/"))
close(g)




