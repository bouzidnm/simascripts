# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#' Test if All Vecter Elements are Equal
#'
#' same(x) is a simple function to check if all of the elements of a vector are equal. NA's are ignored using na.omit().
#' @param x Vector to be tested.
#' @keywords same
#' @export
#' @examples
#' x <- rep(3, times = 10)
#' same(x) # returns TRUE
#' y <- 1:10
#' same (y) # returns FALSE
same <- function(x){
  if (length(na.omit(unique(x)))==1) return(TRUE)
  #if (any(is.na(x))) print("At least one element is NA")
  else(return(FALSE))
}

#' Prepare Input Files for Latent Factor Mixed Modeling (LFMM)
#'
#' same(x) is a simple function to check if all of the elements of a vector are equal. NA's are ignored using na.omit().
#' @param x Vector to be tested.
#' @keywords same
#' @export
#' @examples
#' x <- rep(3, times = 10)
#' same(x) # returns TRUE
#' y <- 1:10
#' same (y) # returns FALSE
str2lfmm <- function(str.data, trait.data, exclude) {
  ## exclude excluded samples
  if (missing(exclude)) {
    str.data <- str.data
  }
  if (!missing(exclude)) {
    str.data <- str.data[-(which(str.data$V1 == exclude)), ]
    trait.data <-
      trait.data[-(which(trait.data$CATALOG_NUM == exclude)), ]
  }

  ## Make empty matrix to store genotype data
  var.data <-
    matrix(
      data = NA,
      nrow = nrow(str.data) / 2,
      ncol = ncol(str.data) - 1
    )

  ## Identify reference alleles
  alleles <- apply(str.data[, -1], 2, function(s) {
    na.omit(unique(s))
  })
  alleles <- as.data.frame(alleles)
  reference <- alleles[1, ]

  ## Populate matrix with variant calls
  IDs <- unique(str.data$V1)
  for (i in 1:length(IDs)) {
    tmp.df <- str.data[which(str.data$V1 == IDs[i]), -1]
    for (j in 1:ncol(str.data) - 1) {
      # if(any(is.na(unique(tmp.df[,j])))){
      #   var.data[i,j] <- NA
      # }
      if (same(tmp.df[, j])) {
        if (unique(tmp.df[, j]) == reference[j]) {
          var.data[i, j] <- 0
        }
        else
          (var.data[i, j] <- 2)
      }
      if (!same(tmp.df[, j]) & !any(is.na(unique((tmp.df[, j]))))) {
        var.data[i, j] <- 1
      }
    }
  }

  ## Subset trait data
  trait.sub <-
    trait.data[which(trait.data$CATALOG_NUM %in% str.data$V1), ]
  all.equal(trait.sub$CATALOG_NUM, unique(str.data$V1)) # sanity check
  trait.out <- trait.sub[, pmatch(trait, colnames(trait.sub))]

  ## Combine var.data and trait.sub
  lfmm.data <- list()
  lfmm.data$genotype <- var.data
  lfmm.data$environment <- trait.out
  return(lfmm.data)
}
