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


#' Value Negating
#'
#' %nin% is the opposite of the binary operator %in%. It returns a logical vector
#' indicating if there is a mismatch (TRUE) or a match (FALSE). In other words, it
#' returns TRUE for the values in the first vector but not in the second.
#' @param x vector or NULL: the values to be matched.
#' @param y vector or NULL: the values to be matched against.
#' @keywords same
#' @export
#' @examples
#' x = 1:6
#' y = 4:10
#' x %nin% y # which values of x are not in y?
#' # TRUE  TRUE  TRUE FALSE FALSE FALSE
'%nin%' = Negate('%in%')

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
#' str2lfmm recodes biallelic SNP loci in STRUCTURE format to a format suitable for LFMM.
#' An arbitrary vector of reference alleles is automatically chosen, though the reference
#' can be specified by the user if known. str2lfmm accepts both diploid and haploid SNP data.
#' Each diploid individual must be coded on two lines. Ambiguity codes must be removed and
#' nucleotides must be coded as integers (e.g. A, T, C, G as 1, 2, 3, 4). Trait data (either
#' phenotype or environmental variables) can also be included. Samples may be excluded by
#' specifying using 'exclude,'
#' @param str.data SNP data in STRUCTURE format
#' @param trait.data Phenotype or environment data
#' @param exclude Vector of sample IDs to exclude
#' @keywords lfmm, structure, recode
#' @export
#' @examples
#' str2lfmm(str.data = str.data,
#' trait.data = trait.data,
#' exclude = NULL)
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

#' Calculate 2D and Cumulative Distance Between Coordinates
#'
#' Appends two columns to a data frame containing Longitude and Latitude coordinates. The first column (Dist) is the
#' incremental change in distance between coordinates. The second column (DistTotal) is the cumulative distance along a
#' 2D transect beginning at 0.
#' @param x data frame containing coordinate data as columns labeled "Longitude" (x) and "Latitude" (y).
#' @param longlat logical for calculating Great Circle (TRUE) or Euclidean (FALSE) distances. Defaults to TRUE.
#' Great Circle distance is returned in km.
#' @keywords great circle, dist
#' @export
#' @details Relies on the spDistsN1 function from the package 'sp.'
#' @examples
#' df_dist <- dist_calc(df, longlat = TRUE)
dist_calc <- function(df, longlat = TRUE){
  ## Test if the package 'sp' is attached
  if("sp" %nin% (.packages()) & "sp" %in% .packages(all.available = TRUE)){
    library(sp)
  }
  if("sp" %nin% (.packages()) & "sp" %nin% .packages(all.available = TRUE)){
    warning("Please install package:sp using install.packages(\"sp\")")
  }
  Dist <- 0
  for(i in 2:length(df[,"Longitude"])) {
    Dist[i] = spDistsN1(as.matrix(df[i,c("Longitude", "Latitude")]),
                        c(df[,"Longitude"][i-1], df[,"Latitude"][i-1]),
                        longlat = TRUE) # longlat so distances will be in km,
    #great circle distance
    #longlat = TRUE) / 1.609 # longlat so distances will be in km,
    #then divide to convert to miles
  }

  DistTotal <- 0
  for(i in 2:length(df[,"Longitude"])) {
    DistTotal[i] = Dist[i] + DistTotal[i-1]
  }
  new_df <- cbind(df, Dist, DistTotal)
  return(new_df)
}
