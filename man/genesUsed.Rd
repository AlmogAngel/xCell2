% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xcell2object-generics.R,
%   R/xcell2object-methods.R
\name{genesUsed}
\alias{genesUsed}
\alias{setGenesUsed<-}
\alias{getGenesUsed,xCell2Object-method}
\alias{getGenesUsed}
\alias{setGenesUsed<-,xCell2Object-method}
\title{Access Genes Used}
\usage{
getGenesUsed(object)

setGenesUsed(object) <- value

\S4method{getGenesUsed}{xCell2Object}(object)

\S4method{setGenesUsed}{xCell2Object}(object) <- value
}
\arguments{
\item{object}{An \linkS4class{xCell2Object}.}

\item{value}{A character vector of genes (for the setter).}
}
\value{
For `getGenesUsed`, a character vector of genes. 
        For `setGenesUsed<-`, the updated \linkS4class{xCell2Object}.

For `getGenesUsed`, a character vector of genes. 
        For `setGenesUsed<-`, the updated \linkS4class{xCell2Object}.
}
\description{
Retrieve or assign the genes used in training the reference for an \linkS4class{xCell2Object}.

Retrieve or assign the genes used in training the reference for an \linkS4class{xCell2Object}.
}
\examples{
data(DICE_demo.xCell2Ref, package = "xCell2")
setGenesUsed(DICE_demo.xCell2Ref) <- c("GeneA", "GeneB", "GeneC")
getGenesUsed(DICE_demo.xCell2Ref)
data(DICE_demo.xCell2Ref, package = "xCell2")
setGenesUsed(DICE_demo.xCell2Ref) <- c("GeneA", "GeneB", "GeneC")
}
\seealso{
\link{xCell2Object-class}

\link{xCell2Object-class}
}
