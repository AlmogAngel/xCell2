% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xcell2object-generics.R,
%   R/xcell2object-methods.R
\name{spillMat}
\alias{spillMat}
\alias{setSpillMat<-}
\alias{getSpillMat,xCell2Object-method}
\alias{getSpillMat}
\alias{setSpillMat<-,xCell2Object-method}
\title{Access Spillover Matrix}
\usage{
getSpillMat(object)

setSpillMat(object) <- value

\S4method{getSpillMat}{xCell2Object}(object)

\S4method{setSpillMat}{xCell2Object}(object) <- value
}
\arguments{
\item{object}{An \linkS4class{xCell2Object}.}

\item{value}{A matrix of spillover correction factors (for the setter).}
}
\value{
For `getSpillMat`, a matrix of spillover correction factors. 
        For `setSpillMat<-`, the updated \linkS4class{xCell2Object}.

For `getSpillMat`, a matrix of spillover correction factors. 
        For `setSpillMat<-`, the updated \linkS4class{xCell2Object}.
}
\description{
Retrieve or assign the spillover correction matrix for an \linkS4class{xCell2Object}.

Retrieve or assign the spillover correction matrix for an \linkS4class{xCell2Object}.
}
\examples{
data(DICE_demo.xCell2Ref, package = "xCell2")
spill_mat <- matrix(c(1, 0.1, 0.1, 1), nrow = 2, byrow = TRUE)
rownames(spill_mat) <- colnames(spill_mat) <- c("T_cells", "B_cells")
setSpillMat(DICE_demo.xCell2Ref) <- spill_mat
getSpillMat(DICE_demo.xCell2Ref)
data(DICE_demo.xCell2Ref, package = "xCell2")
spill_mat <- matrix(c(1, 0.05, 0.05, 1), nrow = 2, byrow = TRUE)
rownames(spill_mat) <- colnames(spill_mat) <- c("T_cells", "B_cells")
setSpillMat(DICE_demo.xCell2Ref) <- spill_mat
}
\seealso{
\link{xCell2Object-class}

\link{xCell2Object-class}
}
