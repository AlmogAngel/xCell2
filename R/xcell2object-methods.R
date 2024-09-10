# Getter and Setter for Signatures
setGeneric("getSignatures", function(object) standardGeneric("getSignatures"))
setMethod("getSignatures", "xCell2Object", function(object) object@signatures)

setGeneric("setSignatures<-", function(object, value) standardGeneric("setSignatures<-"))
setReplaceMethod("setSignatures", "xCell2Object", function(object, value) {
  object@signatures <- value
  object
})

# Getter and Setter for Dependencies
setGeneric("getDeps", function(object) standardGeneric("getDeps"))
setMethod("getDeps", "xCell2Object", function(object) object@dependencies)

setGeneric("setDeps<-", function(object, value) standardGeneric("setDeps<-"))
setReplaceMethod("setDeps", "xCell2Object", function(object, value) {
  object@dependencies <- value
  object
})

# Getter and Setter for Params
setGeneric("getParams", function(object) standardGeneric("getParams"))
setMethod("getParams", "xCell2Object", function(object) object@params)

setGeneric("setParams<-", function(object, value) standardGeneric("setParams<-"))
setReplaceMethod("setParams", "xCell2Object", function(object, value) {
  object@params <- value
  object
})

# Getter and Setter for Spillover Matrix
setGeneric("getSpillMat", function(object) standardGeneric("getSpillMat"))
setMethod("getSpillMat", "xCell2Object", function(object) object@spill_mat)

setGeneric("setSpillMat<-", function(object, value) standardGeneric("setSpillMat<-"))
setReplaceMethod("setSpillMat", "xCell2Object", function(object, value) {
  object@spill_mat <- value
  object
})

# Getter and Setter for Genes Used
setGeneric("getGenesUsed", function(object) standardGeneric("getGenesUsed"))
setMethod("getGenesUsed", "xCell2Object", function(object) object@genes_used)

setGeneric("setGenesUsed<-", function(object, value) standardGeneric("setGenesUsed<-"))
setReplaceMethod("setGenesUsed", "xCell2Object", function(object, value) {
  object@genes_used <- value
  object
})
