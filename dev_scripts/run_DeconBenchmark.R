library(DeconBenchmark)

# Get the list of supported methods
allSupportMethods <- getSupportedMethods()

methodsToRun <- c("EPIC", "dtangle") # Select methods to run (must be in the list of supported methods)
requiredInputs <- getMethodsInputs(methodsToRun, containerEngine = "docker")
print(requiredInputs)


# Get all Tabula Sapiens simulations
dataFiles <- list.files('/bigdata/almogangel/twelve_years_decon_paper/analysis/data/sim')

# Get the first simulation (i=1)
i <- 1
dataFile <- paste0('/bigdata/almogangel/twelve_years_decon_paper/analysis/data/sim/', dataFiles[i])
data <- readRDS(dataFile)
reference <- generateReference(data$singleCellExpr, data$singleCellLabels, type="signature") # Generate reference
deconvolutionResults <- runDeconvolution(methodsToRun, bulk = data$bulk, singleCellExpr = data$singleCellExpr,
                                         singleCellLabels = data$singleCellLabels, signature=reference$signature,
                                         sigGenes = data$sigGenes, cellTypeExpr = data$cellTypeExpr) # Run benchmarking



# test args for xCell2.0
DeconBenchmark::.writeArgs(h5file = "/bigdata/almogangel/twelve_years_decon_paper/docker/xCell2/test_args.h5", bulk = data$bulk, singleCellExpr = data$singleCellExpr)
