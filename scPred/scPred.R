library(scPred)
library(methods)
# Single-cell transcriptomics and cell type information are both curated from literature and pre-processed generating ndata and the corresponding celltype objects.
# ndata represents the normolized dgCMatrix object using Seurat. Each row represents a gene (gene symbol) and each column represents a cell (cell barcode).  
# celltype represents the data.frame object containg two columns, namely cell barcode and cell type.

# ref_ndata represents the normolized dgCMatrix object using Seurat from HCL or MCA. Each row represents a gene (gene symbol) and each column represents a cell (cell barcode).  
# ref_celltype represents the data.frame object containg two columns, namely cell barcode and cell type.

getFeatureSpace_debug <- function(object, pVar, varLim = 0.01, correction = "fdr", sig = 0.05){
  
  
  # Validations -------------------------------------------------------------
  
  if(!is(object, "scPred") & !is(object, "Seurat")){
    stop("Invalid class for object: must be 'scPred' or 'Seurat'")
  }
  
  if(!any(correction %in% stats::p.adjust.methods)){
    stop("Invalid multiple testing correction method. See ?p.adjust function")
  }
  
  if(is(object, "scPred")){
    classes <- metadata(object)[[pVar]]
  }else{
    classes <- object[[pVar, drop = TRUE]]
  }
  
  if(is.null(classes)){
    stop("Prediction variable is not stored in metadata slot")
  }
  
  if(!is.factor(classes)){
    message("Transforming prediction variable to factor object...")
    classes <- as.factor(classes)
  }
  
  # Filter principal components by variance ---------------------------------
  
  if(is(object, "scPred")){ # scPred object
    
    # Get PCA
    i <- object@expVar > varLim
    pca <- getPCA(object)[,i]
    
    # Get variance explained
    expVar <- object@expVar
    
  }else{ # seurat object
    
    # Check if a PCA has been computed
    if(!("pca" %in% names(object@reductions))){
      stop("No PCA has been computet yet. See RunPCA() function")
    }
    
    # Check if available was normalized
    
    assay <- DefaultAssay(object)
    cellEmbeddings <- Embeddings(object)
    
    
    # Subset PCA
    expVar <- Stdev(object)**2/sum(Stdev(object)**2)
    names(expVar) <- colnames(Embeddings(object))
    i <-  expVar > varLim
    
    # Create scPred object
    pca <- Embeddings(object)[,i]
    
  }
  
  uniqueClasses <- unique(classes)
  isValidName <- uniqueClasses == make.names(uniqueClasses)
  
  if(!all(isValidName)){
    
    invalidClasses <- paste0(uniqueClasses[!isValidName], collapse = "\n")
    message("Not all the classes are valid R variable names\n")
    message("The following classes are renamed: \n", invalidClasses)
    classes <- make.names(classes)
    classes <- factor(classes, levels = unique(classes))
    newPvar <- paste0(pVar, ".valid")
    if(is(object, "scPred")){
      object@metadata[[newPvar]] <- classes
    }else{
      object@meta.data[[newPvar]] <- classes
    }
    message("\nSee new classes in '", pVar, ".valid' column in metadata:")
    message(paste0(levels(classes)[!isValidName], collapse = "\n"), "\n")
    pVar <- newPvar
  }
  
  
  
  
  # Select informative principal components
  # If only 2 classes are present in prediction variable, train one model for the positive class
  # The positive class will be the first level of the factor variable
  
  if(length(levels(classes)) == 2){
    
    message("First factor level in '", pVar, "' metadata column considered as positive class")
    res <- .getFeatures(levels(classes)[1], expVar, classes, pca, correction, sig)
    res <- list(res)
    names(res) <- levels(classes)[1]
    
  }else{
    
    res <- pblapply(levels(classes), .getFeatures, expVar, classes, pca, correction, sig)
    names(res) <- levels(classes)
    
  }
  
  
  nFeatures <- unlist(lapply(res, nrow))
  
  noFeatures <- nFeatures == 0
  
  if(any(noFeatures)){
    
    warning("\nWarning: No features were found for classes:\n",
            paste0(names(res)[noFeatures], collapse = "\n"), "\n")
    
    res1<- list()
    for (j in 1:length(res)) {
      d1<- res[[j]]
      if (nrow(d1) > 0) {
        res1[names(res)[j]]<- res[j]
      }
    }
    res<- res1
  }
  
  message("\nDONE!")
  
  
  # Assign feature space to `features` slot
  if(inherits(object, "Seurat")){
    
    # Create scPred object
    scPredObject <- list(expVar = expVar,
                         features = res,
                         pVar = pVar,
                         pseudo = FALSE)
    
    object@misc <- list(scPred = scPredObject)
    
  }else{
    
    
    object@features <- res
    object@pVar <- pVar
  }
  
  object
  
}

scPredict_debug <- function(object, newData = NULL, threshold = 0.9, 
                            returnProj = TRUE, returnData = FALSE, informative = TRUE,
                            useProj = FALSE){
  
  # Function validations ----------------------------------------------------
  
  # Validate if provided object is an scPred object
  if(!is(object, "scPred")){
    stop("'object' must be of class 'scPred'")
  }
  
  # Validate if scPred models have been trained already
  if(!length(object@train)){
    stop("No models have been trained!")
  }
  
  # Predictions are only possible if new data is provided. The new data can be a gene expression
  # matrix or a loading projection that has been computed independently.
  # Validate if new data is provided. If only a projection is found in the @projection slot,
  # this one is used as the prediction/test data
  if(is.null(newData) & (nrow(object@projection) == 0)){ # Neither newData nor projection
    
    stop("No newData or pre-computed projection")
    
  }else if(is.null(newData) & nrow(object@projection)){ # No newData and projection
    
    message("Using projection stored in object as prediction set")
    useProj <- TRUE
    
  }else if(!is.null(newData) & nrow(object@projection)){ # NewData and projection
    
    if (!(is(newData, "matrix") | is(newData, "Matrix"))){
      stop("'newData' object must be a matrix or seurat object")
    }
    message("newData provided and projection stored in scPred object. Set 'useProj = TRUE' to override default projection execution")
  }
  
  # Evaluate if there are identified features stored in the scPred object
  if(!length(object@features)){
    stop("No informative principal components have been obtained yet.\nSee getInformativePCs() function")
  }
  
  # Convert newData to sparse Matrix object
  if(is(newData, "matrix")){
    newData <- Matrix(newData)
  }
  
  # Data projection ---------------------------------------------------------
  
  # By default, projection of training loadings is automatically done. This option can be override
  # with the `useProj` argument to use an independent projection stored directly in the object.  
  if(!useProj){
    projection <- projectNewData(object = object,
                                 newData = newData,
                                 informative = informative, 
                                 seurat = if(!is.null(object@svd$seurat)){TRUE}else{FALSE})
  }else{
    projection <- object@projection
  }
  # Get cell classes used for training
  classes <- names(object@features)
  # Cell class prediction ---------------------------------------------------
  # For all cell classes which a training models has been trained for, get
  # the prediction probability for each cell in the `projection` object
  message("Predicting cell types")
  res <- pbapply::pblapply(classes, .predictClass, object, projection)
  names(res) <- levels(classes)
  # Exclude NA cekk type
  res_nrow <- NULL
  for (i in 1:length(res)) {
    res_nrow1 <- nrow(res[[i]])
    res_nrow <- c(res_nrow,res_nrow1)
    if (res_nrow1 != nrow(projection)) {
      d1 <- data.frame(res_nrow1 = rep(0,nrow(projection)))
      colnames(d1)<- colnames(res[[i]])
      res[[i]]<- d1
    }
  }
  res_nrow1<- which(res_nrow != nrow(projection))
  # Gather results in a dataframe
  res <- as.data.frame(res)
  if (length(res_nrow1) > 0) {
    res <- res[,-res_nrow1]
  }
  row.names(res) <- rownames(projection)
  if(length(classes) == 1){
    # If only one model was trained (positive and negative classes present in prediction
    # variable only), get the probability for the negative class using the complement rule
    # P(negative_class) = 1 - positive_class
    allClasses <-  unique(object@metadata[[object@pVar]])
    positiveClass <- classes
    negativeClass <- as.character(allClasses[!allClasses %in% classes])
    res[[negativeClass]] <- 1 - res[[positiveClass]]
    # Assign cells according to their associated probabilities
    res$predClass <- ifelse(res[,1] > threshold, positiveClass,
                            ifelse(res[,2] > threshold, negativeClass, "unassigned")) %>%
      as.factor()
    # Save results to `@predictions` slot
    object@predictions <- res
  }else{
    # If more than one model was trained (there are 3 classes or more in the prediction variable),
    # obtain the maximum probability for each cell and assign the respective class to that cell.
    i <- apply(res, 1, which.max)
    prob <- c()
    for(j in seq_len(nrow(res))){
      prob[j] <- res[j,i[j]]
    }
    predictions <- data.frame(res, probability = prob, prePrediction = names(res)[i])
    rownames(predictions) <- rownames(projection)
    predictions %>% 
      tibble::rownames_to_column("id") %>% 
      dplyr::mutate(predClass = ifelse(probability > threshold, as.character(prePrediction), "unassigned")) %>%
      dplyr::select(-probability, -prePrediction) %>% 
      tibble::column_to_rownames("id") -> finalPrediction
    
    # Save results to `@predictions` slot
    object@predictions <- finalPrediction
  }
  if(returnProj & !useProj){
    object@projection <- projection
  }
  if(returnData & !is.null(newData)){
    if(object@pseudo){
      object@predData <- log2(newData + 1)
    }else{
      object@predData <- newData
    }
  }
  return(object)
}



# Eigendecomposition
scp <- eigenDecompose(as.matrix(ref_ndata),pseudo = F)
metadata(scp) <- ref_celltype
# Feature selection
scp <- getFeatureSpace_debug(scp, pVar = "celltype")
# Model training
scp <- trainModel(scp)
# Prediction step
scp <- scPredict_debug(scp, newData = as.matrix(ndata))

scp_pred<- getPredictions(scp)
celltype$scPred<- scp_pred$predClass

