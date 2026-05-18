#!/usr/bin/env Rscript
# ==============================================================================
# Cell Type Classification Functions
# ==============================================================================
# Collection of 15 different classifier functions for single-cell annotation
# Each function takes reference counts/labels and query counts
# Returns a list with predictions and optional scores
# ==============================================================================

suppressPackageStartupMessages({
  require(SingleR)
  require(SingleCellExperiment)
  require(scuttle)
  require(glmnet)
  require(xgboost)
  require(nnet)
  require(randomForest)
  require(e1071)
  require(class)
  require(rpart)
  require(MASS)
  require(Matrix)
  require(Seurat)
  require(scPred)
})

# ncores is provided by the calling environment (set in run_label_transfer_classifier via config)

# ============================================================================
# HELPER FUNCTION: LOG-NORMALIZE COUNTS (SPARSE-AWARE)
# ============================================================================
Log_Normalize <- function(mat, scale.factor = 1e4, margin = 2L) {
  # Compute total counts per cell (column) or gene (row)
  total.counts <- if (margin == 2L) Matrix::colSums(mat) else Matrix::rowSums(mat)
  
  # Ensure sparse matrix operations
  xnorm <- mat
  xnorm@x <- xnorm@x / rep(total.counts, diff(xnorm@p))  # Sparse division by column sums
  xnorm@x <- xnorm@x * scale.factor                       # Sparse scaling
  xnorm@x <- log1p(xnorm@x)                               # Sparse log1p
  
  return(xnorm)
}

# ============================================================================
# HELPER FUNCTION: COMPUTE ROW VARIANCE FOR SPARSE MATRICES
# ============================================================================
sparse_row_var <- function(mat) {
  # Compute row variance without converting to dense: var(X) = E[X²] - E[X]²
  row_means <- Matrix::rowMeans(mat)
  row_sq_means <- Matrix::rowMeans(mat^2)
  row_var <- row_sq_means - row_means^2
  return(row_var)
}

# ============================================================================
# HELPER FUNCTION: SELECT HIGHLY VARIABLE GENES
# ============================================================================
select_hvg <- function(counts, n_hvg = 2000) {
  # Log-normalize counts for variance calculation (keeps sparse)
  counts_normalized <- Log_Normalize(counts, scale.factor = 1e4, margin = 2L)
  
  # Calculate variance for each gene using sparse operations
  gene_vars <- sparse_row_var(counts_normalized)
  
  # Select top genes by variance
  n_hvg <- min(n_hvg, nrow(counts))
  top_genes <- order(gene_vars, decreasing = TRUE)[1:n_hvg]
  
  return(top_genes)
}

# ============================================================================
# 1. SINGLERR CLASSIFIER
# ============================================================================
classify_SingleR <- function(ref_counts, ref_labels, query_counts, method = "pearson") {
  res <- tryCatch({
    message("Running classifier: SingleR")
    # SingleR: use all genes (reference-based, needs full gene set)
    sce_ref <- SingleCellExperiment::SingleCellExperiment(list(counts = ref_counts))
    sce_query <- SingleCellExperiment::SingleCellExperiment(list(counts = query_counts))

    workers <- max(1L, as.integer(ncores))
    sce_ref <- scuttle::logNormCounts(sce_ref)
    sce_query <- scuttle::logNormCounts(sce_query)

    bpparam <- BiocParallel::MulticoreParam(workers = workers, progressbar = FALSE)
    message(sprintf("SingleR using %d BiocParallel workers", BiocParallel::bpnworkers(bpparam)))
    
    pred <- SingleR::SingleR(
      test = sce_query,
      ref = sce_ref,
      labels = ref_labels,
      BPPARAM = bpparam
    )
    
    pred <- as.data.frame(pred)
    as.character(pred$pruned.labels)
    
  }, error = function(e) {
    warning("SingleR failed: ", e$message)
    NA
  })
  res
}

# ============================================================================
# 2. LOGISTIC REGRESSION (using glmnet)
# ============================================================================
classify_LogisticRegression <- function(ref_counts, ref_labels, query_counts) {
  res <- tryCatch({
    message("Running classifier: LogisticRegression")
    # Select highly variable genes
    hvg_idx <- select_hvg(ref_counts, n_hvg = 2000)
    
    # Drop classes with <2 cells to satisfy multinomial requirements
    class_counts <- table(ref_labels)
    keep_classes <- names(class_counts)[class_counts >= 5]
    keep_cells <- ref_labels %in% keep_classes
    if (length(unique(ref_labels[keep_cells])) < 2) {
      stop("Not enough classes with >=2 cells for multinomial logistic regression")
    }
    
    # Log-normalize FULL gene set, then subset to HVGs (best practice)
    ref_norm_full <- Log_Normalize(ref_counts, scale.factor = 1e4, margin = 2L)
    query_norm_full <- Log_Normalize(query_counts, scale.factor = 1e4, margin = 2L)
    
    # Subset HVGs and class-filtered cells
    ref_norm <- ref_norm_full[hvg_idx, keep_cells, drop = FALSE]
    query_norm <- query_norm_full[hvg_idx, , drop = FALSE]
    
    X_train <- Matrix::t(ref_norm)
    X_test <- Matrix::t(query_norm)
    # drop unused factor levels
    y_train <- droplevels(as.factor(ref_labels[keep_cells]))
    
    # Filter out zero-variance features to avoid NaN in normalization
    sds_all <- apply(X_train, 2, sd)
    keep_features <- which(sds_all > 0)
    X_train <- X_train[, keep_features, drop = FALSE]
    X_test  <- X_test[, keep_features, drop = FALSE]
    
    # Z-score normalization
    means <- Matrix::colMeans(X_train)
    sds <- sds_all[keep_features]
    X_train <- sweep(sweep(X_train, 2, means, "-"), 2, sds, "/")
    X_test <- sweep(sweep(X_test, 2, means, "-"), 2, sds, "/")
    
    # Replace any remaining NaN/Inf with 0
    X_train[!is.finite(X_train)] <- 0
    X_test[!is.finite(X_test)] <- 0
    
    # Use cv.glmnet for automatic lambda selection (if enough samples per class)
    # cv.glmnet requires nfolds > 3, so only use it if minimum class size >= 4
    min_class_size <- min(table(y_train))
    
    # Use CV when we have enough samples per class
    nfolds <- min(5L, min_class_size)

    cv_fit <- cv.glmnet(
      X_train, y_train,
      family = "multinomial",
      alpha = 1,
      type.measure = "class",
      nfolds = nfolds
    )
    pred_prob <- predict(cv_fit, X_test, s = "lambda.min", type = "class")
    
    as.character(pred_prob[, 1])
    
  }, error = function(e) {
    warning("Logistic Regression failed: ", e$message)
    NA
  })
  res
}

# ============================================================================
# 3. XGBOOST CLASSIFIER
# ============================================================================
classify_XGBoost <- function(ref_counts, ref_labels, query_counts) {
  res <- tryCatch({
    message("Running classifier: XGBoost")
    # Log-normalize FULL gene set, then subset to HVGs
    ref_norm_full <- Log_Normalize(ref_counts, scale.factor = 1e4, margin = 2L)
    query_norm_full <- Log_Normalize(query_counts, scale.factor = 1e4, margin = 2L)
    
    # Use highly variable genes to improve signal-to-noise
    hvg_idx <- select_hvg(ref_counts, n_hvg = 2000)
    ref_counts2 <- ref_norm_full[hvg_idx, ]
    query_counts2 <- query_norm_full[hvg_idx, ]
    
    X_train <- Matrix::t(ref_counts2)
    X_test <- Matrix::t(query_counts2)
    
    label_encoding <- as.factor(ref_labels)
    y_train <- as.numeric(label_encoding) - 1
    label_map <- levels(label_encoding)
    
    dtrain <- xgb.DMatrix(X_train,
                          label = y_train)
    
    params <- list(
      objective = "multi:softprob",
      num_class = length(label_map),
      eta = 0.3,
      max_depth = 6,
      subsample = 0.8,
      colsample_bytree = 0.8,
      nthread = ncores
    )
    
    model <- xgb.train(params, dtrain,
                       nrounds = 100,
                       verbose = 0)
    pred_prob <- predict(model, xgb.DMatrix(X_test))
    pred_labels <- label_map[apply(matrix(pred_prob, ncol = length(label_map), byrow = TRUE), 1, which.max)]
    
    pred_labels
  }, error = function(e) {
    warning("XGBoost failed: ", e$message)
    NA
  })
  res
}

# ============================================================================
# 4. MULTILAYER PERCEPTRON (neural network)
# ============================================================================
classify_MLP <- function(ref_counts, ref_labels, query_counts) {
  res <- tryCatch({
    message("Running classifier: MLP")
    # Log-normalize FULL gene set, then subset to HVGs
    ref_norm_full <- Log_Normalize(ref_counts, scale.factor = 1e4, margin = 2L)
    query_norm_full <- Log_Normalize(query_counts, scale.factor = 1e4, margin = 2L)
    
    # Select highly variable genes (2000)
    hvg_idx <- select_hvg(ref_counts, n_hvg = 2000)
    
    ## 2. Design matrices (cells × genes)
    X_train <- ref_norm_full[hvg_idx, ] %>% 
      Matrix::t() %>% 
      as.matrix()
    X_test  <- query_norm_full[hvg_idx, ] %>% 
      Matrix::t() %>% 
      as.matrix()
    
    ## 3. Remove zero-variance genes (CRITICAL)
    vars <- matrixStats::colVars(X_train)
    keep <- vars > 0
    
    X_train <- X_train[, keep, drop = FALSE]
    X_test  <- X_test[, keep, drop = FALSE]
    
    ## 4. PCA (NO scaling here)
    n_pcs <- min(30L, ncol(X_train) - 1L)
    
    pca_fit <- stats::prcomp(
      X_train,
      center = TRUE,
      scale. = FALSE,
      rank. = n_pcs
    )
    
    X_train_pca <- pca_fit$x
    X_test_pca  <- scale(
      X_test,
      center = pca_fit$center,
      scale = FALSE
    ) %*% pca_fit$rotation[, seq_len(n_pcs), drop = FALSE]
    
    ## 5. Replace any non-finite values
    X_train_pca[!is.finite(X_train_pca)] <- 0
    X_test_pca[!is.finite(X_test_pca)] <- 0
    
    ## 6. Convert to data frames
    X_train_pca <- as.data.frame(X_train_pca)
    X_test_pca  <- as.data.frame(X_test_pca)
    
    # Ensure valid column names
    colnames(X_train_pca) <- make.names(colnames(X_train_pca))
    colnames(X_test_pca) <- make.names(colnames(X_test_pca))
    
    ## 7. Labels (drop unused levels and ensure proper format)
    y_train <- droplevels(as.factor(ref_labels))
    
    ## Convert labels to one-hot numeric matrix
    y_train_matrix <- stats::model.matrix(~ y_train - 1)
    
    ## Train neural network
    model <- nnet::nnet(
      x = X_train_pca,
      y = y_train_matrix,
      size = 25,
      maxit = 500,
      decay = 1e-4,
      MaxNWts = 1e6,
      trace = FALSE
    )
    
    # Predict probabilities
    pred_probs <- predict(model, X_test_pca)
    # Clean column names
    colnames(pred_probs) <- sub("^y_train", "", colnames(pred_probs))
    # Get predicted labels
    pred_labels <- colnames(pred_probs)[max.col(pred_probs)]
    as.character(pred_labels)
  }, error = function(e) {
    warning("MLP failed: ", e$message)
    NA
  })
  res
}

# ============================================================================
# 5. RANDOM FOREST CLASSIFIER
# ============================================================================
classify_RandomForest <- function(ref_counts, ref_labels, query_counts) {
  res <- tryCatch({
    message("Running classifier: RandomForest")
    # Log-normalize FULL gene set, then subset to HVGs
    ref_norm_full <- Log_Normalize(ref_counts, scale.factor = 1e4, margin = 2L)
    query_norm_full <- Log_Normalize(query_counts, scale.factor = 1e4, margin = 2L)
    
    # Use highly variable genes to reduce noise and runtime
    hvg_idx <- select_hvg(ref_counts, n_hvg = 2000)
    ref_counts2 <- ref_norm_full[hvg_idx, ]
    query_counts2 <- query_norm_full[hvg_idx, ]
    
    X_train <- as.data.frame(as.matrix(Matrix::t(ref_counts2)))
    X_test <- as.data.frame(as.matrix(Matrix::t(query_counts2)))
    y_train <- droplevels(as.factor(ref_labels))
    
    model <- randomForest(X_train, y_train, ntree = 100, maxnodes = 30)
    pred <- predict(model, X_test)
    
    as.character(pred)
  }, error = function(e) {
    warning("Random Forest failed: ", e$message)
    NA
  })
  res
}

# ============================================================================
# 6. SUPPORT VECTOR MACHINE (SVM)
# ============================================================================
classify_SVM <- function(ref_counts, ref_labels, query_counts) {
  res <- tryCatch({
    message("Running classifier: SVM")
    # Log-normalize FULL gene set, then subset to HVGs
    ref_norm_full <- Log_Normalize(ref_counts, scale.factor = 1e4, margin = 2L)
    query_norm_full <- Log_Normalize(query_counts, scale.factor = 1e4, margin = 2L)
    
    # Select highly variable genes (2000)
    hvg_idx <- select_hvg(ref_counts, n_hvg = 2000)
    
    ## Design matrices (cells × genes)
    X_train <- ref_norm_full[hvg_idx, ] %>% 
      Matrix::t() %>% 
      as.matrix()
    X_test  <- query_norm_full[hvg_idx, ] %>% 
      Matrix::t() %>% 
      as.matrix()
    
    ## Remove zero-variance genes
    vars <- matrixStats::colVars(X_train)
    keep <- vars > 0
    
    X_train <- X_train[, keep, drop = FALSE]
    X_test  <- X_test[, keep, drop = FALSE]
    
    ## PCA for dimensionality reduction (efficient for kernel methods)
    n_pcs <- min(30L, ncol(X_train) - 1L)
    
    pca_fit <- stats::prcomp(
      X_train,
      center = TRUE,
      scale. = FALSE,
      rank. = n_pcs
    )
    
    X_train_pca <- pca_fit$x
    X_test_pca  <- scale(
      X_test,
      center = pca_fit$center,
      scale = FALSE
    ) %*% pca_fit$rotation[, seq_len(n_pcs), drop = FALSE]
    
    ## Replace non-finite values
    X_train_pca[!is.finite(X_train_pca)] <- 0
    X_test_pca[!is.finite(X_test_pca)] <- 0
    
    ## Convert to data frames
    X_train_pca <- as.data.frame(X_train_pca)
    X_test_pca  <- as.data.frame(X_test_pca)
    
    ## Labels
    y_train <- droplevels(as.factor(ref_labels))
    
    ## Train SVM on PCA space
    model <- svm(X_train_pca, y_train, kernel = "radial", cost = 1, probability = TRUE)
    pred <- predict(model, X_test_pca)
    
    as.character(pred)
  }, error = function(e) {
    warning("SVM failed: ", e$message)
    NA
  })
  res
}

# ============================================================================
# 7. K-NEAREST NEIGHBORS (kNN)
# ============================================================================
classify_kNN <- function(ref_counts, ref_labels, query_counts, k = 15) {
  res <- tryCatch({
    message("Running classifier: kNN")
    # Log-normalize FULL gene set, then subset to HVGs
    ref_norm_full <- Log_Normalize(ref_counts, scale.factor = 1e4, margin = 2L)
    query_norm_full <- Log_Normalize(query_counts, scale.factor = 1e4, margin = 2L)
    
    # Select highly variable genes (2000)
    hvg_idx <- select_hvg(ref_counts, n_hvg = 2000)
    
    ## Design matrices (cells × genes)
    X_train <- ref_norm_full[hvg_idx, ] %>% 
      Matrix::t() %>% 
      as.matrix()
    X_test  <- query_norm_full[hvg_idx, ] %>% 
      Matrix::t() %>% 
      as.matrix()
    
    ## Remove zero-variance genes
    vars <- matrixStats::colVars(X_train)
    keep <- vars > 0
    
    X_train <- X_train[, keep, drop = FALSE]
    X_test  <- X_test[, keep, drop = FALSE]
    
    ## PCA for dimensionality reduction (best practice for distance-based methods)
    n_pcs <- min(30L, ncol(X_train) - 1L)
    
    pca_fit <- stats::prcomp(
      X_train,
      center = TRUE,
      scale. = FALSE,
      rank. = n_pcs
    )
    
    X_train_pca <- pca_fit$x
    X_test_pca  <- scale(
      X_test,
      center = pca_fit$center,
      scale = FALSE
    ) %*% pca_fit$rotation[, seq_len(n_pcs), drop = FALSE]
    
    ## Replace non-finite values
    X_train_pca[!is.finite(X_train_pca)] <- 0
    X_test_pca[!is.finite(X_test_pca)] <- 0
    
    ## Run kNN on PCA space
    pred <- knn(X_train_pca, X_test_pca, as.factor(ref_labels), k = k)
    as.character(pred)
  }, error = function(e) {
    warning("kNN failed: ", e$message)
    NA
  })
  res
}

# ============================================================================
# 8. NAIVE BAYES CLASSIFIER
# ============================================================================
classify_NaiveBayes <- function(ref_counts, ref_labels, query_counts) {
  res <- tryCatch({
    message("Running classifier: NaiveBayes")
    # Log-normalize FULL gene set, then subset to HVGs
    ref_norm_full <- Log_Normalize(ref_counts, scale.factor = 1e4, margin = 2L)
    query_norm_full <- Log_Normalize(query_counts, scale.factor = 1e4, margin = 2L)
    
    # Select highly variable genes (2000)
    hvg_idx <- select_hvg(ref_counts, n_hvg = 2000)
    
    ## Design matrices (cells × genes)
    X_train <- ref_norm_full[hvg_idx, ] %>% 
      Matrix::t() %>% 
      as.matrix()
    X_test  <- query_norm_full[hvg_idx, ] %>% 
      Matrix::t() %>% 
      as.matrix()
    
    ## Remove zero-variance genes
    vars <- matrixStats::colVars(X_train)
    keep <- vars > 0
    
    X_train <- X_train[, keep, drop = FALSE]
    X_test  <- X_test[, keep, drop = FALSE]
    
    ## PCA for dimensionality reduction (orthogonal features satisfy independence assumption)
    n_pcs <- min(30L, ncol(X_train) - 1L)
    
    pca_fit <- stats::prcomp(
      X_train,
      center = TRUE,
      scale. = FALSE,
      rank. = n_pcs
    )
    
    X_train_pca <- pca_fit$x
    X_test_pca  <- scale(
      X_test,
      center = pca_fit$center,
      scale = FALSE
    ) %*% pca_fit$rotation[, seq_len(n_pcs), drop = FALSE]
    
    ## Replace non-finite values
    X_train_pca[!is.finite(X_train_pca)] <- 0
    X_test_pca[!is.finite(X_test_pca)] <- 0
    
    ## Convert to data frames
    X_train_pca <- as.data.frame(X_train_pca)
    X_test_pca  <- as.data.frame(X_test_pca)
    
    ## Labels
    y_train <- droplevels(as.factor(ref_labels))
    
    ## Train Naive Bayes on PCA space
    model <- e1071::naiveBayes(X_train_pca, y_train)
    pred <- predict(model, X_test_pca)
    
    as.character(pred)
  }, error = function(e) {
    warning("Naive Bayes failed: ", e$message)
    NA
  })
  res
}

# ============================================================================
# 9. LINEAR DISCRIMINANT ANALYSIS (LDA)
# ============================================================================
classify_LDA <- function(ref_counts, ref_labels, query_counts) {
  res <- tryCatch({
    message("Running classifier: LDA")
    # Log-normalize FULL gene set, then subset to HVGs
    ref_norm_full <- Log_Normalize(ref_counts, scale.factor = 1e4, margin = 2L)
    query_norm_full <- Log_Normalize(query_counts, scale.factor = 1e4, margin = 2L)
    
    # Select highly variable genes
    hvg_idx <- select_hvg(ref_counts, n_hvg = 2000)
    ref_counts2 <- ref_norm_full[hvg_idx, ]
    query_counts2 <- query_norm_full[hvg_idx, ]
    
    X_train <- as.matrix(Matrix::t(ref_counts2))
    X_test <- as.matrix(Matrix::t(query_counts2))
    y_train <- droplevels(as.factor(ref_labels))
    
    pca_fit <- prcomp(X_train, rank. = min(30, ncol(X_train) - 1), scale. = TRUE)
    X_train_pca <- pca_fit$x
    X_test_pca <- scale(X_test, center = pca_fit$center, scale = pca_fit$scale) %*% pca_fit$rotation[, 1:ncol(X_train_pca)]
    
    df_train <- as.data.frame(X_train_pca)
    df_train$label <- y_train
    model <- MASS::lda(label ~ ., data = df_train)
    
    df_test <- as.data.frame(X_test_pca)
    pred <- predict(model, df_test)$class
    
    as.character(pred)
  }, error = function(e) {
    warning("LDA failed: ", e$message)
    NA
  })
  res
}

# ============================================================================
# 10. SEURAT TRANSFER DATA (CCA projection)
# ============================================================================
classify_SeuratTransfer <- function(ref_counts, ref_labels, query_counts) {
  res <- tryCatch({
    message("Running classifier: SeuratTransfer")
    # SeuratTransfer: use all genes (Seurat does FindVariableFeatures internally)
    ref_seurat <- Seurat::CreateSeuratObject(counts = ref_counts)
    query_seurat <- Seurat::CreateSeuratObject(counts = query_counts)
    
    ref_seurat$cell_type <- ref_labels
    
    ref_seurat <- Seurat::NormalizeData(ref_seurat, verbose = FALSE) %>%
      Seurat::FindVariableFeatures(verbose = FALSE) %>%
      Seurat::ScaleData(verbose = FALSE) %>%
      Seurat::RunPCA(verbose = FALSE)
    
    query_seurat <- Seurat::NormalizeData(query_seurat, verbose = FALSE) %>%
      Seurat::ScaleData(features = Seurat::VariableFeatures(ref_seurat), verbose = FALSE)
    
    anchors <- Seurat::FindTransferAnchors(reference = ref_seurat,
                                           query = query_seurat,
                                           dims = 1:30,
                                           verbose = FALSE)
    
    predictions <- Seurat::TransferData(anchorset = anchors,
                                        refdata = ref_seurat$cell_type,
                                        dims = 1:30,
                                        verbose = FALSE)
    
    pred_labels <- predictions$predicted.id
    as.character(pred_labels)
  }, error = function(e) {
    warning("Seurat Transfer failed: ", e$message)
    NA
  })
  res
}

# ============================================================================
# 11. DECISION TREE CLASSIFIER
# ============================================================================
classify_DecisionTree <- function(ref_counts, ref_labels, query_counts) {
  res <- tryCatch({
    message("Running classifier: DecisionTree")
    # Log-normalize FULL gene set, then subset to HVGs
    ref_norm_full <- Log_Normalize(ref_counts, scale.factor = 1e4, margin = 2L)
    query_norm_full <- Log_Normalize(query_counts, scale.factor = 1e4, margin = 2L)
    
    # Select highly variable genes (2000)
    hvg_idx <- select_hvg(ref_counts, n_hvg = 2000)
    
    ## Design matrices (cells × genes)
    X_train <- ref_norm_full[hvg_idx, ] %>% 
      Matrix::t() %>% 
      as.matrix()
    X_test  <- query_norm_full[hvg_idx, ] %>% 
      Matrix::t() %>% 
      as.matrix()
    
    ## Remove zero-variance genes
    vars <- matrixStats::colVars(X_train)
    keep <- vars > 0
    
    X_train <- X_train[, keep, drop = FALSE]
    X_test  <- X_test[, keep, drop = FALSE]
    
    ## PCA for dimensionality reduction
    n_pcs <- min(30L, ncol(X_train) - 1L)
    
    pca_fit <- stats::prcomp(
      X_train,
      center = TRUE,
      scale. = FALSE,
      rank. = n_pcs
    )
    
    X_train_pca <- pca_fit$x
    X_test_pca  <- scale(
      X_test,
      center = pca_fit$center,
      scale = FALSE
    ) %*% pca_fit$rotation[, seq_len(n_pcs), drop = FALSE]
    
    ## Replace non-finite values
    X_train_pca[!is.finite(X_train_pca)] <- 0
    X_test_pca[!is.finite(X_test_pca)] <- 0
    
    ## Convert to data frames
    X_train_pca <- as.data.frame(X_train_pca)
    X_test_pca  <- as.data.frame(X_test_pca)
    
    colnames(X_train_pca) <- make.names(colnames(X_train_pca))
    colnames(X_test_pca) <- make.names(colnames(X_test_pca))
    
    ## Labels
    y_train <- droplevels(as.factor(ref_labels))
    
    df_train <- X_train_pca
    df_train$label <- y_train
    model <- rpart(label ~ ., data = df_train, method = "class", maxdepth = 15)
    
    pred <- predict(model, X_test_pca, type = "class")
    as.character(pred)
  }, error = function(e) {
    warning("Decision Tree failed: ", e$message)
    NA
  })
  res
}

# ============================================================================
# 12. SCPRED CLASSIFIER
# ============================================================================
classify_scPred <- function(ref_counts, ref_labels, query_counts) {
  res <- tryCatch({
    message("Running classifier: scPred")
    # Create Seurat objects directly from raw counts
    # scPred handles all preprocessing internally (normalization, HVG, scaling, PCA)
    seurat_ref <- Seurat::CreateSeuratObject(counts = ref_counts, min.cells = 0, min.features = 0)
    seurat_query <- Seurat::CreateSeuratObject(counts = query_counts, min.cells = 0, min.features = 0)
    
    # Add cell type metadata to reference
    seurat_ref$cell_type <- ref_labels
    
    tbl <- table(ref_labels)
    keep <- names(tbl)[tbl >= 10]   # I strongly recommend ≥10
    
    seurat_ref <- subset(
      seurat_ref,
      subset = cell_type %in% keep
    )
    
    # Preprocess reference: normalize, find variable features, scale, PCA
    seurat_ref <- seurat_ref %>%
      Seurat::NormalizeData(verbose = F) %>%
      Seurat::FindVariableFeatures(verbose = F) %>%
      Seurat::ScaleData(verbose = F) %>%
      Seurat::RunPCA(npcs = 30,
                     verbose = FALSE)
    
    # Preprocess query: normalize only (scPred will align to reference)
    seurat_query <- seurat_query %>%
      Seurat::NormalizeData()

    # Train scPred models on reference
    seurat_ref <- scPred::getFeatureSpace(seurat_ref, "cell_type")
    seurat_ref <- scPred::trainModel(seurat_ref)
    
    # Predict on query using scPred
    seurat_query <- scPred::scPredict(seurat_query, seurat_ref)
    
    # Extract predictions from query metadata
    pred_labels <- seurat_query$scpred_prediction
    
    # Handle "unassigned" cases by returning them as-is
    as.character(pred_labels)
    
  }, error = function(e) {
    warning("scPred failed: ", e$message)
    NA
  })
  res
}

# ============================================================================
# 13. RANDOM CLASSIFIER (baseline/negative control)
# ============================================================================
classify_Random <- function(ref_counts, ref_labels, query_counts) {
  res <- tryCatch({
    message("Running classifier: Random")
    # Get unique cell types from reference
    cell_types <- unique(ref_labels)
    
    # Randomly assign cell types to query cells
    n_query <- ncol(query_counts)
    random_assignments <- sample(cell_types, size = n_query, replace = TRUE)
    
    as.character(random_assignments)
  }, error = function(e) {
    warning("Random classifier failed: ", e$message)
    NA
  })
  res
}

# ============================================================================
# 14. ENSEMBLE VOTING CLASSIFIER
# ============================================================================
classify_Ensemble <- function(predictions_list) {
  res <- tryCatch({
    message("Running classifier: Ensemble")
    n_cells <- length(predictions_list[[1]])
    
    ensemble_pred <- sapply(1:n_cells, function(i) {
      votes <- sapply(predictions_list, function(pred_vec) {
        if (all(is.na(pred_vec))) NA else pred_vec[i]
      })
      
      votes <- votes[!is.na(votes)]
      
      if (length(votes) == 0) {
        return(NA)
      }
      
      names(sort(table(votes), decreasing = TRUE))[1]
    })
    
    as.character(ensemble_pred)
  }, error = function(e) {
    warning("Ensemble voting failed: ", e$message)
    NA
  })
  res
}
