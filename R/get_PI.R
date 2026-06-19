########## Annotation-Based Variant-Importance Score Estimation ##########
#
# This file implements get_PI() for estimating variant-importance scores
# from functional annotations. The function supports both LASSO and GLM models
# and uses an ensemble approach with multiple models fit on different
# control samples to improve prediction stability.

#################### Main Function for PI Estimation ####################

#' Estimate Variant-Importance Scores from Functional Annotations
#'
#' @description
#' Estimates variant-importance scores (PI) based on functional
#' annotations using supervised learning. The function trains an ensemble of
#' prediction models (LASSO or GLM) on case and control annotations and applies
#' them to testing data to estimate the relative importance of each variant
#' for the phenotype of interest.
#'
#' @param training_caseAnnotation Numeric matrix (N_case x M) of functional
#'   annotation scores for case variants (known causal variants from training
#'   data). Each row represents one variant, each column represents one
#'   annotation feature.
#' @param training_controlAnnotation Numeric matrix (N_control x M) of functional
#'   annotation scores for control variants (known non-causal variants from
#'   training data). Must have the same number of columns as training_caseAnnotation.
#' @param training_control_need_N Integer, number of control variants to sample
#'   for each model in the ensemble. Should be less than or equal to
#'   nrow(training_controlAnnotation). Typically chosen to balance the training
#'   dataset (e.g., equal to or 2-3 times nrow(training_caseAnnotation)).
#' @param model_need_N Integer, number of models to fit in the ensemble. Higher
#'   values increase prediction stability at the cost of computation time.
#'   Typical values: 5-20.
#' @param modelType Character string specifying the model type: "LASSO" or "GLM".
#'   \itemize{
#'     \item "LASSO": L1-regularized logistic regression with automatic feature
#'       selection via cross-validation. Recommended when M is large or features
#'       are correlated.
#'     \item "GLM": Standard logistic regression using all features. Recommended
#'       when M is small relative to N_case and features are relatively independent.
#'   }
#' @param testing_annotations Numeric matrix or data frame (N_test x M) of
#'   functional annotation scores for the variants being analyzed. Must have
#'   the same columns (annotation features) as the training data.
#'
#' @return A numeric vector of length N_test containing estimated
#'   variant-importance scores (PI values) for each variant in testing_annotations. Values
#'   are in the range (0, 1), where higher values indicate higher estimated
#'   importance based on the annotation profile.
#'
#' @details
#' \strong{Methodology:}
#'
#' This function implements annotation-based variant-importance estimation for the GLOW
#' framework. The key idea is to leverage functional annotations (e.g., CADD,
#' PolyPhen, conservation scores) to estimate which variants are more likely
#' to be important, allowing GLOW to upweight likely important variants in the
#' association test.
#'
#' \strong{Algorithm:}
#'
#' \enumerate{
#'   \item \strong{Ensemble Model Training}:
#'     \itemize{
#'       \item For i = 1 to model_need_N:
#'         \enumerate{
#'           \item Sample training_control_need_N controls randomly from
#'                 training_controlAnnotation
#'           \item Combine all cases with sampled controls
#'           \item Fit a logistic regression model (LASSO or GLM) predicting
#'                 case/control status (1 for cases, 0 for controls) from annotations
#'         }
#'       \item This creates model_need_N independent models
#'     }
#'   \item \strong{Prediction}:
#'     \itemize{
#'       \item Apply each of the model_need_N models to testing_annotations
#'       \item Average predictions across models to obtain final PI estimates
#'     }
#' }
#'
#' \strong{Why Ensemble Approach?}
#'
#' Using multiple models with different control samples:
#' \itemize{
#'   \item Reduces dependence on specific control sample
#'   \item Improves prediction stability and robustness
#'   \item Accounts for variability in control annotations
#'   \item Provides implicit uncertainty quantification
#' }
#'
#' \strong{LASSO vs GLM:}
#'
#' \emph{Use LASSO when:}
#' \itemize{
#'   \item Number of annotations M is large (e.g., M > 20)
#'   \item Annotations are correlated
#'   \item You want automatic feature selection
#'   \item Sample size is limited relative to M
#' }
#'
#' \emph{Use GLM when:}
#' \itemize{
#'   \item Number of annotations M is small (e.g., M < 10)
#'   \item Sample size is large relative to M
#'   \item All annotations are believed to be informative
#'   \item You want interpretable coefficients for all features
#' }
#'
#' \strong{Training Data Requirements:}
#'
#' Good training data should have:
#' \itemize{
#'   \item \strong{Cases}: Variants known or strongly suspected to be associated
#'     with the phenotype (e.g., from previous GWAS hits, Mendelian variants,
#'     functional studies)
#'   \item \strong{Controls}: Variants known or strongly suspected to be neutral
#'     (e.g., common synonymous variants, intergenic variants far from genes)
#'   \item \strong{Same trait}: Training and testing data should be for the same
#'     or similar phenotype
#'   \item \strong{Sufficient sample size}: At least 50-100 cases recommended;
#'     controls should be 2-10x more abundant
#' }
#'
#' \strong{Annotation Data Format:}
#'
#' Each annotation should be:
#' \itemize{
#'   \item Numeric (continuous or ordinal)
#'   \item Higher values should generally indicate higher pathogenicity
#'     (if opposite, multiply by -1 before input)
#'   \item Reasonably scaled (extreme ranges may require normalization)
#'   \item Available for all variants (no missing values)
#' }
#'
#' Common annotations include:
#' \itemize{
#'   \item CADD: Combined Annotation Dependent Depletion
#'   \item PolyPhen-2: Prediction of functional effects
#'   \item SIFT: Sorts Intolerant From Tolerant
#'   \item GERP: Genomic Evolutionary Rate Profiling
#'   \item PhyloP/PhastCons: Conservation scores
#'   \item Functional class indicators (missense, LoF, etc.)
#' }
#'
#' @section Computational Complexity:
#' O(model_need_N * cost_per_model + N_test * model_need_N * M) where:
#' \itemize{
#'   \item For LASSO: cost_per_model = O(N_train * M * n_lambda) with cross-validation
#'   \item For GLM: cost_per_model = O(N_train * M^2) for matrix inversion
#'   \item Prediction: O(N_test * model_need_N * M)
#' }
#'
#' Typical runtime: seconds to minutes depending on M and model_need_N.
#'
#' @section Assumptions:
#' \itemize{
#'   \item training_caseAnnotation and training_controlAnnotation have the same
#'     number and types of columns
#'   \item testing_annotations has the same columns as training data
#'   \item training_control_need_N <= nrow(training_controlAnnotation)
#'   \item No missing values in annotation data
#'   \item Annotations are informative for case/control status
#' }
#'
#' @section Output Interpretation:
#' \itemize{
#'   \item PI values close to 1: High predicted importance
#'     (based on annotations similar to known associated variants)
#'   \item PI values close to 0: Low predicted importance
#'     (based on annotations similar to known neutral variants)
#'   \item PI values near 0.5: Uncertain prediction (annotations don't strongly
#'     favor associated or neutral)
#' }
#'
#' The estimated PI values are then used as input to \code{\link{Optimal_Weights_M}}
#' to calculate optimal weights that incorporate biological prior knowledge.
#'
#' @section Cross-Validation:
#' For LASSO models, the optimal lambda parameter is selected using 10-fold
#' cross-validation (default in cv.glmnet) to minimize binomial deviance. Each
#' model in the ensemble performs its own cross-validation, which can lead to
#' different levels of sparsity across ensemble members.
#'
#' @examples
#' # Example 1: Basic usage with LASSO
#' set.seed(123)
#' # Simulate training data: cases have higher annotation scores
#' n_case <- 100
#' n_control <- 500
#' n_anno <- 5
#' training_case <- matrix(rnorm(n_case * n_anno, mean = 0.5, sd = 1),
#'                        ncol = n_anno)
#' training_control <- matrix(rnorm(n_control * n_anno, mean = 0, sd = 1),
#'                           ncol = n_anno)
#' colnames(training_case) <- colnames(training_control) <- paste0("Anno", 1:n_anno)
#'
#' # Simulate testing data
#' n_test <- 200
#' testing_anno <- matrix(rnorm(n_test * n_anno, mean = 0.2, sd = 1),
#'                       ncol = n_anno)
#' colnames(testing_anno) <- paste0("Anno", 1:n_anno)
#'
#' # Estimate PI using LASSO ensemble
#' PI_estimates <- get_PI(
#'   training_caseAnnotation = training_case,
#'   training_controlAnnotation = training_control,
#'   training_control_need_N = 100,
#'   model_need_N = 10,
#'   modelType = "LASSO",
#'   testing_annotations = testing_anno
#' )
#'
#' # Check output
#' summary(PI_estimates)
#' hist(PI_estimates, main = "Distribution of PI Estimates",
#'      xlab = "Estimated Variant Importance")
#'
#' \dontrun{
#' # Example 2: Using GLM instead of LASSO
#' PI_glm <- get_PI(
#'   training_caseAnnotation = training_case,
#'   training_controlAnnotation = training_control,
#'   training_control_need_N = 100,
#'   model_need_N = 10,
#'   modelType = "GLM",
#'   testing_annotations = testing_anno
#' )
#'
#' # Compare LASSO vs GLM predictions
#' plot(PI_estimates, PI_glm,
#'      xlab = "LASSO PI", ylab = "GLM PI",
#'      main = "LASSO vs GLM PI Estimates")
#' abline(0, 1, col = "red")
#' cor(PI_estimates, PI_glm)
#'
#' # Example 3: Using real annotation data (conceptual)
#' # Load variant annotations (CADD, PolyPhen, etc.)
#' # case_anno <- read.csv("causal_variants_annotations.csv")
#' # control_anno <- read.csv("neutral_variants_annotations.csv")
#' # test_anno <- read.csv("study_variants_annotations.csv")
#' #
#' # PI_real <- get_PI(
#' #   training_caseAnnotation = as.matrix(case_anno),
#' #   training_controlAnnotation = as.matrix(control_anno),
#' #   training_control_need_N = 200,
#' #   model_need_N = 20,
#' #   modelType = "LASSO",
#' #   testing_annotations = as.matrix(test_anno)
#' # )
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{model_PI}}: Generate ensemble of PI prediction models
#'   \item \code{LASSOmodel}: Fit single LASSO model (internal)
#'   \item \code{GLMmodel}: Fit single GLM model (internal)
#'   \item \code{\link{Optimal_Weights_M}}: Calculate optimal weights using PI estimates
#' }
#'
#' @references
#' Zhang, H., Liu, M., Landers, J. E., and Wu, Z. Integrated Weighted Association
#' Test with Application to Genetic Association Studies. Annals of Applied
#' Statistics (in revision).
#'
#' Ionita-Laza, I., Lee, S., Makarov, V., Buxbaum, J. D., and Lin, X. (2013).
#' Sequence kernel association tests for the combined effect of rare and common
#' variants. The American Journal of Human Genetics, 92(6), 841-853.
#'
#' Friedman, J., Hastie, T., and Tibshirani, R. (2010). Regularization paths
#' for generalized linear models via coordinate descent. Journal of Statistical
#' Software, 33(1), 1-22.
#'
#' @import stats
#' @import glmnet
#'
#' @export
get_PI <- function(training_caseAnnotation, training_controlAnnotation,
                   training_control_need_N, model_need_N, modelType,
                   testing_annotations) {

  # Input validation
  if (!is.matrix(training_caseAnnotation) && !is.data.frame(training_caseAnnotation)) {
    training_caseAnnotation <- as.matrix(training_caseAnnotation)
  }
  if (!is.matrix(training_controlAnnotation) && !is.data.frame(training_controlAnnotation)) {
    training_controlAnnotation <- as.matrix(training_controlAnnotation)
  }
  if (!is.matrix(testing_annotations) && !is.data.frame(testing_annotations)) {
    testing_annotations <- as.matrix(testing_annotations)
  }

  # Check dimensions
  if (ncol(training_caseAnnotation) != ncol(training_controlAnnotation)) {
    stop("training_caseAnnotation and training_controlAnnotation must have the same number of columns (annotations)")
  }
  if (ncol(testing_annotations) != ncol(training_caseAnnotation)) {
    stop("testing_annotations must have the same number of columns as training data")
  }

  # Check model type
  if (!modelType %in% c("LASSO", "GLM")) {
    stop("modelType must be either 'LASSO' or 'GLM'")
  }

  # Check training_control_need_N
  if (training_control_need_N > nrow(training_controlAnnotation)) {
    stop("training_control_need_N (", training_control_need_N,
         ") exceeds the number of control annotations (",
         nrow(training_controlAnnotation), ")")
  }

  # Check model_need_N
  if (model_need_N < 1) {
    stop("model_need_N must be at least 1")
  }

  # Generate ensemble of PI models
  # This function is already implemented in helpers_optimalWeights.R
  PImodels <- model_PI(
    caseAnnotation = training_caseAnnotation,
    controlAnnotation = training_controlAnnotation,
    modelType = modelType,
    control_need_N = training_control_need_N,
    model_need_N = model_need_N
  )

  # Preallocate matrix for predictions from each model
  # Each column represents predictions from one model
  piArr <- matrix(NA, nrow = nrow(testing_annotations), ncol = length(PImodels))

  # Generate PI predictions using each model in the ensemble
  for (i in 1:length(PImodels)) {
    model <- PImodels[[i]]

    if (modelType == "LASSO") {
      # For LASSO: use predict.glmnet with type = "response" to get probabilities
      # s = model$lambda uses the lambda value from model fitting
      piArr[, i] <- predict(
        model,
        newx = as.matrix(testing_annotations),
        type = "response",
        s = model$lambda
      )

    } else if (modelType == "GLM") {
      # For GLM: use predict.glm with type = "response" to get probabilities
      # This applies the logistic (sigmoid) function to convert log-odds to probabilities
      # Model was fit with y ~ x where x is a matrix in allData list
      # Must provide newdata as a list with 'x' as a matrix to match training
      piArr[, i] <- predict(
        model,
        newdata = list(x = as.matrix(testing_annotations)),
        type = "response"
      )
    }
  }

  # Average predictions across all models in the ensemble
  # This reduces variance and improves prediction stability
  PI <- rowMeans(piArr, na.rm = TRUE)

  return(PI)
}


#################### Helper Functions ####################

#' Get Training Data for PI Model
#'
#' @description
#' Prepares training data for the variant-importance score (PI) model by combining case
#' and control annotation data. This function samples control annotations to
#' balance the training dataset and formats the data for model fitting.
#'
#' @param caseAnnotation Numeric vector or matrix of annotation scores for cases
#' @param controlAnnotation Numeric vector or matrix of annotation scores for controls
#' @param control_need_N Integer, number of control samples to use in training
#'
#' @return A list with two elements:
#'   \item{x}{Matrix of annotation features (predictors)}
#'   \item{y}{Matrix of binary outcomes (1 for cases, 0 for controls)}
#'
#' @details
#' This function is used internally by \code{\link{model_PI}} to create balanced
#' training datasets. It combines all case data with a random sample of control
#' data to create a training set for logistic regression or LASSO models.
#'
#' The sampling is done without replacement to ensure diversity in the control
#' samples used for training.
#'
#' @section Computational Complexity:
#' O(n) where n is control_need_N (sampling and data frame operations)
#'
#' @section Assumptions:
#' \itemize{
#'   \item control_need_N <= nrow(controlAnnotation)
#'   \item caseAnnotation and controlAnnotation have compatible dimensions
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate annotation data
#' set.seed(123)
#' case_anno <- matrix(runif(50 * 5), ncol = 5)
#' control_anno <- matrix(runif(1000 * 5), ncol = 5)
#'
#' # Get training data
#' train_data <- getTrainData(case_anno, control_anno, control_need_N = 100)
#' str(train_data)
#' }
#'
#' @keywords internal
#' @noRd
getTrainData <- function(caseAnnotation, controlAnnotation, control_need_N) {
  # Input validation
  if (control_need_N > nrow(controlAnnotation)) {
    stop("control_need_N (", control_need_N, ") exceeds the number of total control annotations (",
         nrow(controlAnnotation), ")")
  }

  # Combine case data into a data frame
  # All cases are labeled as y = 1
  caseData <- data.frame(x = caseAnnotation, y = 1)

  # Randomly sample control annotations without replacement
  # This ensures diversity in the control samples
  ctrl_indices <- sample(nrow(controlAnnotation), control_need_N, replace = FALSE)
  ctrlAnno <- controlAnnotation[ctrl_indices, , drop = FALSE]
  ctrlData <- data.frame(x = ctrlAnno, y = 0)

  # Combine case and control data
  allData <- rbind(caseData, ctrlData)

  # Extract features (x) and outcomes (y) as matrices
  x <- as.matrix(allData[, names(allData) != "y", drop = FALSE])
  y <- as.matrix(allData[, "y", drop = FALSE])

  return(list(x = x, y = y))
}


#' Fit LASSO Model for PI Estimation
#'
#' @description
#' Fits a LASSO (L1-regularized) logistic regression model for variant-importance
#' score (PI) prediction. The optimal lambda parameter is selected via
#' cross-validation.
#'
#' @param allData A list with elements x (matrix of annotation features) and
#'   y (matrix of binary outcomes)
#'
#' @return A fitted glmnet model object with lambda set to the cross-validated
#'   minimum value
#'
#' @details
#' This function uses the glmnet package to fit a LASSO model for binary
#' classification. The LASSO penalty (alpha = 1) encourages sparsity in the
#' coefficient estimates, which can help identify the most important
#' annotation features for variant-importance prediction.
#'
#' The lambda parameter is selected using 10-fold cross-validation (default
#' in cv.glmnet) to minimize the binomial deviance.
#'
#' @section Computational Complexity:
#' O(n * p * k) where n is sample size, p is number of features, and k is
#' number of lambda values tested in cross-validation
#'
#' @section Dependencies:
#' Requires the glmnet package for LASSO model fitting
#'
#' @examples
#' \dontrun{
#' # Simulate annotation data
#' set.seed(123)
#' n <- 200
#' p <- 10
#' x <- matrix(rnorm(n * p), ncol = p)
#' y <- rbinom(n, 1, plogis(x[,1] + x[,2]))
#' allData <- list(x = x, y = y)
#'
#' # Fit LASSO model
#' model <- LASSOmodel(allData)
#' coef(model)
#' }
#'
#' @keywords internal
#' @noRd
LASSOmodel <- function(allData) {
  # Perform cross-validation to select optimal lambda
  # alpha = 1 specifies LASSO (L1 penalty)
  # family = "binomial" for binary outcome
  cv.lasso <- glmnet::cv.glmnet(
    x = allData$x,
    y = allData$y,
    alpha = 1,
    family = "binomial"
  )

  # Fit final model using the lambda that minimizes CV error
  model <- glmnet::glmnet(
    x = allData$x,
    y = allData$y,
    alpha = 1,
    family = "binomial",
    lambda = cv.lasso$lambda.min
  )

  return(model)
}


#' Fit GLM Model for PI Estimation
#'
#' @description
#' Fits a standard generalized linear model (GLM) with logistic regression for
#' variant-importance score (PI) prediction. This is an alternative to the LASSO
#' approach that does not impose sparsity constraints.
#'
#' @param allData A list with elements x (matrix of annotation features) and
#'   y (matrix of binary outcomes)
#'
#' @return A fitted glm object with binomial family
#'
#' @details
#' This function fits a standard logistic regression model using all provided
#' annotation features. Unlike LASSO, this approach does not perform feature
#' selection or regularization. It may be preferred when:
#' \itemize{
#'   \item The number of features is small relative to sample size
#'   \item All features are believed to be relevant
#'   \item Interpretability of all coefficients is important
#' }
#'
#' @section Computational Complexity:
#' O(n * p^2) where n is sample size and p is number of features (matrix inversion)
#'
#' @section Assumptions:
#' \itemize{
#'   \item Sample size >> number of features to avoid overfitting
#'   \item No perfect multicollinearity among features
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate annotation data
#' set.seed(123)
#' n <- 200
#' p <- 5
#' x <- matrix(rnorm(n * p), ncol = p)
#' y <- rbinom(n, 1, plogis(x[,1] + x[,2]))
#' allData <- list(x = x, y = y)
#'
#' # Fit GLM
#' model <- GLMmodel(allData)
#' summary(model)
#' }
#'
#' @keywords internal
#' @noRd
GLMmodel <- function(allData) {
  # Fit logistic regression model
  # family = "binomial" specifies logistic regression
  # When x is a matrix, glm automatically uses all columns as predictors
  model <- glm(y ~ x, data = allData, family = "binomial")

  return(model)
}


#' Generate Multiple PI Models
#'
#' @description
#' Generates multiple variant-importance score (PI) prediction models by repeatedly
#' sampling control data and fitting models. This ensemble approach can improve
#' prediction stability and account for variability in control sampling.
#'
#' @param caseAnnotation Numeric vector or matrix of annotation scores for cases
#' @param controlAnnotation Numeric matrix of annotation scores for controls
#' @param modelType Character string specifying model type: "LASSO" or "GLM"
#' @param control_need_N Integer, number of control samples to use per model
#' @param model_need_N Integer, number of models to generate
#'
#' @return A list of fitted model objects (length = model_need_N)
#'
#' @details
#' This function creates an ensemble of PI models by:
#' \enumerate{
#'   \item Repeatedly sampling different subsets of control data
#'   \item Fitting a model (LASSO or GLM) to each sampled dataset
#'   \item Returning all fitted models
#' }
#'
#' The ensemble approach helps to:
#' \itemize{
#'   \item Reduce dependence on specific control samples
#'   \item Improve prediction robustness
#'   \item Allow for uncertainty quantification via model averaging
#' }
#'
#' For LASSO models, each model uses cross-validation to select lambda
#' independently, which can lead to different levels of sparsity across models.
#'
#' @section Computational Complexity:
#' O(model_need_N * cost_per_model) where cost_per_model depends on modelType:
#' \itemize{
#'   \item LASSO: O(n * p * k) for CV with k lambda values
#'   \item GLM: O(n * p^2) for matrix inversion
#' }
#'
#' @section Assumptions:
#' \itemize{
#'   \item modelType is either "LASSO" or "GLM"
#'   \item control_need_N <= nrow(controlAnnotation)
#'   \item caseAnnotation and controlAnnotation have compatible dimensions
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate annotation data
#' set.seed(123)
#' case_anno <- matrix(runif(50 * 5), ncol = 5)
#' control_anno <- matrix(runif(1000 * 5), ncol = 5)
#'
#' # Generate 10 LASSO models
#' models <- model_PI(
#'   caseAnnotation = case_anno,
#'   controlAnnotation = control_anno,
#'   modelType = "LASSO",
#'   control_need_N = 100,
#'   model_need_N = 10
#' )
#'
#' length(models)  # Should be 10
#' }
#'
#' @export
model_PI <- function(caseAnnotation, controlAnnotation, modelType,
                     control_need_N, model_need_N) {
  # Input validation for NA/Inf
  # Convert data frames/lists to matrices for validation
  if (is.data.frame(caseAnnotation) || is.list(caseAnnotation)) {
    caseAnnotation <- as.matrix(caseAnnotation)
  }
  if (is.data.frame(controlAnnotation) || is.list(controlAnnotation)) {
    controlAnnotation <- as.matrix(controlAnnotation)
  }

  validate_numeric_input(caseAnnotation, "caseAnnotation")
  validate_numeric_input(controlAnnotation, "controlAnnotation")

  # Additional validation
  if (!modelType %in% c("LASSO", "GLM")) {
    stop("modelType must be either 'LASSO' or 'GLM'")
  }

  # Initialize models list
  models <- vector("list", model_need_N)

  # Generate models by repeatedly sampling and fitting
  for (i in seq_len(model_need_N)) {
    # Get training data with random control sample
    allData <- getTrainData(caseAnnotation, controlAnnotation, control_need_N)

    # Fit model based on specified type
    if (modelType == "LASSO") {
      models[[i]] <- LASSOmodel(allData)
    } else if (modelType == "GLM") {
      models[[i]] <- GLMmodel(allData)
    }
  }

  return(models)
}
