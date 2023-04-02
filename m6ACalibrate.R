#' Calibrate putative m6A maps obtained from peak-calling.
#'
#' @description This function takes in putative m6A maps as input and generates calibrated and high-confidence m6A maps as output.
#'
#' @param x A \code{GRanges} object for the genomic location of m6A sites.
#' @param txdb A \code{\link{TxDb}} object for the transcript annotation.
#'
#' The \code{TxDb} can be obtained from either Bioconductor or from the GFF/GTF files using the function \code{\link{makeTxDbFromGFF}}.
#'
#' @param bsgenome A \code{character} or a \code{\link{BSgenome}} for the reference genome.
#'
#' The character should be the UCSC genome name which is acceptable by \code{\link{getBSgenome}} or/and \code{\link{makeTxDbFromUCSC}}; example: \code{"hg19"}.
#'
#' @param type An \code{character} specifying the type of m6A-specific antibody, can be one of \code{c("Abcam", "NEB", "SYSY", "ensemble")}; default \code{= "ensemble"}.
#'
#' @importFrom randomForest randomForest
#' @import GenomicFeatures
#' @import GenomicRanges
#' @import BSgenome
#' @import IRanges
#' @import AnnotationDbi
#'
#' @return A \code{GRanges} object containing the calibrated m6A sites.
#'
#' @examples
#' library(m6ACalibrateR)
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' # Load example data
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38
#'
#' x <- readRDS(system.file("extdata", "peaks.rds", package = "m6ACalibrateR"))
#'
#' # Calibrate m6A maps with default parameters
#' calibrated_m6A <- m6ACalibrate(x, txdb, bsgenome, "ensemble")
#' calibrated_m6A
#'
#' @export

m6ACalibrate <- function(x, txdb, bsgenome, type = c("Abcam", "NEB", "SYSY", "ensemble")) {

  # Check input validity
  if(!is(x, "GRanges")) {
    stop("Argument 'x' must be a GRanges object.")
  }

  if(!is(txdb, "TxDb")) {
    stop("Argument 'txdb' must be a TxDb object.")
  }

  if(!(is.character(bsgenome)|is(bsgenome, "BSgenome")|is.null(bsgenome))) {
    stop("Argument 'bsgenome' must be a BSgenome object.")
  }

  if(!is.character(type) | !type %in% c("Abcam", "NEB", "SYSY", "ensemble")) {
    stop("Argument 'type' must be a character vector and must be one of 'Abcam', 'NEB', 'SYSY', or 'ensemble'.")
  }

  # Prepare reference genome
  if (is.character(bsgenome)){
    bsgenome <- getBSgenome(bsgenome)
  }

  # Load random forest models
  rf_Abcam <- readRDS(system.file("data", "Abcam.rds", package = "m6ACalibrateR"))
  rf_NEB <- readRDS(system.file("data", "NEB.rds", package = "m6ACalibrateR"))
  rf_SYSY <- readRDS(system.file("data", "SYSY.rds", package = "m6ACalibrateR"))

  # Extract all DRACH motif on exons
  exbtx <- exonsBy(txdb, by = "tx")
  motif_all <- sort(sampleSequence("DRACH", exbtx, bsgenome) - 2)

  # Find DRACH motifs overlapping with input m6A sites
  x_motif <- subsetByOverlaps(motif_all, x)

  # Extract genomic features
  gfeatures <- feature_extraction(x_motif, txdb)

  # Use random forest models to predict m6A sites
  if(type == "Abcam") {
    rf_pred <- predict(rf_Abcam, newdata = gfeatures)
    filter_index <- which(rf_pred == 0)
  } else if(type == "NEB") {
    rf_pred <- predict(rf_NEB, newdata = gfeatures)
    filter_index <- which(rf_pred == 0)
  } else if(type == "SYSY") {
    rf_pred <- predict(rf_SYSY, newdata = gfeatures)
    filter_index <- which(rf_pred == 0)
  } else {
    rf_Abcam_pred <- predict(rf_Abcam, newdata = gfeatures, type = "prob")
    rf_NEB_pred <- predict(rf_NEB, newdata = gfeatures, type = "prob")
    rf_SYSY_pred <- predict(rf_SYSY, newdata = gfeatures, type = "prob")
    predp0 = data.frame(abcam = rf_Abcam_pred[,1], neb = rf_NEB_pred[,1], sysy = rf_SYSY_pred[,1])
    filter_index <- which(rowMeans(predp0) >= 0.5)
  }

  # Filter out low-confidence m6A sites
  motif_x <- subsetByOverlaps(x, motif_all)
  filter_x_motif <- x_motif[filter_index]
  filter_x_index <- unique(queryHits(findOverlaps(motif_x, filter_x_motif)))
  x_retain <- motif_x[-filter_x_index]

  # Check if any m6A sites were filtered out
  num_filtered <- length(filter_x_index)
  if(num_filtered > 0) {
    message(paste0(num_filtered, " m6A sites were filtered out due to low-confidence predictions."))
  }

  return(x_retain)
}
