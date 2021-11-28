#' Predict sex of cells in scRNA-seq data
#' 
#' This function will predict the sex for each cell in scRNA-seq data. The 
#' classifier is based on logistic regression models that have been trained
#' on mouse and human single cell RNA-seq data.
#'
#' For bulk RNA-seq, checking the sex of the samples for mouse and human 
#' experiments is trivial as we can simply check the expression of *Xist/XIST*. 
#' It is not as simple for single cell RNA-seq data as the number of counts 
#' measured per gene and per cell is often quite low. Simply relying on cut-offs
#' on the expression of genes like *Xist* means that many cells are unable to 
#' be classified. Hence we have developed a classifier based on a combination of
#' X- and Y-linked genes in order to accurately predict the sex of each cell.
#' 
#' Cells with zero counts on Xist and the sum of the Y chromosome genes will 
#' not be classified as there is simply not enough information to accurately 
#' classify as Male/Female, and NAs will be returned. In addition, the user has 
#' the option to perform quality control on the data first, by specifying 
#' \code{qc=TRUE}, which will not classify cells that are deemed low-quality.
#' 
#' @aliases classifySex
#' @param x counts matrix, rows correspond to genes and columns correspond to 
#' cells. Row names must be gene symbols. 
#' @param genome the genome the data arises from. Current options are 
#' human: genome = "Hs" or mouse: genome = "Mm".
#' @param qc logical, indicates whether to perform quality control or not. 
#' qc = TRUE will predict cells that pass quality control only and the filtered 
#' cells will not be classified. qc = FALSE will predict every cell except the 
#' cells with zero counts on *XIST/Xist* and the sum of the Y genes. Default is TRUE.
#' 
#' @return a dataframe with predicted labels for each cell
#' 
#' @importFrom stats predict  
#' @export classifySex
#' 
#' @author Xinyi Jin
#' 
#' @examples 
#' 
#' library(speckle)
#' library(SingleCellExperiment)
#' library(CellBench)
#' library(org.Hs.eg.db)
#'
#' sc_data <- load_sc_data()
#' sc_10x <- sc_data$sc_10x
#'
#' counts <- counts(sc_10x)
#' ann <- select(org.Hs.eg.db, keys=rownames(sc_10x),
#'              columns=c("ENSEMBL","SYMBOL"), keytype="ENSEMBL")
#' m <- match(rownames(counts), ann$ENSEMBL)
#' rownames(counts) <- ann$SYMBOL[m]
#'
#' sex <- classifySex(counts, genome="Hs")
#' 
#' table(sex$prediction)
#' boxplot(counts["XIST",]~sex$prediction)
#'    
classifySex<-function(x, genome=NULL, qc = TRUE)
#    Classify cells as male or female
#    Xinyi Jin and Belinda Phipson
#    11 February 2021
#    Modified 11 February 2021
{
    # Perform some checks on the data
    if(is.null(x)) stop("Counts matrix missing")
    x <- as.matrix(x)
    if(is.null(genome)){
        message("Genome not specified. Human genome used. Options are 'Hs' for 
        human and 'Mm' for mouse. We currently don't support other genomes.")
    }
    # Default is Hs
    genome <- match.arg(genome,c("Hs","Mm"))
    
    # pre-process 
    processed.data<-preprocess(x, genome = genome, qc = qc)
  
    # the processed transposed count matrix 
    tcm <-processed.data$tcm.final
  
    # the normalised, scaled transposed count matrix 
    data.df <- processed.data$data.df
  
    # cells that filtered by QC
    discarded.cells <- processed.data$discarded.cells
  
    # cells with zero count on XIST and superY.all
    zero.cells <- processed.data$zero.cells
  
    # store the final predictions 
    final.pred<-data.frame(prediction=rep("NA", ncol(x)))
    row.names(final.pred)<- colnames(x)
  
    # load trained models 
    if(genome == "Mm"){
      model <- Mm.model
    }
    else{
      model <- Hs.model
    }

    preds <- predict(model, newdata = data.df)
    final.pred[row.names(data.df), "prediction"]<- as.character(preds)
  
    final.pred
}