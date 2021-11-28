#' Pre-processing function for sex classification 
#' 
#' The purpose of this function is to process a single cell counts matrix into
#' the appropriate format for the \code{classifySex} function. 
#' 
#' This function will filter out cells that are unable to be classified due to 
#' zero counts on *XIST/Xist* and all of the Y chromosome genes. If 
#' \code{qc=TRUE} additional cells are removed as identified by the 
#' \code{perCellQCMetrics} and \code{quickPerCellQC} functions from the 
#' \code{scuttle} package. The resulting counts matrix is then log-normalised 
#' and scaled. 
#' 
#' @param x the counts matrix, rows are genes and columns are cells. Row names 
#' must be gene symbols. 
#' @param genome the genome the data arises from. Current options are 
#' human: genome = "Hs" or mouse: genome = "Mm".
#' @param qc logical, indicates whether to perform additional quality control on 
#' the cells. qc = TRUE will predict cells that pass quality control only and 
#' the filtered cells will not be classified. qc = FALSE will predict every cell 
#' except the cells with zero counts on *XIST/Xist* and the sum of the Y genes. 
#' Default is TRUE.
#' 
#' @return outputs a list object with the following components
#' \item{tcm.final }{A transposed count matrix where rows are cells and columns 
#' are the features used for classification.}
#' \item{data.df }{The normalised and scaled \code{tcm.final} matrix.} 
#' \item{discarded.cells }{Character vector of cell IDs for the cells that are
#' discarded when \code{qc=TRUE}.}
#' \item{zero.cells }{Character vector of cell IDs for the cells that can not
#' be classified as male/female due to zero counts on *Xist* and all the 
#' Y chromosome genes.}
#' 
#' @importFrom AnnotationDbi select
#' @importFrom stringr str_to_title 
#' @importFrom scuttle perCellQCMetrics
#' @importFrom scuttle quickPerCellQC
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @export preprocess
#'
#' @examples 
#' 
#' library(speckle)
#' library(SingleCellExperiment)
#' library(CellBench)
#' library(org.Hs.eg.db)
#'
#' # Get data from CellBench library
#' sc_data <- load_sc_data()
#' sc_10x <- sc_data$sc_10x
#'
#' # Get counts matrix in correct format with gene symbol as rownames 
#' # rather than ENSEMBL ID.
#' counts <- counts(sc_10x)
#' ann <- select(org.Hs.eg.db, keys=rownames(sc_10x),
#'              columns=c("ENSEMBL","SYMBOL"), keytype="ENSEMBL")
#' m <- match(rownames(counts), ann$ENSEMBL)
#' rownames(counts) <- ann$SYMBOL[m]
#' 
#' # Preprocess data
#' pro.data <- preprocess(counts, genome="Hs", qc = TRUE)
#' 
#' # Look at counts on XIST and superY.all
#' plot(pro.data$tcm.final$XIST, pro.data$tcm.final$superY)
#' 
#' # Cells that are identified as low quality
#' pro.data$discarded.cells
#' 
#' # Cells with zero counts on XIST and all Y genes
#' pro.data$zero.cells
#' 
preprocess<- function(x, genome=genome, qc=qc){
  
  x <- as.matrix(x)
  row.names(x)<- toupper(row.names(x))
  # genes located in the X chromosome that have been reported to escape 
  # X-inactivation
  # http://bioinf.wehi.edu.au/software/GenderGenes/index.html
  Xgenes<- c("ARHGAP4","STS","ARSD", "ARSL", "AVPR2", "BRS3", "S100G", "CHM",
             "CLCN4", "DDX3X","EIF1AX","EIF2S3", "GPM6B", "GRPR", "HCFC1",
             "L1CAM", "MAOA", "MYCLP1", "NAP1L3", "GPR143", "CDK16", "PLXNB3",
             "PRKX", "RBBP7", "RENBP", "RPS4X", "TRAPPC2", "SH3BGRL", "TBL1X",
             "UBA1", "KDM6A", "XG", "XIST", "ZFX", "PUDP", "PNPLA4", "USP9X",
             "KDM5C", "SMC1A", "NAA10", "OFD1", "IKBKG", "PIR", "INE2", "INE1",
             "AP1S2", "GYG2", "MED14", "RAB9A", "ITM2A", "MORF4L2", "CA5B",
             "SRPX2", "GEMIN8", "CTPS2", "CLTRN", "NLGN4X", "DUSP21", "ALG13",
             "SYAP1", "SYTL4", "FUNDC1", "GAB3", "RIBC1", "FAM9C","CA5BP1")
  
  # genes belonging to the male-specific region of chromosome Y (unique genes)
  # http://bioinf.wehi.edu.au/software/GenderGenes/index.html
  Ygenes<-c("AMELY", "DAZ1", "PRKY", "RBMY1A1", "RBMY1HP", "RPS4Y1", "SRY",
            "TSPY1", "UTY", "ZFY","KDM5D", "USP9Y", "DDX3Y", "PRY", "XKRY",
            "BPY2", "VCY", "CDY1", "EIF1AY", "TMSB4Y","CDY2A", "NLGN4Y",
            "PCDH11Y", "HSFY1", "TGIF2LY", "TBL1Y", "RPS4Y2", "HSFY2",
            "CDY2B", "TXLNGY","CDY1B", "DAZ3", "DAZ2", "DAZ4")
  
  # build artificial genes
  Xgene.set <-Xgenes[Xgenes %in% row.names(x)]
  Ygene.set <-Ygenes[Ygenes %in% row.names(x)]
  cm.new<-as.data.frame(matrix(rep(0, 3*ncol(x)), ncol = ncol(x),nrow = 3))
  row.names(cm.new) <- c("XIST","superX","superY")
  colnames(cm.new) <- colnames(x)
  cm.new["XIST", ]<- x["XIST", ]
  cm.new["superX", ] <-colSums(x[Xgene.set,])
  cm.new["superY", ] <-colSums(x[Ygene.set,])
  
#  if (genome == "Mm"){
#    ann <- suppressWarnings(AnnotationDbi::select(org.Mm.eg.db,keys=str_to_title(row.names(x)),
#                                  columns=c("SYMBOL","GENENAME","CHR"),
#                                  keytype="SYMBOL"))
#  }else{
#    ann <- suppressWarnings(AnnotationDbi::select(org.Hs.eg.db,keys=row.names(x),
#                                  columns=c("SYMBOL","GENENAME","CHR"),
#                                  keytype="SYMBOL"))
#  }
#  # create  superY.all
#  Ychr.genes<- toupper(unique(ann[which(ann$CHR=="Y"), "SYMBOL"]))
#  missing <- setdiff(Ychr.genes, row.names(x))
#  if (length(missing) >0){
#    Ychr.genes <- Ychr.genes[-match(missing, Ychr.genes)]
#  }
#  cm.new["superY.all", ] <- colSums(x[Ychr.genes,])
  
  ############################################################################
  # Pre-processing
  # perform simple QC
  # keep a copy of library size
  discarded.cells <- NA
  if (qc == TRUE){
    #data.sce <-SingleCellExperiment(assays = list(counts = x))
    qcstats <- scuttle::perCellQCMetrics(x,subsets=list(Mito=1:100))
    qcfilter <- scuttle::quickPerCellQC(qcstats, 
                                        percent_subsets=c("subsets_Mito_percent"))
    # save the discarded cells
    discarded.cells <- colnames(x[,qcfilter$discard])
    
    # cm.new only contains cells that pass the quality control
    cm.new <-cm.new[,!qcfilter$discard]
  }
  
  tcm.final <- t(cm.new)
  tcm.final <- as.data.frame(tcm.final)
  
  #Do Not Classify
  zero.cells <- NA
  dnc <- tcm.final$superY==0 & tcm.final$superX==0
  if(any(dnc)==TRUE){
      zero.cells <- row.names(tcm.final)[dnc]
      cat(length(zero.cells), "cell/s are unable to be classified due to an 
      abundance of zeroes on X and Y chromosome genes\n")
  }
  tcm.final <- tcm.final[!dnc, ]

  cm.new <- cm.new[,!dnc]
  cm.lib.size<- colSums(x[,colnames(cm.new)], na.rm=TRUE)
  
  # log-normalisation performed for each cell
  # scaling performed for each gene
  normsca.cm <- normSca(cm.new, lib.size=cm.lib.size)
  data.df <- t(normsca.cm)
  data.df <- as.data.frame(data.df)
  
  list(tcm.final=tcm.final, data.df=data.df, discarded.cells=discarded.cells, 
       zero.cells=zero.cells)
}

