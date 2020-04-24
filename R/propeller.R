#' Finding statistically significant differences in cell type proportions
#'
#' Calculates cell type proportions, performs a variance stabilising
#' transformation on the proportions and determines whether the cell type
#' proportions are statistically significant between different groups using
#' linear modelling.
#'
#' This function will take a \code{SingleCellExperiment} or \code{Seurat}
#' object and extract the \code{group}, \code{sample} and \code{clusters} cell
#' information. The user can either state these factor vectors explicitly in
#' the call to the \code{propeller} function, or internal functions will
#' extract them from the relevants objects. The user must ensure that
#' \code{group} and \code{sample} are columns in the metadata assays of the
#' relevant objects (any combination of upper/lower case is acceptable). For
#' \code{Seurat} objects the clusters are extracted using the \code{Idents}
#' function. For \code{SingleCellExperiment} objects, \code{clusters} needs to
#' be a column in the \code{colData} assay.
#'
#' The \code{propeller} function calculates cell type proportions for each
#' biological replicate, performs a variance stabilising transformation on the
#' matrix of proportions and fits a linear model for each cell type or cluster
#' using the \code{limma} framework. Propeller tests whether there is a
#' difference in the cell type proportions between multiple groups. If there
#' are only 2 groups, a t-test is used to calculate p-values, and if there are
#' more than 2 groups, an F-test (ANOVA) is used. Cell type proportions of 1 or
#' 0 are accommodated. Benjamini and Hochberg false discovery rates are
#' calculated to account to multiple testing of cell types/clusters.
#'
#' @aliases propeller
#' @param x object of class \code{SingleCellExperiment} or \code{Seurat}
#' @param clusters a factor specifying the cluster or cell type for every cell.
#' For \code{SingleCellExperiment} objects this should correspond to a column
#' called \code{clusters} in the \code{colData} assay. For \code{Seurat}
#' objects this will be extracted by a call to \code{Idents(x)}.
#' @param sample a factor specifying the biological replicate for each cell.
#' For \code{SingleCellExperiment} objects this should correspond to a column
#' called \code{sample} in the \code{colData} assay and for \code{Seurat}
#' objects this should correspond to \code{x$sample}.
#' @param group a factor specifying the groups of interest for performing the
#' differential proportions analysis. For \code{SingleCellExperiment} objects
#' this should correspond to a column called \code{group} in the \code{colData}
#' assay.  For \code{Seurat} objects this should correspond to \code{x$group}.
#' @param trend logical, if true fits a mean variance trend on the transformed
#' proportions
#' @param robust logical, if true performs robust empirical Bayes shrinkage of
#' the variances
#'
#' @return produces a dataframe of results
#' @export
#'
#' @author Belinda Phipson
#'
#' @seealso \code{\link{lmFit}}, \code{\link{eBayes}},
#' \code{\link{getTransformedProps}}
#'
#' @references Smyth, G.K. (2004). Linear models and empirical Bayes methods
#' for assessing differential expression in microarray experiments.
#' \emph{Statistical Applications in Genetics and Molecular Biology}, Volume
#' \bold{3}, Article 3.
#'
#' Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery
#' rate: a practical and powerful approach to multiple testing. \emph{Journal
#' of the Royal Statistical Society Series}, B, \bold{57}, 289-300.
#'
#' @examples
#'
#'   library(speckle)
#'   library(ggplot2)
#'   library(limma)
#'
#'   # Make up some data
#'   # True cell type proportions for 4 samples
#'   p_s1 <- c(0.5,0.3,0.2)
#'   p_s2 <- c(0.6,0.3,0.1)
#'   p_s3 <- c(0.3,0.4,0.3)
#'   p_s4 <- c(0.4,0.3,0.3)
#'
#'   # Total numbers of cells per sample
#'   numcells <- c(1000,1500,900,1200)
#'
#'   # Generate cell-level vector for sample info
#'   biorep <- rep(c("s1","s2","s3","s4"),numcells)
#'   length(biorep)
#'
#'   # Numbers of cells for each of the 3 clusters per sample
#'   n_s1 <- p_s1*numcells[1]
#'   n_s2 <- p_s2*numcells[2]
#'   n_s3 <- p_s3*numcells[3]
#'   n_s4 <- p_s4*numcells[4]
#'
#'   # Assign cluster labels for 4 samples
#'   cl_s1 <- rep(c("c0","c1","c2"),n_s1)
#'   cl_s2 <- rep(c("c0","c1","c2"),n_s2)
#'   cl_s3 <- rep(c("c0","c1","c2"),n_s3)
#'   cl_s4 <- rep(c("c0","c1","c2"),n_s4)
#'
#'   # Generate cell-level vector for cluster info
#'   clust <- c(cl_s1,cl_s2,cl_s3,cl_s4)
#'   length(clust)
#'
#'   # Assume s1 and s2 belong to group 1 and s3 and s4 belong to group 2
#'   grp <- rep(c("grp1","grp2"),c(sum(numcells[1:2]),sum(numcells[3:4])))
#'
#'   propeller(clusters = clust, sample = biorep, group = grp,
#'   robust = FALSE, trend = FALSE)
#'
propeller <- function(x=NULL, clusters=NULL, sample=NULL, group=NULL,
                      trend=FALSE, robust=TRUE)
#    Testing for differences in cell type proportions
#    Belinda Phipson
#    29 July 2019
#    Modified 22 April 2020
{

    if(is.null(x) & is.null(sample) & is.null(group) & is.null(clusters))
        stop("Please provide either a SingleCellExperiment object or Seurat
        object with required annotation metadata, or explicitly provide
        clusters, sample and group information")

    if((is.null(clusters) | is.null(sample) | is.null(group)) & !is.null(x)){
        # Extract cluster, sample and group info from SCE object
        if(is(x,"SingleCellExperiment"))
            y <- .extractSCE(x)

        # Extract cluster, sample and group info from Seurat object
        if(is(x,"Seurat"))
            y <- .extractSeurat(x)

        clusters <- y$clusters
        sample <- y$sample
        group <- y$group
    }

    # Get transformed proportions
    prop.list <- getTransformedProps(clusters, sample)

    # Calculate baseline proportions for each cluster
    baseline.props <- table(clusters)/sum(table(clusters))

    # Collapse group information
    group.coll <- table(sample, group)

    design <- matrix(as.integer(group.coll != 0), ncol=ncol(group.coll))
    colnames(design) <- colnames(group.coll)

    if(ncol(design)==2){
        message("group variable has 2 levels, t-tests will be performed")
        contrasts <- c(1,-1)
        out <- propeller.ttest(prop.list, design, contrasts=contrasts,
                               robust=robust, trend=trend)
        out <- data.frame(BaselineProp=baseline.props,out)
        out
    }
    else if(ncol(design)>=2){
        message("group variable has > 2 levels, ANOVA will be performed")
        coef <- seq_len(ncol(design))
        out <- propeller.anova(prop.list, design, coef=coef, robust=robust,
                               trend=trend)
        out <- data.frame(BaselineProp=as.vector(baseline.props),out)
        out
  }

}

#' Extract metadata from \code{SingleCellExperiment} object
#'
#' This is an accessor function that extracts cluster, sample and group
#' information for each cell.
#'
#' @param x object of class \code{SingleCellExperiment}
#'
#' @return a dataframe containing clusters, sample and group
#' @export
#' @author Belinda Phipson
#'
.extractSCE <- function(x){
    message("extracting sample information from SingleCellExperiment object")
    colnames(colData(x)) <- toupper(colnames(colData(x)))
    clusters <- factor(colData(x)$CLUSTER)
    sample <- factor(colData(x)$SAMPLE)
    group <- factor(colData(x)$GROUP)
    data.frame(clusters=clusters,sample=sample,group=group)
}

#' Extract metadata from \code{Seurat} object
#'
#' This is an accessor function that extracts cluster, sample and group
#' information for each cell.
#'
#' @param x object of class \code{Seurat}
#'
#' @return a dataframe containing clusters, sample and group
#' @export
#' @author Belinda Phipson
#'
.extractSeurat <- function(x){
    message("extracting sample information from Seurat object")
    colnames(x@meta.data) <- toupper(colnames(x@meta.data))
    clusters <- factor(Idents(x))
    sample <- factor(x$SAMPLE)
    group <- factor(x$GROUP)
    data.frame(clusters=clusters,sample=sample,group=group)
}


#' Calculates and transforms cell type proportions
#'
#' Calculates cell types proportions based on clusters/cell types and sample
#' information and performs a variance stabilising transformation on the
#' proportions.
#'
#' This function is called by the \code{propeller} function and calculates cell
#' type proportions and performs an arcsin-square root transformation.
#'
#' @param clusters a factor specifying the cluster or cell type for every cell.
#' @param sample a factor specifying the biological replicate for every cell.
#' @return outputs a list object with the following components
#' \item{TransformedProps }{A matrix of transformed cell type proportions with
#' the rows corresponding to the clusters/cell types and the columns
#' corresponding to the biological replicates/samples.} \item{Proportions }{A
#' matrix of cell type proportions with the rows corresponding to the
#' clusters/cell types and the columns corresponding to the biological
#' replicates/samples.}
#'
#' @export
#'
#' @author Belinda Phipson
#'
#' @seealso \code{\link{propeller}}
#'
#' @examples
#'
#'   library(speckle)
#'   library(ggplot2)
#'   library(limma)
#'
#'   # Make up some data
#'
#'   # True cell type proportions for 4 samples
#'   p_s1 <- c(0.5,0.3,0.2)
#'   p_s2 <- c(0.6,0.3,0.1)
#'   p_s3 <- c(0.3,0.4,0.3)
#'   p_s4 <- c(0.4,0.3,0.3)
#'
#'   # Total numbers of cells per sample
#'   numcells <- c(1000,1500,900,1200)
#'
#'   # Generate cell-level vector for sample info
#'   biorep <- rep(c("s1","s2","s3","s4"),numcells)
#'   length(biorep)
#'
#'   # Numbers of cells for each of 3 clusters per sample
#'   n_s1 <- p_s1*numcells[1]
#'   n_s2 <- p_s2*numcells[2]
#'   n_s3 <- p_s3*numcells[3]
#'   n_s4 <- p_s4*numcells[4]
#'
#'   cl_s1 <- rep(c("c0","c1","c2"),n_s1)
#'   cl_s2 <- rep(c("c0","c1","c2"),n_s2)
#'   cl_s3 <- rep(c("c0","c1","c2"),n_s3)
#'   cl_s4 <- rep(c("c0","c1","c2"),n_s4)
#'
#'   # Generate cell-level vector for cluster info
#'   clust <- c(cl_s1,cl_s2,cl_s3,cl_s4)
#'   length(clust)
#'
#'   getTransformedProps(clusters = clust, sample = biorep)
#'
getTransformedProps <- function(clusters=clusters, sample=sample)
{
    tab <- table(sample, clusters)
    props <- tab/rowSums(tab)
    prop.trans <- asin(sqrt(props))
    list(TransformedProps=t(prop.trans), Proportions=t(props))
}

#' Performs t-tests of transformed cell type proportions
#'
#' This function is called by \code{propeller} and performs t-tests between two
#' experimental groups or conditions on the transformed cell type proportions.
#'
#' In order to run this function, the user needs to run the
#' \code{getTransformedProps} function first. The output from
#' \code{getTransformedProps} is used as input. The \code{propeller.ttest}
#' function expects that the design matrix is not in the intercept format
#' and a contrast vector needs to be supplied. This contrast vector will
#' identify the two groups to be tested. Note that additional confounding
#' covariates can be accounted for as extra columns in the design matrix after
#' the group-specific columns.
#'
#' The \code{propeller.ttest} function uses the \code{limma} functions
#' \code{lmFit}, \code{contrasts.fit} and \code{eBayes} which has the additional
#' advantage that empirical Bayes shrinkage of the variances are performed.
#'
#' @param prop.list a list object containing two matrices:
#' \code{TransformedProps} and \code{Proportions}
#' @param design a design matrix with rows corresponding to samples and columns
#' to coefficients to be estimated
#' @param contrasts a vector specifying which columns of the design matrix
#' correspond to the two groups to test
#' @param robust logical, should robust variance estimation be used. Defaults to
#' TRUE.
#' @param trend logical, should a trend between means and variances be accounted
#' for. Defaults to FALSE.
#'
#' @return produces a dataframe of results
#' @export
#'
#' @author Belinda Phipson
#'
#' @seealso \code{\link{propeller}}, \code{\link{getTransformedProps}},
#' \code{\link{lmFit}}, \code{\link{contrasts.fit}}, \code{\link{eBayes}}
#'
#' @examples
#'   library(speckle)
#'   library(ggplot2)
#'   library(limma)
#'
#'   # Make up some data
#'
#'   # True cell type proportions for 4 samples
#'   p_s1 <- c(0.5,0.3,0.2)
#'   p_s2 <- c(0.6,0.3,0.1)
#'   p_s3 <- c(0.3,0.4,0.3)
#'   p_s4 <- c(0.4,0.3,0.3)
#'
#'   # Total numbers of cells per sample
#'   numcells <- c(1000,1500,900,1200)
#'
#'   # Generate cell-level vector for sample info
#'   biorep <- rep(c("s1","s2","s3","s4"),numcells)
#'   length(biorep)
#'
#'   # Numbers of cells for each of 3 clusters per sample
#'   n_s1 <- p_s1*numcells[1]
#'   n_s2 <- p_s2*numcells[2]
#'   n_s3 <- p_s3*numcells[3]
#'   n_s4 <- p_s4*numcells[4]
#'
#'   cl_s1 <- rep(c("c0","c1","c2"),n_s1)
#'   cl_s2 <- rep(c("c0","c1","c2"),n_s2)
#'   cl_s3 <- rep(c("c0","c1","c2"),n_s3)
#'   cl_s4 <- rep(c("c0","c1","c2"),n_s4)
#'
#'   # Generate cell-level vector for cluster info
#'   clust <- c(cl_s1,cl_s2,cl_s3,cl_s4)
#'   length(clust)
#'
#'   prop.list <- getTransformedProps(clusters = clust, sample = biorep)
#'
#'   # Assume s1 and s2 belong to group 1 and s3 and s4 belong to group 2
#'   grp <- rep(c("A","B"), each=2)
#'
#'   design <- model.matrix(~0+grp)
#'   design
#'
#'   # Compare Grp A to B
#'   contrasts <- c(1,-1)
#'
#'   propeller.ttest(prop.list, design=design, contrasts=contrasts, robust=TRUE,
#'   trend=FALSE)
#'
propeller.ttest <- function(prop.list=prop.list, design=design,
                            contrasts=contrasts, robust=robust, trend=trend)
{
    prop.trans <- prop.list$TransformedProps
    prop <- prop.list$Proportions

    fit <- lmFit(prop.trans, design)
    fit.cont <- contrasts.fit(fit, contrasts=contrasts)
    fit.cont <- eBayes(fit.cont, robust=robust, trend=trend)

    # Get mean cell type proportions and relative risk for output
    fit.prop <- lmFit(prop, design)
    z <- apply(fit.prop$coefficients, 1, function(x) x^contrasts)
    RR <- apply(z, 2, prod)

    fdr <- p.adjust(fit.cont$p.value, method="BH")

    out <- data.frame(PropMean=fit.prop$coefficients, PropRatio=RR,
                      Tstatistic=fit.cont$t[,1], P.Value=fit.cont$p.value[,1],
                      FDR=fdr)
    o <- order(out$P.Value)
    out[o,]
}

#' Performs F-tests for transformed cell type proportions
#'
#' This function is called by \code{propeller} and performs F-tests between
#' multiple experimental groups or conditions (> 2) on arcsin square root
#' transformed cell type proportions.
#'
#' In order to run this function, the user needs to run the
#' \code{getTransformedProps} function first. The output from
#' \code{getTransformedProps} is used as input. The \code{propeller.anova}
#' function expects that the design matrix is not in the intercept format.
#' This \code{coef} vector will identify the columns in the design matrix that
#' correspond to the groups being tested.
#' Note that additional confounding covariates can be accounted for as extra
#' columns in the design matrix, but need to come after the group-specific
#' columns.
#'
#' The \code{propeller.anova} function uses the \code{limma} functions
#' \code{lmFit} and \code{eBayes} to extract F statistics and p-values.
#' This has the additional advantage that empirical Bayes shrinkage of the
#' variances are performed.
#'
#' @param prop.list a list object containing two matrices:
#' \code{TransformedProps} and \code{Proportions}
#' @param design a design matrix with rows corresponding to samples and columns
#' to coefficients to be estimated
#' @param coef a vector specifying which the columns of the design matrix
#' corresponding to the groups to test
#' @param robust logical, should robust variance estimation be used. Defaults to
#' TRUE.
#' @param trend logical, should a trend between means and variances be accounted
#' for. Defaults to FALSE.
#'
#' @return produces a dataframe of results
#'
#' @export
#'
#' @author Belinda Phipson
#'
#'@seealso \code{\link{propeller}}, \code{\link{getTransformedProps}},
#' \code{\link{lmFit}}, \code{\link{eBayes}}
#'
#' @examples
#'   library(speckle)
#'   library(ggplot2)
#'   library(limma)
#'
#'   # Make up some data
#'
#'   # True cell type proportions for 4 samples
#'   p_s1 <- c(0.5,0.3,0.2)
#'   p_s2 <- c(0.6,0.3,0.1)
#'   p_s3 <- c(0.3,0.4,0.3)
#'   p_s4 <- c(0.4,0.3,0.3)
#'   p_s5 <- c(0.8,0.1,0.1)
#'   p_s6 <- c(0.75,0.2,0.05)
#'
#'   # Total numbers of cells per sample
#'   numcells <- c(1000,1500,900,1200,1000,800)
#'
#'   # Generate cell-level vector for sample info
#'   biorep <- rep(c("s1","s2","s3","s4","s5","s6"),numcells)
#'   length(biorep)
#'
#'   # Numbers of cells for each of 3 clusters per sample
#'   n_s1 <- p_s1*numcells[1]
#'   n_s2 <- p_s2*numcells[2]
#'   n_s3 <- p_s3*numcells[3]
#'   n_s4 <- p_s4*numcells[4]
#'   n_s5 <- p_s5*numcells[5]
#'   n_s6 <- p_s6*numcells[6]
#'
#'   cl_s1 <- rep(c("c0","c1","c2"),n_s1)
#'   cl_s2 <- rep(c("c0","c1","c2"),n_s2)
#'   cl_s3 <- rep(c("c0","c1","c2"),n_s3)
#'   cl_s4 <- rep(c("c0","c1","c2"),n_s4)
#'   cl_s5 <- rep(c("c0","c1","c2"),n_s5)
#'   cl_s6 <- rep(c("c0","c1","c2"),n_s6)
#'
#'   # Generate cell-level vector for cluster info
#'   clust <- c(cl_s1,cl_s2,cl_s3,cl_s4,cl_s5,cl_s6)
#'   length(clust)
#'
#'   prop.list <- getTransformedProps(clusters = clust, sample = biorep)
#'
#'   # Assume s1 and s2 belong to group A, s3 and s4 belong to group B, s5 and
#'   # s6 belong to group C
#'   grp <- rep(c("A","B","C"), each=2)
#'
#'   # Make sure design matrix does not have an intercept term
#'   design <- model.matrix(~0+grp)
#'   design
#'
#'   propeller.anova(prop.list, design=design, coef=c(1,2,3), robust=TRUE,
#'   trend=FALSE)
#'
propeller.anova <- function(prop.list=prop.list, design=design, coef = coef,
                            robust=robust, trend=trend)
{
    prop.trans <- prop.list$TransformedProps
    prop <- prop.list$Proportions

    # get cell type mean proportions ignoring other variables
    # this assumes that the design matrix is not in Intercept format
    fit.prop <- lmFit(prop, design[,coef])

    # Change design matrix to intercept format
    design[,1] <- 1
    colnames(design)[1] <- "Int"

    # Fit linear model taking into account all confounding variables
    fit <- lmFit(prop.trans,design)

    # Get F statistics corresponding to group information only
    # You have to remove the intercept term for this to work
    fit <- eBayes(fit[,coef[-1]], robust=robust, trend=trend)

    # Extract F p-value
    p.value <- fit$F.p.value
    # and perform FDR adjustment
    fdr <- p.adjust(fit$F.p.value, method="BH")

    out <- data.frame(PropMean=fit.prop$coefficients, Fstatistic= fit$F,
                      P.Value=p.value, FDR=fdr)
    o <- order(out$P.Value)
    out[o,]
}


#' Plot cell type proportions for each sample
#'
#' This is a plotting function that shows the cell type composition for each
#' sample as a stacked barplot. The \code{plotCellTypeProps} returns a
#' \code{ggplot2} object enabling the user to make style changes as required.
#'
#' @param x object of class \code{SingleCellExperiment} or \code{Seurat}
#' @param clusters a factor specifying the cluster or cell type for every cell.
#' For \code{SingleCellExperiment} objects this should correspond to a column
#' called \code{clusters} in the \code{colData} assay. For \code{Seurat}
#' objects this will be extracted by a call to \code{Idents(x)}.
#' @param sample a factor specifying the biological replicate for each cell.
#' For \code{SingleCellExperiment} objects this should correspond to a column
#' called \code{sample} in the \code{colData} assay and for \code{Seurat}
#' objects this should correspond to \code{x$sample}.
#'
#' @return a ggplot2 object
#' @export
#'
#' @author Belinda Phipson
#'
#' @examples
#'
#' library(speckle)
#' library(ggplot2)
#' library(limma)
#'
#' # Generate some fake data from a multinomial distribution
#' # Group A, 4 samples, 1000 cells in each sample
#' countsA <- rmultinom(4, size=1000, prob=c(0.1,0.3,0.6))
#' colnames(countsA) <- paste("s",1:4,sep="")
#'
#' # Group B, 3 samples, 800 cells in each sample
#'
#' countsB <- rmultinom(3, size=800, prob=c(0.2,0.05,0.75))
#' colnames(countsB) <- paste("s",5:7,sep="")
#' rownames(countsA) <- rownames(countsB) <- paste("c",0:2,sep="")
#'
#' allcounts <- cbind(countsA, countsB)
#' sample <- c(rep(colnames(allcounts),allcounts[1,]),
#'           rep(colnames(allcounts),allcounts[2,]),
#'           rep(colnames(allcounts),allcounts[3,]))
#' clust <- rep(rownames(allcounts),rowSums(allcounts))
#'
#' plotCellTypeProps(clusters=clust, sample=sample)
#'
plotCellTypeProps <- function(x=NULL, clusters=NULL, sample=NULL)
{
    if(is.null(x) & is.null(sample) & is.null(clusters))
        stop("Please provide either a SingleCellExperiment object or Seurat
            object with required annotation metadata, or explicitly provide
            clusters and sample information")

    if((is.null(clusters) | is.null(sample)) & !is.null(x)){
        # Extract cluster, sample and group info from SCE object
        if(is(x,"SingleCellExperiment"))
            y <- .extractSCE(x)

        # Extract cluster, sample and group info from Seurat object
        if(is(x,"Seurat"))
            y <- .extractSeurat(x)

        clusters <- y$clusters
        sample <- y$sample
    }

    prop.list <- getTransformedProps(clusters, sample)

    Proportions <- as.vector(t(prop.list$Proportions))
    Samples <- rep(colnames(prop.list$Proportions), nrow(prop.list$Proportions))
    Clusters <- rep(rownames(prop.list$Proportions),
                    each=ncol(prop.list$Proportions))

    plotdf <- data.frame(Samples=Samples, Clusters=Clusters,
                         Proportions=Proportions)

    ggplot(plotdf,aes(x=Samples,y=Proportions,fill=Clusters)) +
        geom_bar(stat="identity") +
        theme(axis.text.x = element_text(size=12),
            axis.text.y = element_text(size=12),
            axis.title = element_text(size=14),
            legend.text = element_text(size=12),
            legend.title = element_text(size=14))

}
