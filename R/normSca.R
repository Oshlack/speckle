#' Normalise and scale counts matrix
#' 
#' This function is called by the \code{preprocess} function and performs
#' log-normalisation and scaling of a counts matrix.
#' 
#' @param x a matrix of counts
#' @param lib.size the library size 
#' @param log logical, indicates whether the output should be on the log2 scale
#' or counts scale. Default is TRUE.
#' @param prior.count The prior count to add if the data is log2 normalised.
#' Default is a small count of 0.5.
#' 
#' @return a dataframe of log-normalised and scaled counts
#' 
#' @export normSca
#'
#' @examples
#' 
#' y <- matrix(rnbinom(1000, size=2, mu= 20),ncol=10)
#' colnames(y)<- paste("Cell",1:10, sep="")
#' row.names(y)<-paste("Gene",1:100,sep="")
#' norm.data <- normSca(y,lib.size=colSums(y))
#' 
#' #Visualise the counts vs scaled data
#' boxplot(y)
#' boxplot(norm.data)
#' 
normSca<-function(x, lib.size=lib.size, log = TRUE, prior.count = 0.5)
#    Normalise and scale counts matrix
#    Xinyi Jin and Belinda Phipson
#    17 February 2021
#    Modified 17 February 2021    
{
    x <- as.matrix(x)
    # log normalise
    normalisedVal<- normCounts(x, log = log, prior.count = prior.count, 
                               lib.size=lib.size)
    # scale
    #scaledVal<-t(scale(t(normalisedVal)))
    #scaledVal
    
    normalisedVal
}


