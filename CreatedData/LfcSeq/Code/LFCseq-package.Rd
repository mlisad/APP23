\name{LFCseq-package}
\alias{LFCseq-package}
%\alias{LFCseq}
\docType{package}
\title{
A nonparametric approach for differential expression
analysis of RNA-seq
}
\description{
LFCseq is a package for differential expression analysis of RNA-seq via a nonparametric approach.
The main idea of LFCseq
is that the difference measured by log2 fold change between conditions is larger than
the difference within conditions. For each feature, when calculating the probability
of the feature being not differentially expressed, we only use log2 fold changes within
conditions for features whose mean read counts are close to the mean read counts of the
feature.
}

}
\author{
Lin Bingqing

Maintainer:Lin Bingqing <bqlin@ntu.edu.sg>
}
\references{
Lin B., Zhang L-F., Chen X. (2014). LFCseq: a nonparametirc approach for differential expression  analysis of RNA-seq data. 
}
\keyword{ package }
\seealso{
?LFCseq
}
\examples{
#subset of MAQC dataset.
data(maqc)
condsAB <- c(rep("brain", 7), rep("UHR", 7))
prob <- LFCseq(countsTable = maqc, condsAB = condsAB)
prob
}

