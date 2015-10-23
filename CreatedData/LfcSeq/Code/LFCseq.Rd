\name{LFCseq}
\alias{LFCseq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute probability of being not differentially expressed for each feature
}
\description{
Use this function to normalize read counts and compute probability of
being not differentially expressed for each feature.
}
\usage{
LFCseq(countsTable, condsAB, norm.method = "npseq", nloc = 50, pares = T)
}

\arguments{
  \item{countsTable}{Matrix or data.frame containing read counts for each feature. 
}
  \item{condsAB}{A vector containing experimental condition for each sample.
}
  \item{norm.method}{Normalization method, "rpm", "npseq" and "deseq".
}
  \item{nloc}{Number of features included in the each local set.
}
  \item{pares}{Logical. 
    Default is "TRUE" which means log2 fold changes within conditions are computed by separating
    the replicates in one condition equally. Pares="FAlSE" means log2 fold changes are computed
    for every pair of replicates within the same condition.
}
}

\value{

\item{pval}{Probability of being not differentially expressed  for each feature.}


}
\references{
Lin B., Zhang L-F., Chen X. (2014). LFCseq: a nonparametirc approach for differential expression  analysis of RNA-seq data. 
}
\author{
Lin Bingqing
}

\examples{
#subset of MAQC dataset.
data(maqc)
condsAB <- c(rep("brain", 7), rep("UHR", 7))
prob <- LFCseq(countsTable = maqc, condsAB = condsAB)$pval
prob
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
