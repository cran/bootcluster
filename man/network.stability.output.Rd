% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotcode.R
\name{network.stability.output}
\alias{network.stability.output}
\title{Plot method for objests from threshold.select}
\usage{
network.stability.output(input, optimal.only = FALSE)
}
\arguments{
\item{input}{a \code{list} of results from function \code{threshold.select}}

\item{optimal.only}{a \code{logical} value indicating whether only plot the network with optimal threshold or not. The default is False, generating
all network figures with a large number of nodes could take some time.}
}
\value{
Plot of network figures
}
\description{
Plot method for objests from threshold.select
}
\details{
\code{network.stability.output} is used to generate a series of network plots based on the given \code{threshold.seq},where the nodes are
colored by the level of stability. The network with optimal
threshold value selected by function \code{threshold.select} is colored as red.
}
\examples{
\donttest{
set.seed(1)
data(wine)
x0 <- wine[1:50,]

mytest<-threshold.select(data.input=x0,threshold.seq=seq(0.1,0.5,by=0.05), B=20, 
cor.method='pearson',large.size=0,
PermuNo = 10,
no_cores=1,
scheme_2 = FALSE)
network.stability.output(mytest)
}
}
\references{
A framework for stability-based module detection in correlation graphs.
Mingmei Tian,Rachael Hageman Blair,Lina Mu, Matthew Bonner, Richard Browne and Han Yu.
}
\author{
Mingmei Tian
}
