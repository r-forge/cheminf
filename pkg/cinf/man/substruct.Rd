% File cinf-mol
\name{substruct}
\alias{substruct}
\title{Extract substructure from a molecule}
\description{
  Extracts substructure from a molecule.
}
\usage{
substruct(mol, oldnum)
}
\arguments{
	\item{mol}{molecule}
	\item{oldnum}{list of the atoms defining the substructure to be extracted}
}
\value{Substructure}
\examples{
data("alkan")
mol <- mdb[[2]]
substruct(mol, c(1,2))
}
