% File cinf-mol
\name{mol_get_ct}
\alias{mol_get_ct}
\title{Return connection table for a molecule}
\description{
  Returns connection table for a molecule.
}
\usage{
mol_get_ct(mol)
}
\arguments{
	\item{mol}{molecule}
}
\value{Connection table}
\examples{
data("alkan")
mol <- mdb[[1]]
mol_get_ct(mol)
}
