% File cinf-getprops.R
\name{get_prop}
\alias{get_prop}
\title{Get data frame with property values}
\description{
  Returns a data frame for a specified property name
}
\usage{
get_prop(moldbase, propname)
}
\arguments{
	\item{moldbase}{Molecular database}
	\item{propname}{Property name}
}
\value{Data frame with property values}
\examples{
data("alkan")
get_prop(mdb, "bp")
}
