% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_drugbank.R
\name{parse_drugbank}
\alias{parse_drugbank}
\title{Parse DrugBank}
\usage{
parse_drugbank(drugbank_xml)
}
\arguments{
\item{drugbank_xml}{Full XML database file for DrugBank ([from here](https://www.drugbank.ca/releases/latest))}
}
\value{
A tibble for DrugBank information
}
\description{
Wrapper function to parse DrugBank XML file and generate a tibble for usage in DrugMineR. Built on the `dbparser` package. The data that is provided with `DrugMineR` was built using DrugBank v5.1.4.
}
\examples{
drugbank_parsed <- parse_drugbank("path/to/drugbank.xml")
}
