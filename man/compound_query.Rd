% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compound_query.R
\name{compound_query}
\alias{compound_query}
\title{Compound query}
\usage{
compound_query(query_string)
}
\arguments{
\item{query_string}{A character string with a compound/drug name}
}
\value{
A list of tibbles that contain general information, target information, pathways where compound is involved, diseases associated with compound/drug, common uses for a drug of interest, and KEGG drug groups for compound of interest.
}
\description{
This function will perform a query of a given compound/drug in KEGG or PubChem and return available data based on queries to the REST APIs for either database.
}
\examples{
compound_query("aspirin")
}
