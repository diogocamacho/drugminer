% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kegg_query.R
\name{kegg_query}
\alias{kegg_query}
\title{KEGG drug query}
\usage{
kegg_query(query_string)
}
\arguments{
\item{query_string}{A character string with compound or drug name.}
}
\value{
The content of the POST request.
}
\description{
Extracts information from KEGG for a given drug or compound of interest using the KEGG REST API. The POST request is done using the `httr` package. Given that this is a character vector that is unformated, it is likely unreadable and it is recommended that you run \code{\link{extract_info}} on the output of this function.
}
\examples{
Query a drug:
aspirin <- kegg_query("aspirin") 

Query a compound:
nad_syn <- kegg_query("dpnh")
}
