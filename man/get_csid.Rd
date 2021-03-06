% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chemspider.R
\name{get_csid}
\alias{get_csid}
\title{ChemSpider ID from compound name, formula, SMILES, InChI or InChIKey}
\usage{
get_csid(query, from = "name", apikey = NULL, control = cs_control())
}
\arguments{
\item{query}{character; search term.}

\item{from}{character; the type of the identifier to convert from. Valid
values are \code{"name"}, \code{"formula"}, \code{"smiles"}, \code{"inchi"},
\code{"inchikey"}. The default value is \code{"name"}.}

\item{apikey}{character; your API key. If NULL (default),
\code{cs_check_key()} will look for it in .Renviron or .Rprofile.}

\item{control}{list; see details.}
}
\value{
Returns a data frame.
}
\description{
Query one or more compunds by name, formula, SMILES, InChI or InChIKey and
return a vector of ChemSpider IDs.
}
\details{
Queries by SMILES, InChI or InChiKey do not use \code{cs_control}
options. Queries by name use \code{order_by} and \code{order_direction}.
Queries by formula also use \code{datasources}. See \code{cs_control()} for a
full list of valid values for these control options.

\code{formula} can be expressed with and without LaTeX syntax.
}
\note{
An API key is needed. Register at \url{https://developer.rsc.org/}
for an API key. Please respect the Terms & conditions:
\url{https://developer.rsc.org/terms}.
}
\examples{
\dontrun{
get_csid("triclosan")
get_csid(c("carbamazepine", "naproxene","oxygen"))
get_csid("C2H6O", from = "formula")
get_csid("C_{2}H_{6}O", from = "formula")
get_csid("CC(O)=O", from = "smiles")
get_csid("InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)", from = "inchi")
get_csid("QTBSBXVTEAMEQO-UHFFFAOYAR", from = "inchikey")
}
}
\references{
\url{https://developer.rsc.org/compounds-v1/apis}
}
\author{
Eduard Szoecs, \email{eduardszoecs@gmail.com}

Tamas Stirling, \email{stirling.tamas@gmail.com}
}
