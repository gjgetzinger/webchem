#' Get Wikidata Item ID
#'
#' @import jsonlite httr
#'
#' @param query character; The searchterm
#' @param language character; the language to search in
#' @param first logical; If TRUE return only first result.
#' @param verbose logical; print message during processing to console?
#'
#' @return A character vector with the item ID and the additional attribute \code{matched}  (the matched
#' label).
#'
#' @note Only matches in labels are returned.
#'
#' @author Eduard Szoecs, \email{eduardszoecs@@gmail.com}
#' @export
#' @examples
#' \dontrun{
#' get_wdid('Triclosan', language = 'en')
#' get_wdid('DDT', language = 'en')
#' get_wdid('DDT', language = 'en', first = TRUE)
#'
#' # multiple inpus
#' comps <- c('Triclosan', 'Glyphosate')
#' sapply(comps, get_wdid, language = 'en')
#' }
get_wdid <- function(query, language, first = FALSE, verbose = TRUE){
  # language <-  'en'
  # query <- 'Triclosan'
  if (length(query) > 1) {
    stop('Cannot handle multiple input strings.')
  }
  limit <-  50
  qurl <- paste0("wikidata.org/w/api.php?action=wbsearchentities&format=json&type=item")
  qurl <- paste0(qurl, "&language=", language, "&limit=", limit, "&search=", query)
  if (verbose)
    message('Querying ', qurl)
  Sys.sleep(0.3)
  cont <- fromJSON(content(GET(qurl, user_agent('webchem (https://github.com/ropensci/webchem)')), 'text'))
  search <- cont$search
  if (length(search) == 0) {
    if (verbose)
      message('Substance not found! Returing NA. \n')
    return(NA)
  }
  # use only matches on label
  search <- search[search$match$type == 'label', ]
  search <- search[tolower(search$match$text) == tolower(query), ]

  if (first) {
    search <- search[1, ]
  }

  out <- search$id
  attr(out, "matched") <- search$label
  return(out)
}

#! Use SPARQL to search of chemical compounds (P31)?! For a finer / better search?



#' Retrieve Indentifiers from wikidata
#'
#' @import jsonlite
#'
#' @param id character; identifier, as returned by \code{\link{get_wdid}}
#' @param verbose logical; print message during processing to console?
#'
#' @return A list of identifiers. Currently these are 'smiles', 'cas', 'cid', 'einecs', 'csid', 'inchi', 'inchikey',
#' 'drugbank', 'zvg', 'chebi', 'chembl', 'unii'
#'
#' @note Only matches in labels are returned.
#'
#' @seealso \code{\link{get_wdid}}
#'
#' @author Eduard Szoecs, \email{eduardszoecs@@gmail.com}
#' @export
#' @examples
#' \dontrun{
#'  id <- "Q408646" # Triclosan
#'  wd_ident(id)
#' }
wd_ident <- function(id, verbose = TRUE){
  # id <- 'Q408646'
  if (length(id) > 1) {
    stop('Cannot handle multiple input strings.')
  }
  baseurl <- 'https://query.wikidata.org/sparql?format=json&query='
  props <- c('P233', 'P231', 'P662', 'P232', 'P661', 'P234', 'P235', 'P715', 'P679',
             'P683', 'P592', 'P652')
  names <- c('smiles', 'cas', 'cid', 'einecs', 'csid', 'inchi', 'inchikey',
             'drugbank', 'zvg', 'chebi', 'chembl', 'unii')

  sparql_head <- paste('PREFIX wd: <http://www.wikidata.org/entity/>',
    'PREFIX wdt: <http://www.wikidata.org/prop/direct/>',
    'SELECT * WHERE {')
  sparql_body <- paste(paste0('OPTIONAL{wd:', id, ' wdt:', props, ' ?', names, ' .}'),
        collapse = ' ')
  sparql <- paste(sparql_head, sparql_body, '}')
  qurl <- paste0(baseurl, sparql)
  qurl <- URLencode(qurl)
  Sys.sleep(0.3)
  if (verbose)
    message('Querying ', qurl)
  tmp <- fromJSON(qurl)
  out <- tmp$results$bindings
  if (length(out) == 0) {
    if (verbose)
      message('Not found! Returing NA. \n')
    return(NA)
  }

  out <- lapply(out, '[[', 'value')
  return(out)
}
