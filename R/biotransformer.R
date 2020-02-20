#' BioTransformer
#'
#' BioTransformer is a software
#' tool that predicts small molecule metabolism in mammals, their gut
#' microbiota, as well as the soil/aquatic microbiota. BioTransformer also
#' assists scientists in metabolite identification, based on the metabolism
#' prediction.
#'
#' Download the jar file from: \url{https://bitbucket.org/djoumbou/biotransformerjar/src/master/}
#' Cite: Djoumbou-Feunang Y, Fiamoncini J, de la Fuente AG, Manach C, Greiner R,
#' and Wishart DS; BioTransformer: A Comprehensive Computational Tool for Small
#' Molecule Metabolism Prediction and Metabolite Identification; Journal of
#' Cheminformatics 201911:2; DOI: 10.1186/s13321-018-0324-5.
#'
#' @param jar Path to the biotransformer jar file
#' @param annotate (T/F) Search PubChem for each product, and store with CID and
#'   synonyms, when available.
#' @param btType The type of description: Type of biotransformer
#'   human super transformer (`superbio`, or `allHuman`), Environmental
#'   microbial(`env`).
#' @param formulas Semicolon-separated list of formulas of compounds to identify
#' @param ismiles The input, which can be a SMILES string
#' @param task The task to be permed: pred for prediction, or cid for compound
#'   identification
#' @param masses Semicolon-separated list of masses of compounds to identify
#' @param nsteps The number of steps for the prediction. This option can be set
#'   by the user for the EC-based, CYP450, Phase II, and Environmental microbial
#'   biotransformers. The default value is 1.
#' @param mTolerence Mass tolerance for metabolite identification (default is
#'   0.01).
#'
#' @return A tibble containing predicted metabolite data
#' @export
#'
biotransformer <-
  function(jar = NULL,
           annotate = FALSE,
           btType = NULL,
           ismiles = NULL,
           task = NULL,
           formulas = NULL,
           masses = NULL,
           nsteps = 1,
           mTolerence = 0.01) {
    # Check input arguments
    if (is.null(ismiles)) {
      stop('No input SMILES string provided')
    }
    if(is.null(task) | !any(task %in% c('pred', 'cid'))){
      stop('tak must be either pred or cid')
    }
    if(is.null(btType) | !any(btType %in% c('env', 'allHuman', 'superbio'))){
      stop('btType must be one of: env, allHuman or superbio')
    }
    if(task == 'cid' & all(is.null(c(masses, formulas)))){
      stop('Must provide masses or formulas when task = cid')
    }
    if(length(ismiles)>1){
      stop('Multiple SMILES strings input, only single queries currenty supported.')
    }
    if (is.null(jar)) {
      # TODO add REST API interface
      stop('POST API interface not yet implemented, must use properly installed jar file.')
      # x_url <- 'http://biotransformer.ca/queries.json'
      # task_type <- switch(task, pred = 'PREDICTION', cid = 'IDENTIFICATION')
      # biotransformer_option <- switch(btType, env ='ENVMICRO', allHuman = 'ALLHUMAN', superbio = 'SUPERBIO')
      # query_input <- paste0('x\t', ismiles)
      #
      # body <- jsonlite::toJSON(list(
      #   biotransformer_option = biotransformer_option,
      #   number_of_steps = nsteps,
      #   query_input = query_input,
      #   task_type = task_type
      # ), auto_unbox = T)
      # r <- httr::POST(x_url, body = body, encode = 'json', verbose())
      #
    } else {
      jar <- normalizePath(path = jar, mustWork = T)
      csvoutput <- tempfile(tmpdir = getwd(), fileext = '.csv')

      arguments <- paste(
        '-jar',jar,
        '-k',task,
        '-b',btType,
        '-ismi',
        shQuote(string = ismiles),
        '-ocsv',csvoutput,
        '-s',nsteps
      )
      if (!is.null(masses)) {
        arguments <-
          paste(arguments, paste0('-m "', paste0(masses, collapse = ';'), '"'), '-t', mTolerence)
      }
      if (!is.null(formulas)) {
        arguments <-
          paste(arguments, paste0('-f "', paste0(formulas, collapse = ';'), '"'), '-t', mTolerence)
      }
      if (annotate) {
        arguments <- paste(arguments, '-a')
      }
      std_out <- system2(command = 'java',
                         args = arguments, stdout = T)
      if(any(grepl('^The results were saved to the following file', std_out))) {
        rst <-
            readr::read_csv(file = csvoutput, col_types = "ccccicdddddccccccccccdd")
        } else {
        stop(paste0('Biotransformer failed for input: ', ismiles, sep = ''))
      }
    }
  return(rst)
  }