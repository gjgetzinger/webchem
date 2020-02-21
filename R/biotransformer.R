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
#' @param molinput A filename for the input mol file (default = NULL).
#' @param sdfinput A filename for the input sdf file (default = NULL).
#' @param csvoutput A filename for the output csv file (default = NULL).
#' @param sdfoutput A filename for the output sdf file (default = NULL).
#' @param task The task to be permed: pred for prediction, or cid for compound
#'   identification
#' @param masses Semicolon-separated list of masses of compounds to identify
#' @param nsteps The number of steps for the prediction. This option can be set
#'   by the user for the EC-based, CYP450, Phase II, and Environmental microbial
#'   biotransformers. The default value is 1.
#' @param mTolerence Mass tolerance for metabolite identification (default is
#'   0.01).
#'
#' @return A path to the output file containing predicted metabolite data
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
           mTolerence = 0.01,
           csvoutput = NULL,
           sdfoutput = NULL,
           molinput = NULL,
           sdfinput = NULL
           ) {
    # Check input arguments
    if (is.null(ismiles) & is.null(molinput) & is.null(sdfinput)) {
      stop('No input SMILES string or structure file provided.')
    }
    if(is.null(task) | !any(task %in% c('pred', 'cid'))){
      stop('task must be either pred or cid')
    }
    if(is.null(btType) | !any(btType %in% c('env', 'allHuman', 'superbio'))){
      stop('btType must be one of: env, allHuman or superbio')
    }
    if(task == 'cid' & all(is.null(c(masses, formulas)))){
      stop('Must provide masses or formulas when task = cid')
    }
    if(!is.null(ismiles) & length(ismiles)>1){
      stop('Multiple SMILES strings input, only single queries currenty supported.')
    }
    if( is.null(csvoutput) & is.null(sdfoutput)){
      stop('Must provide an output file name (.csv or .sdf). ')
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
      out_type <- ifelse(is.null(csvoutput), '-osdf', '-ocsv')
      out_file <- switch(
        out_type,
        '-osdf' = normalizePath(sdfoutput, mustWork = F),
        '-ocsv' = normalizePath(csvoutput, mustWork = F)
      )
      in_type <-
        ifelse(!is.null(sdfinput),
               '-isdf',
               ifelse(!is.null(molinput), '-imol', ifelse(!is.null(ismiles), '-ismi', stop(
                 'Input type not recognized'
               ))))
      in_arg <- switch(in_type,
                       '-ismi' = shQuote(string = ismiles),
                       '-imol' = normalizePath(molinput, mustWork = T),
                       '-isdf' = normalizePath(sdfinput, mustWork = T)
                       )
      arguments <- paste(
        '-jar',jar,
        '-k',task,
        '-b',btType,
        in_type, in_arg,
        out_type, out_file,
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
      system2(command = 'java',
              args = arguments,
              stdout = T)
      message(out_file)
      return(out_file)
    }
  }