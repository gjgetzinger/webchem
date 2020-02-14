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
#' @param annotate (T/F) Search PuChem for each product, and store with CID and
#'   synonyms, when available.
#' @param btType The type of description: Type of biotransformer - EC-based
#'   (ecbased), CYP450 (cyp450), PhaseII (phaseII), Human gut microbial (hgut),
#'   human super transformer (superbio, or allHuman), Environmental
#'   microbial(envimicro).If option -m is enabled, the only valid biotransformer
#'   types are allHuman, superbio and env.
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
#' @return
#' @export
#'
biotransformer <-
  function(jar = NULL,
           annotate = FALSE,
           btType = c('ecbased',
                      'cyp450',
                      'phaseII',
                      'hgut',
                      'superbio',
                      'allHuman',
                      'envimicro'),
           ismiles = NULL,
           task = c('pred', 'cid'),
           formulas = NULL,
           masses = NULL,
           nsteps = 1,
           mTolerence = 0.01) {
    # Check input arguments
    if (is.null(ismiles)) {
      stop('No input SMILES string provided')
    }
    if(task == 'cid' & all(is.null(c(masses, formulas)))){
      stop('Must provide masses or formulas when task = cid')
    }
    if (is.null(jar)) {
      # TODO add REST API interface
      stop('REST API interface not yet implemented, must use properly installed jar file.')
    } else {
      jar <- normalizePath(path = jar, mustWork = T)
      csvoutput <- tempfile(pattern = 'file')
      on.exit(unlink(csvoutput))
      args <- paste(
        '-jar',
        jar,
        '-k',
        task,
        '-b',
        btType,
        '-ismi',
        paste0('"', ismiles, '"'),
        '-ocsv',
        csvoutput,
        '-s',
        nsteps
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
                         args = arguments,
                         stdout = T)

      message(paste0(std_out, '\n'))
      if (any(grepl('Successfully completed metabolism', std_out))) {
        rst <- readr::read_csv(file = csvoutput)
      } else {
        stop(paste0('Biotransformer failed for input: ', ismiles, sep = ''))
      }
    }
  return(rst)
  }