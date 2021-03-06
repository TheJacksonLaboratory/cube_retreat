#####
## AUTHOR: BRIAN M. BOT
## ORGANIZATION: SAGE BIONETWORKS
##
## ALLOW USERS TO LEVERAGE KNITR WHEN CONSTRUCTING SYNAPSE WIKI CONTENT
#####
## PARAMETERS:
##   file: path to a local .Rmd file which to knit
##   owner: a Synapse object which will own the resulting WikiPage (usually a Project, Folder, or File)
##   parentWikiId (optional): if the resulting WikiPage is to be a subpage of another WikiPage, this is the id for the parent WikiPage (NOTE: owner is still required)
##   wikiName (optional): a title for the resulting WikiPage - will default to the file name without the .Rmd extension
##   overwrite (default = FALSE): only if owner specified and parentWikiId is NULL - flag for whether or not to overwrite the previous root WikiPage (if it exists)
#####
## VALUE:
##   a WikiPage object as defined in the synapseClient
#####

knit2synapse <- function(file, owner, parentWikiId=NULL, wikiName=NULL, overwrite=FALSE){
  require(knitr)
  require(RCurl)
  require(stringr)
  require(tools)
  require(synapser)
  
  ## CHECK TO MAKE SURE FILE EXISTS
  file <- path.expand(file)
  if( !file.exists(file) ){
    stop(sprintf("file %s does not exist at this location:\n", basename(f), f))
  }
  
  ## IF NO WIKI NAME GIVEN, DEFAULT TO FILE NAME WITHOUT EXTENSION
  fName <- basename(file_path_sans_ext(file))
  if( is.null(wikiName) ){
    wikiName <- fName
  }
  
  ## IF OWNER IS CHARACTER, TRY TO GET FROM SYNAPSE
  if( is.character(owner) & length(owner) == 1 ){
    owner <- synGet(owner, downloadFile=FALSE)
  }
  
  #####
  ## SET SYNAPSE-SPECIFIC MARKDOWN HOOKS
  #####
  render_markdown()
  
  ## PLOTS
  hook_synapseMdSyntax_plot <- function(x, options){
    synPlotMdOpts <- character()
    
    ## SET URL ENCODING STRINGS
    urlEncodings <- c('{' = "%7B",
                      '}' = "%7D",
                      '-' = "%2D",
                      '_' = "%5F",
                      '.' = "%2E",
                      '!' = "%21",
                      '~' = "%7E",
                      '*' = "%2A",
                      '`' = "%60",
                      '\'' = "%27",
                      '(' = "%28",
                      ')' = "%29",
                      '[' = "%5B",
                      ']' = "%5D",
                      ':' = "%3A",
                      ';' = "%3B",
                      '\n' = "%0A",
                      '\r' = "%0D",
                      '/' = "%2F",
                      '?' = "%3F",
                      '&' = "%26",
                      '=' = "%3D",
                      '+' = "%2B",
                      ',' = "%2C",
                      '#' = "%23",
                      '$' = "%24")
    
    ## CHECK FOR ALIGN OPTION BEING SET
    if( any(names(options) == "align") ){
      ## CHECKS FOR ALIGN OPTION
      if( !is.character(options$align) ){
        stop("align must be one of none, left, right, or center")
      }
      if( !(options$align %in% c("none", "left", "right", "center")) ){
        stop("align must be one of none, left, right, or center")
      }
      synPlotMdOpts <- paste(synPlotMdOpts, "&align=", options$align, sep="")
    } else{
      synPlotMdOpts <- paste(synPlotMdOpts, "&align=none", sep="")
    }
    
    ## CHECK FOR SCALE OPTION BEING SET
    if( any(names(options) == "scale") ){
      ## RANGE CHECKS FOR SCALE OPTION
      if( !is.numeric(options$scale) ){
        stop("scale option must be numeric")
      }
      if( options$scale <= 0 | options$scale > 500 ){
        stop("scale option must be greater than 0 and less than 500")
      }
      
      synPlotMdOpts <- paste(synPlotMdOpts, "&scale=", options$scale, sep="")
    } else{
      synPlotMdOpts <- paste(synPlotMdOpts, "&scale=100", sep="")
    }
    
    paste("${image?fileName=", curlPercentEncode(basename(paste(x, collapse=".")), codes=urlEncodings), synPlotMdOpts, "}\n", sep="")
  }
  knit_hooks$set(plot=hook_synapseMdSyntax_plot)
  opts_chunk$set(tidy=FALSE)
  opts_chunk$set(error=FALSE)
  
  
  ## CREATE TEMPORARY OUTPTU DIRECTORY FOR MD AND PLOTS
  knitDir <- tempfile(pattern="knitDir")
  dir.create(knitDir)
  knitPlotDir <- file.path(knitDir, "plots/")
  dir.create(knitPlotDir)
  opts_chunk$set(fig.path = knitPlotDir)
  
  mdName <- file.path(knitDir, paste(fName, ".md", sep=""))
  
  mdFile <- knit(file,
                 envir=parent.frame(n=2),
                 output=mdName)
  att <- list.files(knitPlotDir, full.names=TRUE)
  
  if( is.null(parentWikiId) ){ ## doesn't have a parentWiki
    if( length(att) > 0 ){ ## has attachments
      w <- WikiPage(owner=owner, 
                    title=wikiName, 
                    attachments=as.list(att),
                    markdown=readChar(mdFile, file.info(mdFile)$size))
    } else{ ## doesn't have attachments
      w <- WikiPage(owner=owner, 
                    title=wikiName, 
                    markdown=readChar(mdFile, file.info(mdFile)$size))
    }
    
    if( overwrite ){
      ## TRY TO STORE
      tmp <- try(synStore(w), silent=TRUE)
      if( class(tmp) == "try-error" ){
        tmp <- synGetWiki(owner)
        tmp <- synDelete(tmp)
        w <- synStore(w)
      } else{
        w <- tmp
      }
    } else{
      w <- synStore(w)
    }
    
  } else{ ## has a parentWiki
    if( length(att) > 0 ){ ## has attachments
      w <- Wiki(owner=owner, 
                    title=wikiName, 
                    attachments=as.list(att),
                    markdown=readChar(mdFile, file.info(mdFile)$size),
                    parentWikiId=parentWikiId)
    } else{ ## doesn't have attachments
      w <- Wiki(owner=owner, 
                    title=wikiName, 
                    markdown=readChar(mdFile, file.info(mdFile)$size),
                    parentWikiId=parentWikiId)
    }
    w <- synStore(w)
  }
  cat(paste("built wiki: '", wikiName, "'\n", sep=""))
  return(w)
}
