# Run one original workflow with external private inputs and isolated outputs.
# Rscript --vanilla run_figure.R list
# BMB_DATA_ROOT must identify the original project/data tree for execution.
argv <- commandArgs(trailingOnly=TRUE)
full_args <- commandArgs(trailingOnly=FALSE)
this <- sub('^--file=','',grep('^--file=',full_args,value=TRUE)[1])
code_root <- dirname(normalizePath(this,winslash='/',mustWork=TRUE))
catalog <- read.csv(file.path(code_root,'figure_code_index.csv'),stringsAsFactors=FALSE)
key <- if(length(argv)) argv[1] else 'list'
if(key=='list') {
  print(catalog[,c('key','figures','role')],row.names=FALSE)
  quit(status=0)
}
job <- catalog[catalog$key==key,,drop=FALSE]
if(nrow(job)!=1) stop('Unknown module. Run with list to see available keys.')
data_root <- Sys.getenv('BMB_DATA_ROOT')
if(!nzchar(data_root)) stop('Set BMB_DATA_ROOT to the private original project/data folder. This repository contains no datasets.')
data_root <- normalizePath(data_root,winslash='/',mustWork=TRUE)
output_base <- Sys.getenv('BMB_OUTPUT_ROOT',file.path(code_root,'outputs'))
dir.create(output_base,recursive=TRUE,showWarnings=FALSE)
output_base <- normalizePath(output_base,winslash='/',mustWork=TRUE)
output_root <- file.path(output_base,key)
dir.create(output_root,recursive=TRUE,showWarnings=FALSE)
working_dir <- normalizePath(file.path(data_root,job$working_directory),winslash='/',mustWork=TRUE)
setwd(working_dir)
ctl <- new.env(parent=emptyenv())
ctl$current_source <- file.path(code_root,job$original_source)
aliases <- if(requireNamespace('jsonlite',quietly=TRUE)) {
  jsonlite::fromJSON(file.path(code_root,'input_aliases.json'))
} else list()

slash <- function(p) gsub('\\','/',p,fixed=TRUE)
collapse <- function(p) {
  parts <- strsplit(slash(p),'/',fixed=TRUE)[[1]]; keep <- character()
  for(part in parts) {
    if(part %in% c('','.')) next
    if(part=='..') {if(length(keep)) keep <- head(keep,-1)} else keep <- c(keep,part)
  }
  paste(keep,collapse='/')
}
relative <- function(p) {
  p <- slash(p)
  for(prefix in c(paste0(code_root,'/original_sources/'),paste0(code_root,'/'),paste0(data_root,'/'))) {
    if(startsWith(p,prefix)) return(collapse(substring(p,nchar(prefix)+1)))
  }
  marker <- regexpr('/BMB_AUT_MS/',p,fixed=TRUE)
  if(marker[1]>0) return(collapse(substring(p,marker[1]+nchar('/BMB_AUT_MS/'))))
  if(grepl('^([A-Za-z]:/|/)',p)) return(basename(p))
  cwd_relative <- if(identical(slash(getwd()),data_root)) '' else substring(slash(getwd()),nchar(data_root)+2)
  collapse(file.path(cwd_relative,p))
}
read_path <- function(p,must_exist=TRUE) {
  if(!is.character(p) || length(p)!=1 || !nzchar(p)) return(p)
  rel <- relative(p)
  tail_curated <- sub('^.*Manuscript_Figure_Code/','',rel)
  candidates <- c(file.path(output_root,basename(rel)),p,file.path(data_root,rel))
  if(rel %in% names(aliases)) candidates <- c(candidates,file.path(data_root,aliases[[rel]]))
  if(tail_curated %in% names(aliases)) candidates <- c(candidates,file.path(data_root,aliases[[tail_curated]]))
  candidates <- c(candidates,
    file.path(data_root,'0_Submission/4_NatComm/Manuscript_Figure_Code',tail_curated),
    file.path(Sys.getenv('BMB_CURATED_ROOT',file.path(data_root,'github_release_package')),tail_curated))
  if(basename(rel)=='pcnm_VPA.rdata') candidates <- c(candidates,file.path(data_root,'10_PCNM_VPA/pcnm_VPA.rdata'))
  existing <- candidates[file.exists(candidates)]
  if(length(existing)) return(normalizePath(existing[1],winslash='/',mustWork=TRUE))
  if(must_exist) stop('Missing external input: ',p,'\nLooked relative to ',data_root)
  NA_character_
}
write_path <- function(p) {
  if(!is.character(p) || length(p)!=1 || !nzchar(p)) return(p)
  # Isolate legacy output names per workflow. Never overwrite external inputs.
  dest <- file.path(output_root,basename(slash(p)))
  dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
  dest
}
if('--check' %in% argv) {
  inputs <- read.csv(file.path(code_root,'required_input_paths.csv'),stringsAsFactors=FALSE)
  inputs <- inputs[inputs$key==key,,drop=FALSE]
  inputs$resolved_path <- vapply(inputs$literal_path,read_path,character(1),must_exist=FALSE)
  inputs$found <- !is.na(inputs$resolved_path)
  write.csv(inputs,file.path(output_root,'input_check.csv'),row.names=FALSE)
  print(inputs[,c('literal_path','found')],row.names=FALSE)
  message('This checks literal inputs only; dynamically constructed file paths also require the private data tree.')
  quit(status=0)
}

# These bindings change file locations only. The original analysis statements
# are evaluated in a child environment, so original rm(list=ls()) commands do
# not remove the I/O configuration or wrappers.
io <- new.env(parent=globalenv())
io$commandArgs <- function(trailingOnly=FALSE) {
  if(trailingOnly) return(character())
  a <- base::commandArgs(trailingOnly=FALSE)
  a[grepl('^--file=',a)] <- paste0('--file=',ctl$current_source)
  a
}
io$setwd <- function(dir) base::setwd(read_path(dir))
io$readRDS <- function(file,...) base::readRDS(read_path(file),...)
io$load <- function(file,envir=parent.frame(),...) base::load(read_path(file),envir=envir,...)
io$read.csv <- function(file,...) {
  actual <- read_path(file); dots <- list(...)
  # Preserve BOM-bearing CSV column names across locales, without changing values.
  if(is.character(actual) && is.null(dots$fileEncoding)) {
    con <- file(actual,'rb'); head <- readBin(con,'raw',n=3); close(con)
    if(identical(head,as.raw(c(239,187,191)))) dots$fileEncoding <- 'UTF-8-BOM'
  }
  do.call(utils::read.csv,c(list(file=actual),dots))
}
io$read.table <- function(file,...) utils::read.table(read_path(file),...)
io$read.delim <- function(file,...) utils::read.delim(read_path(file),...)
io$read_csv <- function(file,...) readr::read_csv(read_path(file),...)
io$read_q2biom <- function(file,...) qiime2R::read_q2biom(read_path(file),...)
io$read_excel <- function(path,...) readxl::read_excel(read_path(path),...)
io$write.csv <- function(x,file,...) utils::write.csv(x,write_path(file),...)
io$write.table <- function(x,file='',...) utils::write.table(x,write_path(file),...)
io$saveRDS <- function(object,file,...) base::saveRDS(object,write_path(file),...)
io$save <- function(...,file,envir=parent.frame()) base::save(...,file=write_path(file),envir=envir)
io$ggsave <- function(filename,...,path=NULL) ggplot2::ggsave(write_path(if(is.null(path)) filename else file.path(path,filename)),...)
io$pdf <- function(file='Rplots.pdf',...) grDevices::pdf(write_path(file),...)
io$png <- function(filename='Rplot%03d.png',...) grDevices::png(write_path(filename),...)
io$dir.create <- function(path,...) base::dir.create(write_path(path),...)
io$source <- function(file,local=parent.frame(),...) {
  actual <- read_path(file); previous <- ctl$current_source
  ctl$current_source <- actual
  on.exit(ctl$current_source <- previous)
  base::source(actual,local=local,...)
}
analysis <- new.env(parent=io)
# Original GWAS preparation notebooks assume tidyverse/vegan already attached.
suppressPackageStartupMessages({library(tidyverse);library(vegan)})
exprs <- parse(file.path(code_root,job$script),keep.source=FALSE)
grDevices::pdf(file.path(output_root,'notebook_plots.pdf'))
for(expr in exprs) {
  result <- withVisible(eval(expr,envir=analysis))
  if(result$visible && inherits(result$value,c('ggplot','patchwork','gg'))) print(result$value)
}
grDevices::dev.off()
capture.output(sessionInfo(),file=file.path(output_root,'session_info.txt'))
message('Completed ',key,'. Outputs: ',output_root)
