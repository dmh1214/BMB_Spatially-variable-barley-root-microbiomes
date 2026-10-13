# Figures: Fig5
# Original source: 0_Submission/4_NatComm/Manuscript_Figure_Code/scripts/fig5d_plot_archived_statistics.R
# Analysis statements retained; input/output locations are handled by run_figure.R.

# Replot using original manuscript statistics for labels.
# This is a plotting-only mode; it does not claim to recompute archived rho.
args <- commandArgs(trailingOnly=FALSE)
this <- sub('^--file=', '', grep('^--file=',args,value=TRUE)[1])
script_dir <- dirname(normalizePath(this,winslash='/',mustWork=TRUE))
project_root <- normalizePath(file.path(script_dir,'..'),winslash='/',mustWork=TRUE)

# Reuse the original count-to-CV data preparation and regression calculations,
# stopping before plotting/saving. Suppress writes from the source workflow.
expressions <- parse(file.path(script_dir,'fig5d_individual_stability.R'),keep.source=FALSE)
for (expr in expressions) {
  operation <- if(is.call(expr)) as.character(expr[[1]])[1] else ''
  if(operation=='<-' && identical(as.character(expr[[2]]),'p')) {
    original_plot_expression <- expr
    break
  }
  if(operation %in% c('write.csv','ggsave')) next
  eval(expr,envir=.GlobalEnv)
}
recalculated_individuals <- indiv_slopes
historical_individuals <- read.csv(
  file.path(project_root,'reference_results','fig5d_manuscript_individual_statistics.csv'),
  check.names=FALSE
)
historical_aggregates <- read.csv(
  file.path(project_root,'reference_results','fig5d_manuscript_aggregate_statistics.csv'),
  check.names=FALSE
)
stopifnot(!anyDuplicated(historical_individuals[c('Organism','Trial','Biomarker')]))
stopifnot(setequal(top_biomarkers$Biomarker,historical_individuals$Biomarker))

# Labels are read from existing historical CSVs, never presented as new tests.
top_biomarkers <- historical_individuals
plot_indiv <- plot_indiv %>%
  select(-biomarker_label) %>%
  left_join(historical_individuals %>% select(Organism,Trial,Biomarker,biomarker_label),
            by=c('Organism','Trial','Biomarker'))
stopifnot(!anyNA(plot_indiv$biomarker_label))
agg_stats <- historical_aggregates %>%
  transmute(
    panel=paste(Organism,Trial,sep=' | '),
    slope=slope_cv,
    p_value=p_lm_cv,
    rho=rho_cv,
    label=paste0('Aggregate: rho = ',sprintf('%.2f',rho_cv),
                 '\nP = ',format.pval(p_lm_cv,digits=2,eps=1e-3))
  )
eval(original_plot_expression,envir=.GlobalEnv)

archive_output <- file.path(project_root,'outputs','fig5d_archived_plot')
dir.create(archive_output,recursive=TRUE,showWarnings=FALSE)
ggsave(file.path(archive_output,'fig5d_original_statistics_plot.pdf'),p,
       width=10.5,height=5.2,bg='white')
ggsave(file.path(archive_output,'fig5d_original_statistics_plot.png'),p,
       width=10.5,height=5.2,dpi=300,bg='white')
write.csv(historical_individuals %>% mutate(statistic_source='archived manuscript result'),
          file.path(archive_output,'individual_statistics_used_for_labels.csv'),row.names=FALSE)
write.csv(historical_aggregates %>% mutate(statistic_source='archived manuscript result'),
          file.path(archive_output,'aggregate_statistics_used_for_labels.csv'),row.names=FALSE)
write.csv(recalculated_individuals,
          file.path(archive_output,'recalculated_individual_statistics_for_comparison.csv'),row.names=FALSE)
writeLines(c(
  'Plot labels use the original archived manuscript statistic tables.',
  'Points and regression lines use the original count-to-CV preparation and regression workflow.',
  'Archived rho values are not described as newly recomputed results.',
  'Final article composition and individually added legend text are manual.'
),file.path(archive_output,'STATISTIC_SOURCES.txt'))
message('Completed original-statistics plotting mode: ',archive_output)
