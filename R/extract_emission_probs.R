#! /usr/bin/Rscript

library(mHMMbayes)

args <- commandArgs(trailingOnly = TRUE)
model_fname<-args[1]
out_fname<-args[2]

args = commandArgs()

scriptName = args[substr(args,1,7) == '--file=']

if (length(scriptName) == 0) {
  scriptName <- rstudioapi::getSourceEditorContext()$path
} else {
  scriptName <- substr(scriptName, 8, nchar(scriptName))
}

pathName = substr(
  scriptName, 
  1, 
  nchar(scriptName) - nchar(strsplit(scriptName, '.*[/|\\]')[[1]][2])
)

# Load utility functions
source(paste0(pathName, "/obtain_emiss_pois.R"))

loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}

out<-loadRData(model_fname)
emissions<-obtain_emiss_pois(out, level="group")

State<-c()
Electrode<-c()
Mean<-c()
Var<-c()
for(j in 1:nrow(emissions[1]$el1)) {
  for(i in 1:length(emissions)) {
    State<-c(State, j)
    Electrode<-c(Electrode, i)
    Mean<-c(Mean, emissions[i][[paste0('el',i)]][j,1])
    Var<-c(Var, emissions[i][[paste0('el',i)]][j,2])
  }
}
df=data.frame(State, Electrode, Mean, Var)
write.csv(df, out_fname, row.names = FALSE)
