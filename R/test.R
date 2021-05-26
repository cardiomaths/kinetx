
#Load all source code
source("kinetxProcess.R")
source("kinetxProcessCalcium.R")
source("kinetxSummary.R")

#print("processing functional data")
kinetxProcess('../data/','../output/','FL1')
kinetxProcess('../data/','../output/','FL4')

#print("summarising data")
kinetxSummary('../output/summaryFL1FileV1.csv','../output/','FL1')
kinetxSummary('../output/summaryFL4FileV1.csv','../output/','FL4')

#print("processing calcium data")
kinetxProcessCalcium('../data/','../output/','FL1','ALL')

