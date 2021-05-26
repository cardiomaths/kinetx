# kinetx - R code to process real-time flow cytometry data

This code is associated with the paper 

Pre-requisites

The following R packages are pre-requisites for running kinetx:

minpack

psych

kinetx dataset:
10 raw datafiles produced by the flow cytometry assay described in the above assay are included under 'data'. These allow the code to be tested.

Summary file format:

Format of functional data summary file:

donor - unique donor number taken from raw datafile directory name

output - function

measure - FL1 or FL4

agonist - 

concentration -

drug -

0Secs - level at 10 seconds

10Secs - level at 10 seconds

20Secs - level at 20 seconds

30Secs - level at 30 seconds

40Secs - level at 40 seconds

60Secs - level at 60 seconds

90Secs - level at 90 seconds

120Secs - level at 120 seconds

180Secs - level at 180 seconds

10SecRoc - rate of change between 0 and 10 seconds

20SecRoc - rate of change between 10 and 20 seconds

30SecRoc - rate of change between 20 and 30 seconds

40SecRoc - rate of change between 30 and 40 seconds

60SecRoc - rate of change between 40 and 60 seconds

90SecRoc - rate of change between 60 and 90 seconds

120SecRoc - rate of change between 90 and 120 seconds

180SecRoc - rate of change between 120 and 180 seconds

IsotypeMedian - median level of isotype

zeroPoint - level at zero

maxRateChange -

maxRateChangePosition -

minRateChange -

minRateChangePosition -

maxRateAccel - maximum rate of acceleration

maxRateAccelPosition -

Label -

FL1: pred120 < 4000 low responder 
FL4: pred120 < 1000 low responder 

LabelMetric: decreasing 0.9 < linear < 2 increasing

LabelMetric - Roc120/(Roc10 + Roc20)/2 

RocLabel -

FL1 RoC: slow 40 < medium < 80 fast

FL4 RoC: slow 5 < medium < 10 fast

Roc - Rate of change between 0 and 120 seconds

EarlyRateLabel - 

FL1: slow 40 < medium < 80 fast

FL4: slow 4 < medium < 12 fast

EarlyRateMetric - Rate of change in first 20 seconds
