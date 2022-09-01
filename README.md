# ArchiveSize
Paper "Effects of Archive Size on Computation Time and Solution Quality for Multi-Objective Optimization"

Here the explanation of the important files are listed:

**generateData.m**: run NSGA-II, MOEA/D-PBI and NSGA-III on DTLZ1-4 and their minus versions, and save offspring and population.

**LazyHVSelection.m**: lazy greedy inclusion hypervolume subset selection (LGI-HSS)

**DSS.m**: greedy distance-based inclusion

**CDR.m**: greedy crowding distance-based removal

"NDSort.m": T-ENS is used to remove dominated solutions in the archive 

**ArchiveStrategy.m**: contain various archiving strategies mentioned in the paper

**ParArcExp.m**: maintain the archive using various strategies and output the final archive

**ParFinalSelection.m**: select the final solution set from the final archive

**ParDRSRemove.m**: remove some of DRSs in the final archive and select the final solution set from it
