# assembloid_LFP_Zourrayetal2025
These MATLAB functions perform the analysis of assembloid local field potentials performed for the paper Zourray et al. 2025.


This set of functions perform the LFP analyses shown in Figure 1 of Zourray et al 2025.  They are performed on local field potential recordings from >d200 brain assembloids.  For more information about the assembloids and recording protocol, please see the Methods section of this paper.

To generate the analyses from scratch, run analyseOrganoidLFP(xlsfile,params) where xlsfile is the absolute filepath of an Excel sheet containing metadata about each recording, and params gives modifying parameters for the function.  The function will create and save a struct Organoid for each organoid recording (row) listed in the xls, alongside a summary figure of the activity through the recording.

To generate similar figures to the paper, run organoidAnalysisReport55mM_20241005 with xlsfile defined similarly.  This will iterate through and load the pregenerated Organoid structures, extract relevant data, and plot these on raincloud plots.  This function generates similar plots to the paper, alongside other common parameters used to assess activity and signal spectral structure.  Note that the figures in the paper used the values generated by this function, but plots were made in GraphPad for consistency with the remainder of the paper.
