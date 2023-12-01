# Fall23D2P
Repository to hold code for Fall23 Data to policy presentation

* Note: pop.rda was manually entered from the source: https://cohealthviz.dphe.state.co.us/t/HealthInformaticsPublic/views/ColoradoPopulationEstimates/PopulationEstimates?iframeSizedToWindow=true&%2C%3Aembed=y&%2C%3AshowAppBanner=false&%2C%3Adisplay_count=no&%2C%3AshowVizHome=no

file descriptions
* jobLossVisualsion.R: File used to create a rate map as well as creating COMap.rda
* name.rda: holds an array that was used to help construct pop.rda. Not currently used but kept due to convenience
* ScanResults.txt: Stores the output of spacial scan test calls
* CEPPResults.txt: Stores the output of CEPP test calls
    * Note: This is not used in the presentation, but provides a similar conclusion as BN
* BNResults.txt: Stores the output of BN test calls
analysis.R: File that performs the analysis and generates all plots made. The two in the presentation were generated from
    1. BN with c* = 5000
    2. Spatial Scan

