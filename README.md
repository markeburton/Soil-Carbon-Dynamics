# Soil-Carbon-Dynamics
CarbonSolve 

This repository includes R code for SNAP (Ritchie 2014) and SNPGRAZE (Ritchie 2020) soil carbon models with additional versions allowing for updated functionality, customizability, and  additional calculations. 

SNAP v1.0 is a mass balance model originally published in Ritchie (2014) using field data collected in Serengeti National Park, Tanzania to predict soil carbon dynamics in a tropical grassland in response to variations in grazing intensity. It examines the response of leaf area index (LAI) to a gradient of grazing intensities. SNAP v1.0 includes five input variables: mean annual precipitation, fire frequency, percent sand, mean proportion of plant as lignin and cellulose, and grazing intensity. 

SNAPGRAZE v.1.0 builds from SNAP by explicitly accounting for management options. This new aspect was addded by incorporating an Episodic Herbivory Model (EHM) tracking plant growth prior to, during, and following grazing events. SNAPGRAZE v1.0 introduces several new input variables: mean annual temperature, a correction factor for systems dominated by annuals or shrubs (APC), maximim relative growth rate (rgr), soil sampling depth, and a series of grazing pattern parameters (period of stay, time prior to grazing episode, number of animals and pastures in the grazing system, mean mass of the grazing animals, and an optional value to specify the default daily consumption per animal. 

## Documentation of model parameter, inputs, and calibration changes in models

The currently available code is described below, with additional models scheduled to be added. 

### SNAPGRAZE v2.3
1. The original SNAPGRAZE model has been updated with new parameter coefficients and allows for calibration to fit a specific project area.
  New updates to SNAPGRAZE v1.0 include inputs for the dominant vegetation class (bare, annuals, perennials, mixture, or forbs) to replace the plant lignin and cellulose percentage and       APC. Similiarly, the number of rotations in a grazing system is an input that replaces the period of stay and time prior to a grazing episode.  The SNAPGRAZE v2.3 also allows for           calibrations to be used to update parameter coeffiecents related to microbial respiration and acvtivity (relateded to days with sufficient soil moisture). 
   
2. Additionally, we added a monte carlo simulation of the project scenario to better understand uncertainty to SOC predictions.



## Resources

The code in this repository is built from the original SNAP and SNAPGRAZE models published in:

Ritchie ME (2014) Plant compensation to grazing and soil carbon dynamics in a tropical grassland. PeerJ. 2:e233 

Ritchie ME (2020) Grazing management, forage production and soil carbon dynamics. Resources. 9(4):49


## Contact

Mark Burton

Email: mark.burton@briwildlife.org

Biodiversity Research Institute; CarbonSolve


Mark Ritchie 

Email: markritchie@soilsfuture.com

Soils for the Future; CarbonSolve


## License

GNU GENERAL PUBLIC LICENSE v3

