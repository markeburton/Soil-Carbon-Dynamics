

#### SNAPGRAZE 2.3 ####

# script metadata ==========================================================================================
# Description: 
    # This version of SNAPGRAZE is based on the Rithie et al. (2020) manuscript. 
    # Citation: Ritchie ME (2020) Grazing management, forage production and soil carbon dynamics. Resources 9(4): 49
# Script Author: Mark Burton
# Script Date (last update): 2025-04-11
# ====================================================================================================

# Version 2.3 has several important changes to SNAPGRAZE v2.0
  # Note that previous versions of SNAPGRAZE include a correction factor of 0.9 for ANPP max to better reflect experimental data from the Serengeti where the model was developed. 
    # This correction is no longer the default. However, the model can be calibrated using a correction factor (ANPP.Correction) to best reflect the applicable system. 
  # The version below uses calibrated parameter coefficients from current project data as described below:

# Calibration Information (other) - Some parameter coefficients can also be calibrated to fit project specific data.
# These are not calculated in this script but input from a separate analysis.


# * Coeff.Resp.Temp: Default is 0.94198
Coeff.Resp.Temp <-	0.91  # From Calibration tab is Carbon Calculation Spreadsheet.
# * Coeff.wetsoil: default is 0.00043
Coeff.wetsoil   <- 0.00051 # From Calibration tab in Carbon Calculation Spreadsheet.
# * Coeff_Respiration: default is 0.00044
Coeff_Respiration <- 0.0005 # From Calibration tab in Carbon Calculation Spreadsheet.
# RGRmax - This is a required input in the input file, so changing here will not impact calculations
# * 0.033 used here but default is 0.036



##### Import input parameters.
#First, read in excel template with needed inputs.
SnapGrazeInputs <- read_excel(paste0(source_dir, "CarbonModellingInputs.xlsx"), sheet = "SNAPGRAZE_2_3_Inputs")
#Lets remove any that do not have VegClass or SOC measurements data right off the bat
SnapGrazeInputs <- SnapGrazeInputs %>%
  filter(!is.na(VegClass) & VegClass != "" &
           !is.na(SOC.measured) & SOC.measured != "")
dim(SnapGrazeInputs)

# Now we start building out the model and intermediate equations.
SNAPGRAZEv2.3 = function(SnapGrazeInputs) {
  
  #Now, lets define the parameters for use in the model.
  YEARS                 <-  SnapGrazeInputs$"Years"[1] + 10        # How long will be project last in years [40 is default] -- Added an extra 10 years beyond project to look at kback calculation 
  SOC.measured          <-  SnapGrazeInputs$"SOC.measured"         # SOC data for the points - model will accrue SOC using this as the starting point under Project scenarios
  A                     <-  SnapGrazeInputs$"Project Area"         # Size of the area that will be in the project.
  # Environmental
  FIRE.base             <-  SnapGrazeInputs$"FIRE.base"            # This is a fire frequency (number of fires per year (0-1)) in the baseline scenario
  FIRE.project          <-  SnapGrazeInputs$"FIRE.project"         # This allows for a reduction in fire freq based on project scenarios.
  MAT                   <-  SnapGrazeInputs$"MAT"                  # Mean Annual Temperature
  MAP                   <-  SnapGrazeInputs$"MAP"                  # Mean Annual Precipitation
  SAND                  <-  SnapGrazeInputs$"SAND"                 # % of soil that is sand
  # Vegetation
  VEGCLASS              <-  SnapGrazeInputs$"VegClass"             # Selection of possibilities for dominant vegetation class as baseline.
  r                     <-  SnapGrazeInputs$"r"                    # Maximum relative growth rate
  # Sampling
  DEPTH                 <-  SnapGrazeInputs$"DEPTH"                # Depth we are using for calculations 
  # Grazing Information
  W                     <-  SnapGrazeInputs$"W" 
  N.base                <-  SnapGrazeInputs$"Animals.base"         # Number of animals in the baseline scenario
  n.base                <-  SnapGrazeInputs$"Pasture.base"         # Number of pastures in the baseline scenario
  ROTATIONS.base        <-  SnapGrazeInputs$"Rotations.base"       # Number of rotations in the baseline scenario
  N.project             <-  SnapGrazeInputs$"Animals.project"      # Number of animals in the project scenario
  n.project             <-  SnapGrazeInputs$"Pasture.project"      # Number of pastures in the project scenario
  ROTATIONS.project     <-  SnapGrazeInputs$"Rotations.project"    # Number of rotations in the project scenario
  
  #### First, use the BASELINE DATA to predict SOC for each point
  # Starts with Episodic Herbivory Model (EHM)
  # Equation 15 from Richie et al. 2020 - Calculate annual aboveground production in absence of grazing (ANPPmax)
  # 0.00834 kJ mol/k is the gas constant
  ANPPmax =  exp(12.039 + (0.718*log(MAP)) - (25.18 / (0.00834 * (273.15 + MAT)))) * (1.33 - (0.0075 * SAND))
  
  # Biomass in the absence of grazing (Sk)  g^m2 year^-1
  # 0.9 is a slight correction is better match the data from Serengeti - it can calibrated to better match exclosure experiments for a given system
  # Published version of SNAPGRAZE used a 0.9 correction factor - This is removed in the default but you can calibrate to better match exclosure experiments for a given system
  ANPP.Correction = 1.0  #enter 0.9 to match published version. 
  Sk = (ANPPmax/ANPP.Correction)
  # Length of plant growing season
  G = (22.993 * MAT) - (Coeff.Resp.Temp * MAT^2) + (0.073 * MAP)
  # Per Animal Daily Consumption
  Cg =  2 * (5300 + 770*log(W))
  # Biomass at the onset of the growing season, produced from resource reserves
  S0 = Sk * 0.1
  # Stocking density (animals per ha)
  d = (N.base/A) 
  # L0
  L0 = (Cg /2) * (365 - G) * d * 10^(-4)
  # Period of stay
  D  = G / (n.base * ROTATIONS.base)
  D2 = ifelse(ROTATIONS.base >= 2,  D, 0)
  D3 = ifelse(ROTATIONS.base >= 3, D2, 0)
  D4 = ifelse(ROTATIONS.base >= 4, D3, 0)
  # Time prior to grazing episode
  E = (G-ROTATIONS.base*D)/(2*ROTATIONS.base)
  # Days from end of grazing period of stay to end of the growing season
  Fd = ifelse(ROTATIONS.base >= 2, 2*E, E)
  F2 = ifelse(ROTATIONS.base == 2, Fd/ROTATIONS.base, ifelse(ROTATIONS.base >= 3, Fd, 0))
  F3 = ifelse(ROTATIONS.base == 3, E, ifelse(ROTATIONS.base > 3, 2*E, 0))
  F4 = ifelse(ROTATIONS.base >= 4, G-E-D-F1-D2-F2-D3-F3-D4,0)

  # Biomass at start of grazing episode
  Se = (Sk * S0) / ((Sk * exp(-r * E)) + S0 * (1 - exp(-r * E)))
  # Relative loss rate of biomass to grazing 
  g = (d * Cg * n.base * 10^(-4)) / Se
  # Biomass removed during grazing period of stay.
  Lg = D * d * Cg * n.base * 10^(-4)
  # Biomass at end of grazing episode
  Sg = ((Sk * Se) / ((Sk * exp((-r + g ) * D)) + (Se * (1-exp((-r + g )*D)))))
  # Biomass at end of growing season
  Sf = (Sk * Sg) / ((Sk * exp(-(r * Fd))) + Sg * (1 - exp(-(r * Fd))))
  #
  Se2 = (Sk*Sg) / (Sk * exp((-r)*Fd) +Sg  *(1 - exp((-r) * Fd)))
  Sg2 = (Sk*Se2)/ (Sk*exp((-r+g)*D2) +Se2 *(1 - exp((-r + g) * D2)))
  Se3 = (Sk*Sg2)/ (Sk*exp((-r) * F2)+Sg2 *(1 - exp((-r) * F2)))
  Sg3 = (Sk*Se3)/ (Sk*exp((-r+g)*D3) +Se3 *(1 - exp((-r + g) * D3)))
  Se4 = (Sk*Sg3)/ (Sk*exp((-r)  *F3)+Sg3 *(1 - exp((-r) * F3)))
  Sg4 = (Sk*Se4)/ (Sk*exp((-r+g)*D4) +Se4 *(1 - exp((-r + g) * D4)))
  # Consumption calculations (including off season consumption)
  CONSUMED = Cg*d*(D)*10^-4
  CONSUMED.OFF = ifelse(L0 < Se2, L0, Se2)
  
  # ANPP intermediate calculation
  ProdG = ifelse(Sg>Se, Sf-S0, Se-S0+Se2-Sg+Se3-Sg2+Se4-Sg3+Sf-Sg4)
  
  # We now need to calculate a new annual correction factor for BNPP
  BNPP.Correction = ifelse(VEGCLASS == "ANNUALS", 0.291179866899812, ifelse(VEGCLASS == "BARE", 0.291179866899812, ifelse(VEGCLASS == "MIXTURE", 0.645589933449906, 1)))

  # Annual aboveground production under grazing
  ANPPest = ProdG * BNPP.Correction
  
  # Other Productivity Calculations
  # Annual belowground productivity under grazing
  # Increase in biomass from initial amount (S0) to a biomass Sk after G days
  Pu = Sk - S0
  
  # Equation 16 - Annual below ground productivity under grazing
  BNPPest = ((0.602 * MAP) - (0.00038 * MAP^2) + (5.888 * MAT)) * (ANPPest/ANPPmax) * BNPP.Correction
  
  # Soil Organic Carbon (SOC) calculations.
  # Calculate mean proportion of plant as lignin and cellulose using dominant vegetation class.
  LIGCELL = ifelse(VEGCLASS == "ANNUALS", 23.248751199377, ifelse(VEGCLASS == "BARE", 23.248751199377, ifelse(VEGCLASS == "MIXTURE", 29.4626011324048, 35.6764510654327)))
  # Plant-derived soil carbon input in year y
  PDSOC.base = (0.45 * (((LIGCELL/100) * ANPPest * (1 - FIRE.base)) + (((LIGCELL/100) + 0.05) * BNPPest))) * (-0.3559+0.3914 * log(DEPTH))
  # Dung-derived soil carbon input in year y
  DDSOC.base = ((LIGCELL/100) * (CONSUMED + CONSUMED.OFF) * 0.45 * (0.4718 - 0.0009*(D))) * (-0.3559+0.3914* log(DEPTH))
  
  # WETDAYS (eq 19) 
  WETDAYS = (Coeff.wetsoil * MAP - 0.025) * G
 
  # SOCeq
  SOCeq.large.base  = (((PDSOC.base+DDSOC.base) / (WETDAYS * (0.7 + 0.3 *(SAND/100)) * Coeff_Respiration)) + (0.579/Coeff_Respiration)) / 100    
  SOCeq.small.base  = ((PDSOC.base+DDSOC.base) / (WETDAYS * (0.7+0.3*(SAND/100)) * exp(-10.872)))^(1/1.296) / 100  
  SOCeq.base        =  pmin(SOCeq.large.base, SOCeq.small.base, na.rm = TRUE)
  
  

  ####
  #### Next, use PROJECT parameters to predict SOC for each point 
  ####
  # Per Animal Daily Consumption
  Cg = ifelse(is.numeric(Cg), Cg, 2 * (5300 + 770*log(W)))
  Cg = ifelse(is.numeric(Cg), Cg, 2 * (5300 + 770*log(W)))
  #Stocking density (animals per ha)
  d = (N.project/A) 
  # L0
  L0 = (Cg /2) * (365 - G) * d * 10^(-4)
  # Period of stay
  D  = G / (n.project * ROTATIONS.project)
  D2 = ifelse(ROTATIONS.project >= 2,  D, 0)
  D3 = ifelse(ROTATIONS.project >= 3, D2, 0)
  D4 = ifelse(ROTATIONS.project >= 4, D3, 0)
  # Time Prior to grazing episode
  E = (G-ROTATIONS.project*D)/(2*ROTATIONS.project)
  # Days from end of grazing period of stay to end of the growing season
  Fd = ifelse(ROTATIONS.project >= 2, 2*E, E)
  F2 = ifelse(ROTATIONS.project == 2, Fd/ROTATIONS.project, ifelse(ROTATIONS.project >= 3, Fd, 0))
  F3 = ifelse(ROTATIONS.project == 3, E, ifelse(ROTATIONS.project > 3, 2*E, 0))
  F4 = ifelse(ROTATIONS.project >= 4, G-E-D-F1-D2-F2-D3-F3-D4,0)
  
  # Biomass at start of grazing episode
  Se = (Sk * S0) / ((Sk * exp(-r * E)) + S0 * (1 - exp(-r * E)))
  # Relative loss rate of biomass to grazing 
  g = (d * Cg * n.project * 10^(-4)) / Se
  # Biomass removed during grazing period of stay.
  Lg = D * d * Cg * n.project * 10^(-4)
  # Biomass at end of grazing episode
  Sg = ((Sk * Se) / ((Sk * exp((-r + g ) * D)) + (Se * (1-exp((-r + g )*D)))))
  # Biomass at end of the growing season
  Sf = (Sk * Sg) / ((Sk * exp(-(r * Fd))) + Sg * (1 - exp(-(r * Fd))))
  #
  Se2 = (Sk*Sg) / (Sk * exp((-r)*Fd) +Sg  *(1 - exp((-r) * Fd)))
  Sg2 = (Sk*Se2)/ (Sk*exp((-r+g)*D2) +Se2 *(1 - exp((-r + g) * D2)))
  Se3 = (Sk*Sg2)/ (Sk*exp((-r) * F2)+Sg2 *(1 - exp((-r) * F2)))
  Sg3 = (Sk*Se3)/ (Sk*exp((-r+g)*D3) +Se3 *(1 - exp((-r + g) * D3)))
  Se4 = (Sk*Sg3)/ (Sk*exp((-r)  *F3)+Sg3 *(1 - exp((-r) * F3)))
  Sg4 = (Sk*Se4)/ (Sk*exp((-r+g)*D4) +Se4 *(1 - exp((-r + g) * D4)))
  # Consumption calculations
  CONSUMED = Cg*d*(D)*10^-4
  CONSUMED.OFF = ifelse(L0 < Se2, L0, Se2)
  
  #ANPP intermediate calculation
  ProdG = ifelse(Sg>Se, Sf-S0, Se-S0+Se2-Sg+Se3-Sg2+Se4-Sg3+Sf-Sg4)
  
  # For project veg class, we assume a one tier upgrade over baseline scenario
  VEGCLASS.project <- ifelse(VEGCLASS == "BARE", "MIXTURE",ifelse(VEGCLASS == "ANNUALS", "MIXTURE",ifelse(VEGCLASS == "MIXTURE", "PERENNIALS",ifelse(VEGCLASS == "PERENNIALS", "PERENNIALS", "FLAG")))) 
  # We now need to calculate a new annual correction factor for BNPP
  BNPP.Correction = ifelse(VEGCLASS.project == "ANNUALS",0.291179866899812, ifelse(VEGCLASS.project == "BARE",0.291179866899812, ifelse(VEGCLASS.project == "MIXTURE", 0.645589933449906, 1)))
  
  # Annual aboveground production under grazing
  ANPPest = ProdG * BNPP.Correction
  
  # Other Productivity Calculations
  # Annual belowground productivity under grazing
  # Increase in biomass from initial amount (S0) to a biomass Sk after G days
  Pu = Sk - S0
  # Equation 16 - Annual below ground productivity under grazing
  BNPPest = ((0.602 * MAP) - (0.00038 * MAP^2) + (5.888 * MAT)) * (ANPPest/ANPPmax) * BNPP.Correction
  
  # Soil Organic Carbon (SOC) calculations.
  # Calculate Mean proportion of plant as lignin and cellulose using dominant vegetation class.
  LIGCELL = ifelse(VEGCLASS.project == "ANNUALS",23.248751199377, ifelse(VEGCLASS.project == "BARE",23.248751199377, ifelse(VEGCLASS.project == "MIXTURE", 29.4626011324048, 35.6764510654327)))
  # Plant-derived soil carbon input in year y
  PDSOC.project = (0.45 * (((LIGCELL/100) * ANPPest * (1 - FIRE.project)) + (((LIGCELL/100) + 0.05) * BNPPest))) * (-0.3559 + 0.3914 * log(DEPTH))
  # Dung-derived soil carbon input in year y
  DDSOC.project = ((LIGCELL/100) * (CONSUMED + CONSUMED.OFF) * 0.45 * (0.4718 - 0.0009*(D))) * (-0.3559 + 0.3914 * log(DEPTH))
  
  # SOCeq
  SOCeq.large.proj  = (((PDSOC.project+DDSOC.project) / (WETDAYS * (0.7 + 0.3 *(SAND/100)) * Coeff_Respiration)) + (0.579/Coeff_Respiration)) / 100 
  SOCeq.small.proj  = ((PDSOC.project+DDSOC.project) / (WETDAYS * (0.7+0.3*(SAND/100)) * exp(-10.872)))^(1/1.296) / 100   
  SOCeq.project     = pmin(SOCeq.large.proj, SOCeq.small.proj, na.rm = TRUE)
  

  # Now lets use the baseline soil equilibrium data and loop through creating calculations for each year of the project.  
  # Create an empty dataframe to store the results with a new column for each year. 
  #Add the SOCeq into the dataframe as an intermediate step
  SnapGrazeOutputs <- SnapGrazeInputs %>%
    mutate(
      WETDAYS = WETDAYS,
      SAND = SAND,
      DEPTH = DEPTH,
      LIGCELL = LIGCELL,
      CONSUMED = CONSUMED,
      CONSUMED.OFF = CONSUMED.OFF,
      ANPPest = ANPPest,
      BNPP.Correction = BNPP.Correction,
      BNPPest = BNPPest,
      SOCeq.base = SOCeq.base,
      PDSOC.y = PDSOC.project,
      DDSOC.y = DDSOC.project,
      SOCeq.large.proj = SOCeq.large.proj,
      SOCeq.small.proj = SOCeq.small.proj,
      SOCeq.project = SOCeq.project
    )
  
  #This creates a seq to account for years, 40 is the default project length but we expand it out to 50 to calculate the kback below.
  year.seq <- paste0("SOC_Year", 1:YEARS)
  
  # Create new columns for each year and initialize with NA
  for (i in 1:YEARS) {
    SnapGrazeOutputs[[paste0("SOC_Year", i)]] <- NA
    SnapGrazeOutputs[[paste0("dSOC_Year", i)]] <- NA
  }
  
  # Loop through each row
  for (row_index in 1:nrow(SnapGrazeOutputs)) {
    row_values <- SnapGrazeOutputs[row_index, ]
    
    # Loop through each year
    for (year in 0:YEARS) {
      if (year == 0) {
        # Initialize the value for the first year based on SOC.measured
        SnapGrazeOutputs[row_index, paste0("SOC_Year", year)] <- row_values["SOC.measured"]
      } else {
        # Calculate the value based on the previous year's value and other parameters
        prev_year_value <- SnapGrazeOutputs[row_index, paste0("SOC_Year", year - 1)]
        new_SOC <- prev_year_value +
                        ((row_values["PDSOC.y"] + row_values["DDSOC.y"] - pmax(
                                (row_values["WETDAYS"] * (0.7 + 0.3 * (row_values["SAND"] / 100)) * (-0.579 + Coeff_Respiration * (prev_year_value * 100))),
                                (row_values["WETDAYS"] * (0.7 + 0.3 * (row_values["SAND"] / 100)) * (exp(-10.872) * (prev_year_value * 100)^(1.296))),
                                 na.rm = TRUE ) ) / 100)
      ##### CHECK 0.579 #####
        
        # Assign the new value to the dataframe for the current year and row
        SnapGrazeOutputs[row_index, paste0("SOC_Year", year)] <- new_SOC
      }
    }
  }
  
  # Loop through each row
  for (row_index in 1:nrow(SnapGrazeOutputs)) {
    row_values <- SnapGrazeOutputs[row_index, ]
    
    # Loop through each year
    for (year in 0:YEARS) {
      if (year == 0) {
        # Initialize the value for the first year based on SOC.measured
        SnapGrazeOutputs[row_index, paste0("dSOC_Year", year)] <- row_values["SOC.measured"]
      } else {
        # Calculate the value based on the previous year's value and other parameters
        prev_year_value <- SnapGrazeOutputs[row_index, paste0("SOC_Year", year - 1)]
        new_dSOC <- (row_values["PDSOC.y"] + row_values["DDSOC.y"] - pmax(
                        (row_values["WETDAYS"] * (0.7 + 0.3 * (row_values["SAND"] / 100)) * (-0.579 + Coeff_Respiration * (prev_year_value * 100))),
                        (row_values["WETDAYS"] * (0.7 + 0.3 * (row_values["SAND"] / 100)) * (exp(-10.872) * (prev_year_value * 100)^(1.296))),
                              na.rm = TRUE ) ) / 100
        ##### CHECK 0.579 #####
        # Assign the new value to the dataframe for the current year and row
        SnapGrazeOutputs[row_index, paste0("dSOC_Year", year)] <- new_dSOC
      }
    }
  }
  
  #Lets also tack on some mean values for dSOC to the df 
  # Find column indices for columns starting with "dSOC" but we dont need the initital columns because they are already included.
  SnapGrazeOutputs <- subset(SnapGrazeOutputs, select = -c(dSOC_Year0,SOC_Year0))

  # Dynamically calculate the mean dSOC by row, we extended the years out beyond the 40 year project span to calculate kback, but we only want to take the mean dSOC for the length of the project (default is 40)
  #  dSOC_Year1 through dSOC_Year[max_year]
  dSOC_columns <- grep(paste0("^dSOC_Year([1-9]$|[1-9][0-9]$|", ifelse(SnapGrazeInputs$Years[1] >= 10, paste0("1-", SnapGrazeInputs$Years[1]), SnapGrazeInputs$Years[1]), "$)"), 
                       names(SnapGrazeOutputs))
  dSOC_columns <- dSOC_columns[names(SnapGrazeOutputs)[dSOC_columns] %in% paste0("dSOC_Year", 1:SnapGrazeInputs$Years[1])]
  # Calculate row-wise average for dSOC columns
  SnapGrazeOutputs$Mean_dSOC <- rowMeans(SnapGrazeOutputs[, dSOC_columns], na.rm = TRUE)
   
  # New Technique of using half life calculations
  SnapGrazeOutputs <- SnapGrazeOutputs %>%
    mutate( Mean_CO2e         = Mean_dSOC  * (44/12),  # This uses the mean annual dSOC from the first 40 years
            k                 = log(dSOC_Year1  / dSOC_Year40) / 40, 
            kback             = log(dSOC_Year40 / dSOC_Year50) / 10 , 
            SOC.consbase      = SOC.measured + dSOC_Year50 * exp(kback * 10),
            Dyears            = 2 * log(1-(k * 0.5)*(SOCeq.project - SOC.consbase) / dSOC_Year1) / -k,
            CO2e_removals     = (44/12) * (SOCeq.project - SOC.consbase) / Dyears,  # This is annual removals calculated using the 1/2 life approach
            Mean_CO2e         = Mean_dSOC  * (44/12)  # This uses the mean annual dSOC from the first 40 years
            )

  # Print the updated dataframe with calculated values for each year
  print(SnapGrazeOutputs)
  # Convert SOC vector to dataframe and transpose it
  return(SnapGrazeOutputs)
}
SNAPGRAZEv23.df <- SNAPGRAZEv2.3(SnapGrazeInputs)
names(SNAPGRAZEv23.df)

# Export as an excel workbook, so different sheets can be added with results tables.
SNAPGRAZE.wb <- createWorkbook()
addWorksheet(SNAPGRAZE.wb, "SNAPGRAZE_data")
writeData(SNAPGRAZE.wb, "SNAPGRAZE_data", SNAPGRAZEv23.df)
saveWorkbook(SNAPGRAZE.wb, paste0(out_dir, "SNAPGRAZEv23_", Project_initials, "_", date, ".xlsx"), overwrite = TRUE)


#### * Summary Table #### 
Removals.mean <- mean(SNAPGRAZEv23.df$Mean_CO2e, na.rm = T)
Removals.se <- ifelse(is.na(se(SNAPGRAZEv23.df$Mean_CO2e, na.rm = T)), "Sample size = 1", se(SNAPGRAZEv23.df$Mean_CO2e, na.rm = T))
Removals.Project <- prettyNum(Removals.mean * SNAPGRAZEv23.df$`Project Area`[1],  big.mark = ",", scientific = FALSE )

Statistic        <- c('Removals.mean', 'Removals.se', 'Removals.Project')
Value.PROJECT   <- c(Removals.mean, Removals.se, Removals.Project)
SNAPGRAZEv23_Removals.df <- data.frame(Statistic = Statistic, Values = Value.PROJECT)

# Lets also create an stratified summary statistics table by community
SNAPGRAZEv23_CommunityRemovals.df <- SNAPGRAZEv23.df %>%
  group_by(Community) %>%
  summarise(
    n                   = sum(!is.na(CO2e_removals)),  # Number of sites in each strata that are not NA
    mean_CO2e_removals  = mean(CO2e_removals, na.rm=TRUE),  # Mean removals tonnes CO2e ha-1 yr-1 by strata
    sd_CO2e_removals    = sd(CO2e_removals, na.rm=TRUE),        # sd removals tonnes CO2e ha-1 yr-1 by strata
    Uncertainty.perc    = ((qnorm(0.95)*sd_CO2e_removals/(n^0.5))/mean_CO2e_removals)*100,
    CommArea_ha         = mean(`Project Area`, na.rm=TRUE),
    TotalRemovals_tCO2e = CommArea_ha * mean_CO2e_removals, 
    Uncertainty.Weightsq=((Uncertainty.perc*TotalRemovals_tCO2e)^2),
    .groups = "drop"
  ) %>%
# also add a total row at the bottom of the table. 
    bind_rows(
    summarise(., 
              Community = "Total",
              n = sum(n),
              mean_CO2e_removals = NA,
              sd_CO2e_removals = NA,
              CommArea_ha = sum(CommArea_ha, na.rm = TRUE),
              TotalRemovals_tCO2e = sum(TotalRemovals_tCO2e, na.rm = TRUE),
              Uncertainty.perc = (sqrt(sum(Uncertainty.Weightsq))/TotalRemovals_tCO2e),
    )     
  )%>%
  select(Community,n,mean_CO2e_removals,sd_CO2e_removals,Uncertainty.perc,CommArea_ha,TotalRemovals_tCO2e)
SNAPGRAZEv23_CommunityRemovals.df

#Now lets write the summary statistics table for export. 
addWorksheet(SNAPGRAZE.wb, "Removals_StrataSummary")
writeData(SNAPGRAZE.wb, "Removals_StrataSummary", SNAPGRAZEv23_CommunityRemovals.df)
saveWorkbook(SNAPGRAZE.wb, paste0(out_dir, "SNAPGRAZEv23_", Project_initials, "_", date, ".xlsx"), overwrite = TRUE)



#### MONTE CARLO ####

# Confirm input file read in at beginning of script
dim(SnapGrazeInputs) 

#Set the number of iterations for the Monte Carlo analysis here
num_iterations <- 100
    #Check the number of rows, the output of the monte carlo iterations is set to print 100 rows, it can be change below, but will get large quickly
    length(unique(SnapGrazeInputs$Site)) * num_iterations


# use set seed so Monte Carlo analysis will be consistent from run to run if desired
set.seed(1234)

# First, lets rename a few columns from the input file 
SnapGrazeInputs <- SnapGrazeInputs %>%
  rename(
    A                  = "Project Area",      # Size of the area that will be in the project.
    VEGCLASS           = "VegClass",          # Dropdown selection of possibilities for dominant vegetation class as baseline
    N.base             = "Animals.base",      # Number of animals in the baseline scenario
    n.base             = "Pasture.base",      # Number of pastures in the baseline scenario
    ROTATIONS.base     = "Rotations.base",    # Number of rotations in the baseline scenario
    ROTATIONS.project  = "Rotations.project", # Number of rotations in the project scenario
    N.project          = "Animals.project",   # Number of animals in the project scenario
    n.project          = "Pasture.project"    # Number of pastures in the project scenario
  )

# We will also need the calculate the mean and se for temperature and precipitation for each strata.
Clim.df <- SnapGrazeInputs %>%
  group_by(Community)  %>%
  summarise(
    MAP.strata.mean = mean(MAP, na.rm=TRUE),
    MAP.strata.se = se(MAP, na.rm=TRUE),
    MAT.strata.mean = mean(MAT, na.rm=TRUE),
    MAT.strata.se = se(MAT, na.rm=TRUE)
  )

# Initialize a results list
results_list <- list()

# Loop over each site (row in the input file)
for (site in 1:nrow(SnapGrazeInputs)) {
  
  # Extract row-specific input values pulled directly from the file
  site_data <- SnapGrazeInputs[site, ]
  
  # Get the actual site name
  site_name <- site_data$Site
  Community_name <- site_data$Community 
  SOC.measured <- site_data$SOC.measured
  
  # Create a df to store site results 
  site_results <- data.frame(
    Site      = rep(site_name, num_iterations),
    Community = rep(Community_name, num_iterations),
    Iteration = 1:num_iterations,
    SOC.measured =  rep(SOC.measured, num_iterations)
  )
  
  # Monte Carlo iterations
  for (i in 1:num_iterations) {
    
    # Generate random parameter coefficients 
    # The coefficients for the following parameters are taken from `SNAPGRAZE Values Parameters Coefficients Uncertainty` 
    # The first number is coefficient, the second is the uncertainty (SE was used in this case because we care about variation around the mean) 
    r                 <- rnorm(1, 0.033, 0.004)    #  From Calibration tab is Carbon Calculation Spreadsheet. Default: rnorm(1, 0.036, 0.004)
    G.MAT.B0          <- rnorm(1, 22.993, 1.727)
    G.MAT.B1          <- rnorm(1, 0.91, 0.099)     # From Calibration tab is Carbon Calculation Spreadsheet. Default: rnorm(1, 0.942, 0.099) 
    G.MAP.B0          <- rnorm(1, 0.073, 0.014)
    Cg.B0             <- rnorm(1, 5300, 0)
    Cg.B1             <- rnorm(1, 770, 0)
    LigCell.Annual    <- rnorm(1, 23.248751199377, 1.55999049699112)
    LigCell.Mixture   <- rnorm(1, 29.4626011324048, 1.55999049699112)
    LigCell.Perennial <- rnorm(1, 35.6764510654327, 1.42118061262896)
    ANPPmax.Sand.B0   <- rnorm(1, 1.33, 0.05)
    ANPPmax.Sand.B1   <- rnorm(1, 0.0075, 0.0008)
    ANPPmax.MAP.B0    <- rnorm(1, 0.718, 0.145)
    ANPPmax.MAT.B0    <- rnorm(1, 25.18, 6.767)
    ANPPmax.MAT.B1    <- rnorm(1, 12.039, 3.533)
    APC.Annual        <- rnorm(1, 0.291179866899812, 0.155142849652921)
    APC.Mixture       <- rnorm(1, 0.645589933449906, 0.155142849652921)
    APC.Perennial     <- rnorm(1, 1, 0)
    BNPP.MAP.B0       <- rnorm(1, 0.602, 0.183)
    BNPP.MAT.B0       <- rnorm(1, 5.888, 3.123)
    BNPP.MAT.B1       <- rnorm(1, 0.00038, 0.00022)
    WETDAYS.Bo        <- rnorm(1, 0.025, 0.033)
    WETDAYS.B1        <- rnorm(1, 0.00051, 0.00006) # From Calibration tab is Carbon Calculation Spreadsheet. Default: rnorm(1, 0.00043, 0.00006)
    SOC.B0            <- rnorm(1, -10.872, 0) 
    SOC.B1            <- rnorm(1, 1.296, 0)
    SOC.B2            <- rnorm(1, 0.579, 0.45)
    SOC.B3            <- rnorm(1, 0.0005, 0.00007) # From Calibration tab is Carbon Calculation Spreadsheet. Default: rnorm(1, 0.00044,0.00007)
    
    # These calculations also generate a value for temperature and precipitation
    # First we need to get the data for MAP and MAT from the Clim.df that matches the community for the community in the current iteration
    clim_row <- Clim.df %>% filter(Community == Community_name)
    # Also add a NA as a flag in case Community data is missing from a row.
    if (nrow(clim_row) == 1) {
      MAP.strata <- rnorm(1, clim_row$MAP.strata.mean, clim_row$MAP.strata.se)
      MAT.strata <- rnorm(1, clim_row$MAT.strata.mean, clim_row$MAT.strata.se)
    } else {
      MAP.strata <- NA
      MAT.strata <- NA
    }
    
    # Now that we have generated the random parameters and pulled in the input values from the input file, lets run the model
    # Starts with Episodic Herbivory Model (EHM)
    # Equation 15 - Calculate annual aboveground production in absence of grazing (ANPPnax)
    # 0.00834 kJ mol/k is the gas constant
    ANPPmax =  exp(ANPPmax.MAT.B1 + (ANPPmax.MAP.B0*log(MAP.strata)) - (ANPPmax.MAT.B0 / (0.00834 * (273.15 + MAT.strata)))) * (ANPPmax.Sand.B0 - (ANPPmax.Sand.B1 * site_data$SAND))
    
    # Biomass in the absence of grazing (Sk) g^m2 year^-1
    # 0.9 is a slight correction is better match the data from Serengeti 
    # Published version of SNAPGRAZE used the 0.9 correction factor - This is removed in the default but you can calibrate to better match exclosure experiments for a given system
    ANPP.Correction = 1.0  #enter 0.9 if matching the published version is desired. 
    Sk = (ANPPmax/ANPP.Correction)
    # Length of plant growing season
    G = (G.MAT.B0 * MAT.strata) - (G.MAT.B1 * MAT.strata^2) + (G.MAP.B0 * MAP.strata)
    # Per Animal Daily Consumption
    Cg = 2 * (Cg.B0 + Cg.B1*log(site_data$W))
    # Biomass at the onset of the growing season, produced from resource reserves
    S0 = Sk * 0.1
    #Stocking density (animals per ha)
    d = (site_data$N.project/site_data$A) 
    # L0
    L0 = (Cg /2) * (365 - G) * d * 10^(-4)
    # Period of stay
    D  = G / (site_data$n.project * site_data$ROTATIONS.project)
    D2 = ifelse(site_data$ROTATIONS.project >= 2,  D, 0)
    D3 = ifelse(site_data$ROTATIONS.project >= 3, D2, 0)
    D4 = ifelse(site_data$ROTATIONS.project >= 4, D3, 0)
    # Time Prior to grazing episode
    E = (G-site_data$ROTATIONS.project*D)/(2*site_data$ROTATIONS.project)
    # Days from end of grazing period of stay to end of the growing season
    Fd = ifelse(site_data$ROTATIONS.project >= 2, 2*E, E)
    F2 = ifelse(site_data$ROTATIONS.project == 2, Fd/site_data$ROTATIONS.project, ifelse(site_data$ROTATIONS.project >= 3, Fd, 0))
    F3 = ifelse(site_data$ROTATIONS.project == 3, E, ifelse(site_data$ROTATIONS.project > 3, 2*E, 0))
    F4 = ifelse(site_data$ROTATIONS.project >= 4, G-E-D-F1-D2-F2-D3-F3-D4,0)
    # Biomass at start of grazing episode
    Se = (Sk * S0) / ((Sk * exp(-r * E)) + S0 * (1 - exp(-r * E)))
    # Relative loss rate of biomass to grazing 
    g = (d * Cg * site_data$n.project * 10^(-4)) / Se
    # Biomass removed during grazing period of stay.
    Lg = D * d * Cg * site_data$n.project * 10^(-4)
    # Biomass at end of grazing episode
    Sg = ((Sk * Se) / ((Sk * exp((-r + g ) * D)) + (Se* (1-exp((-r + g)*D)))))
    # Biomass at end of the growing season
    Sf = (Sk * Sg) / ((Sk * exp(-(r * Fd))) + Sg * (1 - exp(-(r * Fd))))
    #
    Se2 = (Sk*Sg) / (Sk*exp((-r)     *Fd) +  Sg *(1 - exp((-r) * Fd)))
    Sg2 = (Sk*Se2)/ (Sk*exp((-r+g)*D2) + Se2 *(1 - exp((-r + g) * D2)))
    Se3 = (Sk*Sg2)/ (Sk*exp((-r)     *F2) + Sg2 *(1 - exp((-r) * F2)))
    Sg3 = (Sk*Se3)/ (Sk*exp((-r+g)*D3) + Se3 *(1 - exp((-r + g) * D3)))
    Se4 = (Sk*Sg3)/ (Sk*exp((-r)     *F3) + Sg3 *(1 - exp((-r) * F3)))
    Sg4 = (Sk*Se4)/ (Sk*exp((-r+g)*D4) + Se4 *(1 - exp((-r + g) * D4)))
    # Consumption calculations
    CONSUMED = Cg*d*(D)*10^-4
    CONSUMED.OFF = ifelse(L0 < Se2, L0, Se2) # offseason consumption
    # ANPP intermediate calculation
    ProdG = ifelse(Sg>Se, Sf-S0, Se-S0+Se2-Sg+Se3-Sg2+Se4-Sg3+Sf-Sg4)
    
    # For project veg class, we assume a one tier upgrade over baseline scenario
    VEGCLASS.project <- ifelse(site_data$VEGCLASS == "BARE", "MIXTURE",ifelse(site_data$VEGCLASS == "ANNUALS", "MIXTURE",ifelse(site_data$VEGCLASS == "MIXTURE", "PERENNIALS",ifelse(site_data$VEGCLASS == "PERENNIALS", "PERENNIALS", "FLAG")))) 
    # We now need to calculate a new annual correction factor for BNPP
    BNPP.Correction = ifelse(VEGCLASS.project == "ANNUALS", APC.Annual, ifelse(VEGCLASS.project == "BARE",APC.Annual, ifelse(VEGCLASS.project == "MIXTURE",  APC.Mixture, APC.Perennial )))
    
    # Annual aboveground production under grazing
    ANPPest = ProdG * BNPP.Correction
    
    # Other Productivity Calculations
    # Annual belowground productivity under grazing
    # Increase in biomass from initial amount (S0) to a biomass Sk after G days
    Pu = Sk - S0
    # Equation 16 - Annual below ground productivity under grazing
    BNPPest = ((BNPP.MAP.B0 * MAP.strata) - (BNPP.MAT.B1 * MAP.strata^2) + (BNPP.MAT.B0 * MAT.strata)) * (ANPPest/ANPPmax) * BNPP.Correction
    
    # Soil Organic Carbon (SOC) calculations.
    # Calculate Mean proportion of plant as lignin and cellulose using dominant vegetation class.
    LIGCELL = ifelse(VEGCLASS.project == "ANNUALS", LigCell.Annual, ifelse(VEGCLASS.project == "BARE",LigCell.Annual, ifelse(VEGCLASS.project == "MIXTURE", LigCell.Mixture, LigCell.Perennial)))
    
    # Plant-derived soil carbon input 
    PDSOC.project = (0.45 * (((LIGCELL/100) * ANPPest * (1 - site_data$FIRE.project)) + (((LIGCELL/100) + 0.05) * BNPPest))) * (-0.3559 + 0.3914 * log(site_data$DEPTH))
    # Dung-derived soil carbon input 
    DDSOC.project = ((LIGCELL/100) * (CONSUMED + CONSUMED.OFF) * 0.45 * (0.4718 - 0.0009*(D))) * (-0.3559 + 0.3914 * log(site_data$DEPTH))
    
    # WETDAYS (eq 19)
    WETDAYS = (WETDAYS.B1 * MAP.strata - WETDAYS.Bo) * G
    
    # SOCeq
    SOCeq.large.proj  = (((PDSOC.project+DDSOC.project) / (WETDAYS * (0.7 + 0.3 *(site_data$SAND/100)) *  SOC.B3)) + (SOC.B2/ SOC.B3)) / 100 
    SOCeq.small.proj  = ((PDSOC.project+DDSOC.project) / (WETDAYS * (0.7+0.3*(site_data$SAND/100)) * exp(SOC.B0)))^(1/SOC.B1) / 100  # technically the ms states SOC <4600 g/m2
    # this looks at the two SOC calculations and selects the lower one, which is the one using the correct piece wise curve
    SOCeq.project       = pmin(SOCeq.large.proj, SOCeq.small.proj, na.rm = TRUE)

    # Lets calculate the delta - new SOCeq.project - SOC.predicted by modeal
    SOC.diff <- SOCeq.project - site_data$SOC.measured
    
    # Store results in a dataframe
    # First, it will output parameters to confirm monte carlo analysis is randomly drawing. 
    site_results[i, "r"]                   <- r
    site_results[i, "G.MAT.B0"]            <- G.MAT.B0
    site_results[i, "G.MAT.B1"]            <- G.MAT.B1
    site_results[i, "G.MAP.B0"]            <- G.MAP.B0    
    site_results[i, "Cg.B0"]               <- Cg.B0         
    site_results[i, "Cg.B1"]               <- Cg.B1         
    site_results[i, "LigCell.Annual"]      <- LigCell.Annual   
    site_results[i, "LigCell.Mixture"]     <- LigCell.Mixture  
    site_results[i, "LigCell.Perennial"]   <- LigCell.Perennial
    site_results[i, "ANPPmax.Sand.B0"]     <- ANPPmax.Sand.B0 
    site_results[i, "ANPPmax.Sand.B1"]     <- ANPPmax.Sand.B1  
    site_results[i, "ANPPmax.MAP.B0"]      <- ANPPmax.MAP.B0
    site_results[i, "ANPPmax.MAT.B0"]      <- ANPPmax.MAT.B0  
    site_results[i, "ANPPmax.MAT.B1"]      <- ANPPmax.MAT.B1   
    site_results[i, "APC.Annual"]          <- APC.Annual   
    site_results[i, "APC.Mixture"]         <- APC.Mixture   
    site_results[i, "APC.Perennial"]       <- APC.Perennial
    site_results[i, "BNPP.MAP.B0"]         <- BNPP.MAP.B0   
    site_results[i, "BNPP.MAT.B0"]         <- BNPP.MAT.B0      
    site_results[i, "BNPP.MAT.B1"]         <- BNPP.MAT.B1     
    site_results[i, "WETDAYS.Bo"]          <- WETDAYS.Bo     
    site_results[i, "WETDAYS.B1"]          <- WETDAYS.B1 
    site_results[i, "SOC.B0"]              <- SOC.B0         
    site_results[i, "SOC.B1"]              <- SOC.B1         
    site_results[i, "SOC.B2"]              <- SOC.B2      
    site_results[i, "SOC.B3"]              <- SOC.B3 
    site_results[i, "SOC.B2"]              <- SOC.B2      
    site_results[i, "MAP.strata"]          <- MAP.strata 
    site_results[i, "MAT.strata"]          <- MAT.strata 
    
    # Note that this is not every intermediate output, but a selection of important steps to facilitate cross-checking 
    site_results[i, "ANPPmax"]             <- ANPPmax
    site_results[i, "VEGCLASS.project"]    <- VEGCLASS.project
    site_results[i, "WETDAYS"]             <- WETDAYS
    site_results[i, "ANPPest"]             <- ANPPest
    site_results[i, "BNPPest"]             <- BNPPest
    site_results[i, "PDSOC.project"]       <- PDSOC.project
    site_results[i, "PDSOC.project"]       <- DDSOC.project
    site_results[i, "SOCeq.large.proj"]    <- SOCeq.large.proj
    site_results[i, "SOCeq.small.proj"]    <- SOCeq.small.proj
    site_results[i, "SOCeq.project"]       <- SOCeq.project
    site_results[i, "SOC.diff"]            <- SOC.diff
  }
  
  # Append site-specific results to the list
  results_list[[site]] <- site_results
}

# Combine results into a single data frame
MC.Iterations.all <- bind_rows(results_list)


# The number of iterations and sites will quickly get out of hand for excel. 
# So here, we print the first 100 rows for each Site for troubleshooting. 
# Simple change the n = in the slice_head() function below to print desired number of iterations per site.
MC.Iterations.preview <- MC.Iterations.all %>%
  group_by(Site) %>%
  slice_head(n = 100) %>%  
  ungroup() 

# Now save to workbook
addWorksheet(SNAPGRAZE.wb, "MonteCarlo_100Iterations")
writeData(SNAPGRAZE.wb, "MonteCarlo_100Iterations", MC.Iterations.preview)
saveWorkbook(SNAPGRAZE.wb, paste0(out_dir, "SNAPGRAZEv23_", Project_initials, "_", date, ".xlsx"), overwrite = TRUE)

#### * Summary Table ####
# The MC.Iterations.all table is long, a row for each iteration and each station, 
# this aggregates summary stats for each site.
MC.Site <- MC.Iterations.all %>%
  group_by(Site, Community) %>%
  summarise(
    mean.SOC.diff = mean(SOC.diff, na.rm=TRUE),
    se.SOC.diff = se(SOC.diff, na.rm=TRUE),
    Uncertainty.perc = (qnorm(0.975)^2)*se.SOC.diff/abs(mean.SOC.diff)
  )%>%
  # This is specific for a project to reorder sites so that XX 01 is followed by XX 03 rather than XX 011, etc. 
  # Extract prefix and number (e.g., "XXX", 4)
  mutate(
    Site_prefix = str_extract(Site, "^[A-Z]+"),
    Site_number = as.numeric(str_extract(Site, "\\d+$"))
  ) %>%
  arrange(Site_prefix, Site_number) %>%
  select(-Site_prefix, -Site_number)
MC.Site

# And the summary by site can be added to the excel workbook
addWorksheet(SNAPGRAZE.wb, "MonteCarlo_Site")
writeData(SNAPGRAZE.wb, "MonteCarlo_Site", MC.Site)
saveWorkbook(SNAPGRAZE.wb, paste0(out_dir, "SNAPGRAZEv23_", Project_initials, "_", date, ".xlsx"), overwrite = TRUE)

2+2
