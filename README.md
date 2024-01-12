
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EdenMicrobes

<!-- badges: start -->
<!-- badges: end -->

The goal of **EdenMicrobes** is to document and track the progress of
the datasets and analyses related to microbes at the Eden Project.

## Outline

Outline the project here!

## Scripts

- **visualise_climate.R** looks at the data from the temperature and
  moisture probes that were placed in the habitats. Produces
  **temperature_moisture.png**.

## Datasets

We have sequencing datasets of the biomes contained in `data/sequencing`
and data from monitoring of abiotic variables present in `climate_data`.

- **Eden_GPS_data.csv** contains the GPS coordinates of each site. A
  breakdown of the columns are below:

  - Sample - Sample name that corresponds to a sequencing file and a
    sampling point.

  - Ecosystem - What ecosystem did the sample come from.

  - Diversity - Whether the sample has high or low plant diversity.

  - Biome - Which biome did the sample come from (Humid = rainforest,
    Dry = mediterranean).

  - GPS -

  - X -

  - Y -

### Sequencing data (in data/sequencing)

### Climate data (in data/climate)

- **Eden_microclimare**

<!-- -->

- **Eden_soil_chem_data_nov23.csv** contains nutrient and chemical composition data for each sample. A breakdown of all the column names and their units is below, for further information see Duley et al. (2023) and supplementary materials here (https://onlinelibrary.wiley.com/doi/abs/10.1111/rec.13831):


  - Sample - Sample name that corresponds to a sequencing file and a
    sampling point.

  - Ecosystem - What ecosystem did the sample come from.

  - Diversity - Whether the sample has high or low plant diversity.

  - Biome - Which biome did the sample come from (Humid = rainforest,
    Dry = mediterranean).

  - Soil_OM / GWC -Loss on ignition was used as a proxy for soil moisture and organic matter content (Heiri et al. 2001), following methods developed by Jensen et al. (2018) and Joy et al. (2021). Approximately 40g of soil from each pooled sample was placed in a foil tray and dried in an oven at 105°C for 24 hours, these were reweighed to determine gravimetric moisture content (GWC). From this, subsamples of 5g were placed in crucibles and transferred to a muffle furnace at 600°C for 4 hours. Samples were cooled by being left in the furnace overnight after being turned off, before being reweighed.

  - Root_biomass -Roots were extracted from the remainder of pooled soil samples using methods developed by Frasier et al. (2016). To separate the roots from the soil the samples were washed through a submerged 250 μm sieve with running tap water, and larger soil aggregates were broken down by hand. Roots were then collected from the surface of the sieve using tweezers, placed in foil trays, weighed and oven dried at 60°C for 12 hours, then reweighed to give a representative root biomass for each plot.

  - pH - Measured in water (15 g fresh weight soil suspended in 20ml of deinonised water) using a Jensen desktop probe, calibrated with standard buffer solutions of pH 7 and pH 4.


  - Soil_respiration - measured using a TARGAS-1 CO2/H2O infrared gas analyser with soil respiration chamber (PP systems 2016).

  - pH2 - Here and below, methods conducted by NRM Cawood Laboratories following their set methodology. Measured in water (15 g fresh weight soil suspended in 20ml of deinonised water) using a Jensen desktop probe, calibrated with standard buffer solutions of pH 7 and pH 4.

  - Phosphorus --	Phosphorus - Sodium Bicarbonate Extractable - Olsens - reported as mg/l dry basis.

  - Potassium - Ammonium Nitrate Extractable - reported as mg/l dry basis.

  - Magnesium - Ammonium Nitrate Extractable - reported as mg/l dry basis.

  - Sand - Textural Classification Reported as % w/w dry matter basis.

  - Silt - Textural Classification Reported as % w/w dry matter basis.

  - Clay - Textural Classification Reported as % w/w dry matter basis.

  - Nitrogen - Dumas - reported as % w/w dry basis.

  - Manganese - DTPA Extractable Reported as mg/l dry basis.

  - Iron - DTPA Extractable Reported as mg/l dry basis.

  - Phosphorus2 - Reported as mg/kg or % w/wdry basis.

  - Carbon - Dumas - reported as % w/w dry basis.

  - C:N - Dumas - reported as % w/w dry basis.

## Contact

This project is primarily a collaboration between Daniel Padfield
(<d.padfield@exeter.ac.uk>) at the University of Exeter and Julian
Donald (formerly of the Eden Project).
