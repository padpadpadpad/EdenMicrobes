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

- **Sequencing Methods**
At each sample plot, from each corner of a 2 x 2 m quadrat, four soil samples were collected using sterile soil augers. Any leaf litter was removed from the surface, before approximately 200g of soil was collected, typically representing 3 auger cores worth of material, from within the first 10 cm of the topsoil. Auger blades were “cleaned” by immersion in soil adjacent to the collection point prior to each sampling event. Blades were changed between sampling of the humid and dry biomes. The soil cores were sealed in sterile plastic bags and transported to the onsite laboratory for immediate processing.

Two eDNA extractions were performed from a 15 g subsample of each of the 128 samples following methods developed by Taberlet et al. (2012b) and modified by Zinger et al. (2016). DNA extractions were conducted in the field lab less than 2 hours after collection using a NucleoSpin® Soil kit (Machery Nagel, Duren, Germany). Eight negative controls were included for a total of 256 DNA extractions. The last elution step of the DNA extraction protocol was not carried out on site, with DNA on the column instead stored with silica gel and transported back to the EDB lab in Toulouse for subsequent steps. 

PCRs were performed in triplicate, meaning that each sample was extracted twice, and each extract was amplified 3 times, resulting in 6 replicates for each sample in total. Each PCR reaction was performed in a total volume of 20 μl and comprised 10 μl of AmpliTaq Gold Master Mix (Life Technologies, Carlsbad, CA, USA), 5.84 μl of Nuclease-Free Ambion Water (Thermo Fisher Scientific, Massachusetts, USA), 0.25 μM of each primer, 3.2 μg of BSA (Roche Diagnostic, Basel, Switzerland), and 2 μl DNA template that was 10-fold diluted to reduce PCR inhibition. PCR was conducted by targeting a range of barcode regions to sequence for Eukaryotes, Fungi and Bacteria under the following conditions : 
- the V7 region of the 18S rRNA gene as a diagnostic marker for Eukaryotes with the following universal primers : forward (5’-TCACAGACCTGTTATTGC-3’), and reverse (5’-TTTGTCTGCTTAATTSCG-3’) (Guardiola et al. 2015). PCR was conducted with 35 cycles, denaturation at 95°C for 30s, annealing at 45°C for 30s, elongation at 72°C for 60s, with the final elongation for 7 mins. 
- the V5-V6 regions of the 16S rRNA gene as a diagnostic marker for Bacteria with the following universal primers : forward (5’-GGATTAGATACCCTGGTAGT-3’), and reverse (5’-CACGACACGAGCTGACG-3’) (Fliegerova et al. 2014). PCR was conducted with 30 cycles, denaturation at 95°C for 30s, annealing at 57°C for 30s, elongation at 72°C for 90s, with the final elongation for 7 mins. 
- the ITS1 region of the nuclear ribosomal RNA genes as a diagnostic marker for Fungi with the following universal primers : forward (5’-CAAGAGATCCGTTGTTGAAAGTK-3’), and reverse (5’-GGAAGTAAAAGTCGTAACAAGG-3’) (Epp et al 2012; Taberlet et al 2018). PCR was conducted with 35 cycles, denaturation at 95°C for 30s, annealing at 55°C for 30s, elongation at 72°C for 60s, with the final elongation for 7 mins.

Three negative PCR controls per plate were amplified and sequenced in parallel with the regular samples. Three positive controls were also included and consisted of XX TO CONFIRM XX.  Six wells per PCR plate were left empty (non-used tag combinations) to control for tag jumps which can occur during amplification and sequencing (see below for downstream data curation). All PCR products were pooled and the library was constructed using the Illumina TruSeq NanoPCRFree kit following the supplier’s instructions (Illumina Inc., San Diego, California, USA). Sequencing was performed on a Hiseq run (Illumina platform,San Diego, CA, USA) at the GeT-Plage platform (Toulouse, France).

- **Bioinformatics Methods**
All bioinformatic analyses were performed using the GenoToul bioinformatics platform, with sequence reads processed using the OBITOOLS package (Boyer et al. 2016) and R scripts (R Core Team, 2020) following the procedure described in Zinger et al. (2019). To consolidate ASVs, tools from Phyloseq to merge files were deployed...

Bioinformatic analyses were performed on the GenoToul bioinformatics platform (Toulouse, France), with the OBITOOLS package (Boyer et al. 2016). First, ‘illuminapairedend’ was used to assemble paired-end reads. This algorithm is based on an exact alignment algorithm that considers the quality scores at all positions during the assembly process. Subsequently, we used the ‘ngsfilter’ command to identify and remove the primers and tags on each read, and assign reads to their respective samples. This program was used with its default parameters tolerating two mismatches for each of the two primers and no mismatch for the tags. Following this, sequencing reads were dereplicated using the ‘obiuniq’ command. Sequences of low quality (containing Ns or with paired-end alignment scores below 50) were excluded using the ‘obigrep’ command. The same command was used to exclude sequences represented by only one read (singletons) as they are more likely to be molecular artefacts (Taberlet et al. 2018). Sequences outside of the preset range were also discarded (90-200 in length for Eukaryotes; 30-400 for Bacteria; 30-900 for Fungi). 

To assign a taxon to fungal ASVs, we built a reference sequence database using the ecoPCR programme (Ficetola et al. 2010) and the fungi specific markers on the European Molecular Biology Laboratory (EMBL; release 141). ASVs were then assigned a taxonomy, using OBITOOL’s ecotag programme (Boyer et al. 2016), which performs a global alignment of each ASV sequence (the query) against each reference. The reference taxon assigned to each ASV corresponds to the Last Common Ancestor of all the best-match sequences for the query. For taxonomic assignment of bacteria and eukaryote OTUs, the SILVA taxonomic database was used (version 1.3; Quast et al., 2012). Classification was performed by a local nucleotide BLAST search against the non-redundant version of the SILVA SSU Ref dataset (release 132; http://www.arb-silva.de) using blastn (version 2.2.30+; http://blast.ncbi.nlm.nih.gov/Blast.cgi) with standard settings (Camacho et al., 2009). 

Datasets were subsequently filtered to remove contaminants as well as artefacts such as PCR chimeras and remaining sequencing errors, following Zinger et al. (2019) and using the metabaR R package (Zinger et al 2020b), in R version 3.6.1 (R Development Core Team, 2013). The filtering process consisted of four steps: (i) a negative control-based filtering. ASVs whose maximum abundance was found in extraction/PCR negative controls were removed from the dataset, as they were likely to be reagent/aerosol contaminants, better amplified in the absence of competing DNA fragments as it is the case in biological samples. (ii) a reference-based filtering. ASVs which are too dissimilar from sequences available in reference databases are potential chimeras generated during sequencing and amplification. In this study, we chose to set similarity thresholds at 80% for bacteria and eukaryotes and due to the marker being more polymorphic, 65% for fungi. In addition, we removed all taxa that are not targeted by the primer used. (iii) an abundance-based filtering. This procedure targets incorrect assignment of a few numbers of sequences corresponding to true ASVs occurring to the wrong sample, a phenomenon called “tag-switching” (Esling et al. 2015), “tag jumps” (Schnell et al. 2015) or “cross-talk” (Edgar 2018). It consists in setting ASVs abundances to 0 in samples where their abundance represents < 0.03% of the total OTU abundance in the entire dataset. (iv) Finally, we conducted a PCR-based filtering by considering any PCR reaction that yielded less than 1000 reads for fungi, bacteria and eukaryotes as non-functional, and removed them from the dataset. The number of reads, ASVs and PCRs removed at each stage for each marker are detailed in the table below.

| uniq    |          | no_singletons |          | PCRs above threshold (>1000) | Extraction Contaminent Removal | PCR Contaminent Removal | Sequencing Contaminent Removal |
| --------- | -------- | -------- | ------------ | -------- | -------- | -------- | ------- | -------- | ------------- | -------- | ---------------------------- | ------------------------------ | ----------------------- | ------------------------------ |
| Stage     | ASV      | Reads    | ASV          | Reads    | ASV      | Reads    | ASV     | Reads    | ASV           | Reads    | ASV                          | ASV                            | ASV                     | ASV                            |
| ITS1-Rep1 | 3205177  | 3205177  | 2255445      | 2255445  | 2243828  | 2243828  | 469368  | 2243828  | 76499         | 1850959  | 176                          | 14                             | 51                      | 12                             |
| ITS1-Rep2 | 2754254  | 2754254  | 1932598      | 1932598  | 1919187  | 1919187  | 416437  | 1919187  | 63812         | 1566562  | 128                          | 21                             | 3                       | 20                             |
| ITS1-Rep3 | 4054254  | 4054254  | 2800343      | 2800343  | 2786030  | 2786030  | 582040  | 2786030  | 94062         | 2298052  | 244                          | 15                             | 88                      | 21                             |
| 16S-Rep1  | 4539219  | 4539219  | 2697244      | 2697244  | 2682047  | 2682047  | 1504427 | 2682047  | 139331        | 1316951  | 237                          | 287                            | 71                      | 244                            |
| 16S-Rep2  | 4109968  | 4109968  | 2421281      | 2421281  | 2407686  | 2407686  | 1370039 | 2407686  | 125006        | 1162653  | 231                          | 260                            | 38                      | 229                            |
| 16S-Rep3  | 8175582  | 8175582  | 4813103      | 4813103  | 4786542  | 4786542  | 2575508 | 4786542  | 246502        | 2457536  | 244                          | 618                            | 109                     | 470                            |
| 18S-Rep1  | 16837584 | 16837584 | 14364998     | 14364998 | 14281333 | 14281333 | 677100  | 14281333 | 152495        | 13756728 | 250                          | 67                             | 37                      | 270                            |
| 18S-Rep2  | 12133280 | 12133280 | 10412433     | 10412433 | 10345662 | 10345662 | 503842  | 10345662 | 120786        | 9962606  | 249                          | 61                             | 20                      | 172                            |
| 18S-Rep3  | 9110898  | 9110898  | 7732534      | 7732534  | 7694307  | 7694307  | 408830  | 7694307  | 95583         | 7381060  | 244                          | 65                             | 19                      | 150                            |


Upon filtering completion, remaining PCRs per technical replicate were summed and the read count of technical replicates was normalised to reduce potential bias caused by PCR stochasticity and differential sequencing efforts. Standardization consisted in randomly resampling (with replacement) a number of reads that corresponded to the first quartile of the total read number for reads per samples. This returns samples with a read count equal across all samples, whilst maintaining sample specific OTU relative abundances. To do so, each OTU in each sample was resampled with replacement a thousand times, following an approach detailed by Veresoglou et al. (2019). Finally, in order to reduce stochastic variation of taxa from one soil core to another, and to match DNA sequencing data with the soil chemistry ones, the four replicate samples within each subplot were aggregated by summing reads (after normalisation).


### Climate data (in data/climate)

- **Eden_microclimare** To characterise variation in soil temperature and moisture across both biomes, Aranet Substrate Sensors (QH21142) were deployed at each of the sample locations twice for a period of 1 week between 28th of April – 30th of June for the Rainforest Biome and October the 6th – 15th December for the Mediterranean biome. The three metal probe spikes of each logger were buried, with their tips at approx 5cm depth.These took a measure of temperature (degree C) and moisture (volumetric water content, VMC) every 5 minutes, before connecting to the Eden 5G cloud and storing data. Since only 14 probes were available, deployment was rotated across plots, with two probes deployed in one of the 32 plots at any one time.  Upon data processing, the dataset was reduced in dimension by calculating an hourly average from the collected dataset.

- **Eden_soil_chem_data_nov23.csv** contains nutrient and chemical composition data for each sample. Bulk soil samples were processed by NRM laboratory (Berkshire, UK) in 2019 to determine pH, available Phosphorus, available Potassium, available Magnesium, total Nitrogen, total Carbon, total Phosphorus, Manganese, Iron, and the relative percentage of Sand, Silt and Clay.  Further soil analysis conducted in 2021, was used to characterise sample site soil respiration, soil moisture, soil organic matter content, and total root biomass.
For further information see Duley et al. (2023) and supplementary materials here (https://onlinelibrary.wiley.com/doi/abs/10.1111/rec.13831). A breakdown of all the column names and their units is below:

  - Sample - Sample name that corresponds to a sequencing file and a
    sampling point.

  - Ecosystem - What ecosystem did the sample come from.

  - Diversity - Whether the sample has high or low plant diversity.

  - Biome - Which biome did the sample come from (Humid = rainforest,
    Dry = mediterranean).

  - Soil_OM / GWC - Performed on samples aggregated at the quadrat level. Loss on ignition was used as a proxy for soil moisture and organic matter content (Heiri et al. 2001), following methods developed by Jensen et al. (2018) and Joy et al. (2021). Approximately 40g of soil from each pooled sample was placed in a foil tray and dried in an oven at 105°C for 24 hours, these were reweighed to determine gravimetric moisture content (GWC). From this, subsamples of 5g were placed in crucibles and transferred to a muffle furnace at 600°C for 4 hours. Samples were cooled by being left in the furnace overnight after being turned off, before being reweighed.

  - Root_biomass - Performed on samples aggregated at the quadrat level. Roots were extracted from the remainder of pooled soil samples using methods developed by Frasier et al. (2016). To separate the roots from the soil the samples were washed through a submerged 250 μm sieve with running tap water, and larger soil aggregates were broken down by hand. Roots were then collected from the surface of the sieve using tweezers, placed in foil trays, weighed and oven dried at 60°C for 12 hours, then reweighed to give a representative root biomass for each plot.

  - pH - Performed on samples aggregated at the quadrat level. Measured in water (15 g fresh weight soil suspended in 20ml of deinonised water) using a Jensen desktop probe, calibrated with standard buffer solutions of pH 7 and pH 4.


  - Soil_respiration - Performed on samples aggregated at the quadrat level. Measured using a TARGAS-1 CO2/H2O infrared gas analyser with soil respiration chamber (PP systems 2016).

Here and below, methods conducted by NRM Cawood Laboratories following their set methodology, with samples aggregated at the habitat level.

  - pH2 -  Measured in water (15 g fresh weight soil suspended in 20ml of deinonised water) using a Jensen desktop probe, calibrated with standard buffer solutions of pH 7 and pH 4.

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
