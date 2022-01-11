[![DOI](https://zenodo.org/badge/190763363.svg)](https://zenodo.org/badge/latestdoi/190763363)

# Supporting data and code for: Report of a new resistance to carbendazim in *Fusicoccum amygdali*, the causal agent of constriction canker of peach and almond trees

![alt text](https://y9xjpa.db.files.1drv.com/y4mpjcgFHlzSPQzhQNV-dAkypEx7uWgJhPxL6MraJDeaIrM-2RuAwYTpqu6gxKuqMRb9SuyCU4qQntQUYMi3VF9MkovXWgC_P5rG5ny6lLdKjXJBjVdPkGUuy3TLdyJC6jCzGDkk2SOaGWwlQCHLn-uyF7zSROc9SEcCcHqbh9OSCh2vv6eWVhUQ-vt5FdOrxCRxSXRJ_c1xKVKZtnYyOViXw?width=1584&height=588&cropmode=none)


## Context
The constriction canker of peach and almond trees, caused by *Fusicoccum amygdali*, is an important disease which is mainly controled using fungicides. In this study, we report a first case of resistance to carbendazim for this pathogen in France. We first screened a set of populations using bioassay. After identifying a population displaying evidence of resistance, isolates were cultivated on medium and tested with a bioassay, confirming the presence of highly resistant phenotypes. A fragment of the beta-tubulin gene were sequenced from isolates displaying resistance or sensitive phenotypes. A mutation leading to the E198K substitution was identified in resistant isolates. The population from which resistant isolates have been identified is located in Corsica Island.


## Datasets
There are five data sets used in this study. The files can be found in the "data" folder. 

+ **COM_SHP.RDATA:** the first data set contains geographical information to plot the french administrative layer 'commune'

+ **DEP_SHP.RDATA:** the first data set contains geographical information to plot the french administrative layer 'departement'

+ **REG_SHP.RDATA:** the second data set contains geographical information to plot the french administrative layer 'regions'

These three geographical data files were obtained using the data from the [IGN website](http://professionnels.ign.fr/adminexpress). The version of the data used is the "Edition Novembre 2017". The shapefile were simplified using [mapshaper](https://mapshaper.org/) online tool in order to reduce the size of the files. 

+ **fusicodat.txt:** The third dataset contains the results of the dose-response bioassays conducted on populations or isolates. 
  + *strain_ID*: the ID of the population or the strain
  + *species*: species name
  + *strain_type*: type of biological sample. Possible values are 'population', 'individual' or 'reference' (for the reference strain)
  + *active_substance*: the tested active substance
  + *sampling_year*: the year the sample was sampled in the field
  + *nb_day*: number of day post inoculation on plate with medium
  + *dose*: the concentration of the active substance in the medium (mg/L)
  + *perc_croiss*: the ratio between the measured mycelial growth for a given dose and the measured mycelial growth at dose '0 mg/L' (*ie* without active substance) for individuals OR the ratio between the average elongation of the germ tube for a given dose and the the average elongation of the germ tube at dose '0 mg/L' (*ie* without active substance)

+ **fusipopdat.txt:** the fourth dataset contains the information to plot the studied populations on the map. 
  + *pop_ID*: the ID of the sampled population
  + *carbend_R*: the resistance status of the population, at least confirmed by founding a strain resistant to a mycelial growth bioassay
  + *sampling_year*: the year the population was sampled
  + *departement*: the french 'departement' which the population was sampled from


## R scripts
+ **load_fusico_data.R:** the script to load the different datasets in the R environment. 
+ **fusico_bioassay_pop.R:** the script to perform the regression analysis on the populations germination inhibition bioassay results, with carbendazim. It also contains the code for producing Figure 1. 
+ **fusico_bioassay_ind.R:** the script to perform the regression analysis on the individuals growth inhibition bioassay results, with both carbendazim and diethofencarb. It also contains the code for producing Figure 3. 
+ **fusico_map.R:** the script for producing the map (Figure 2). 


## Citation
You can cite the related study as follow: 
+ Fontaine S, Caddoux L, Remuson F, Barrès B. [Report of a new resistance to carbendazim in *Fusicoccum amygdali*, the causal agent of constriction canker of peach and almond trees. *Accepted for publication in Plant Pathology*.](https://bsppjournals.onlinelibrary.wiley.com/doi/10.1111/ppa.13525)

If you want to use (some of) the code found on this page or if you want to cite this repository: 
+ Benoit Barrès. [Supporting data and code for: Report of a new resistance to carbendazim in *Fusicoccum amygdali*, the causal agent of constriction canker of peach and almond trees. Zenodo; 2022.](https://zenodo.org/badge/latestdoi/190763363)
