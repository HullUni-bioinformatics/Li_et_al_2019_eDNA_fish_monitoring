# Li_et_al_2018_eDNA_fish_monitoring

Data processing workflow and supplementary data for Li et al. 2018 - Development of an environmental DNA method for monitoring fish communities: ground truthing in diverse lakes with characterised fish faunas

Release 1.0 of this repository has been archived: 

##Contents
 
  - SRA accession numbers for 12S raw Illumina data ([here](https://github.com/HullUni-bioinformatics/Li_et_al_2018_eDNA_fish_monitoring/blob/master/12S/supplementary_data/Sample_accessions.tsv))
  
  - SRA accession numbers for Cytb raw Illumina data ([here](https://github.com/HullUni-bioinformatics/Li_et_al_2018_eDNA_fish_monitoring/blob/master/Cytb/supplementary_data/Sample_accessions.tsv))
  
  - Jupyter notebook for download the 12S raw read data ([here](https://github.com/HullUni-bioinformatics/Li_et_al_2018_eDNA_fish_monitoring/blob/master/12S/How_to_download_12S_Rawdata_from_SRA.ipynb)).
  
  - Jupyter notebook for download the Cytb raw read data ([here](https://github.com/HullUni-bioinformatics/Li_et_al_2018_eDNA_fish_monitoring/blob/master/Cytb/How_to_download_Cytb_Rawdata_from_SRA.ipynb)).
  
   - 12S curated reference databases used in analyses in Genbank format ([here](https://github.com/HullUni-bioinformatics/Li_et_al_2018_eDNA_fish_monitoring/tree/master/12S/supplementary_data/reference_DBs))
   
  - Cytb curated reference databases used in analyses in Genbank format ([here](https://github.com/HullUni-bioinformatics/Li_et_al_2018_eDNA_fish_monitoring/tree/master/Cytb/supplementary_data/reference_DBs))
  
  - Jupyter notebook for fully rerun/reproduce 12S analyses ([here](https://github.com/HullUni-bioinformatics/Li_et_al_2018_eDNA_fish_monitoring/blob/master/12S/LP_12S.ipynb)).
  
  - Jupyter notebook for fully rerun/reproduce Cytb analyses ([here](https://github.com/HullUni-bioinformatics/Li_et_al_2018_eDNA_fish_monitoring/blob/master/Cytb/LP_Cytb.ipynb)).
    
  - R scripts used to produce the figures in the paper ([here]())
  
  - Appendix S3: Read counts of Cytb OTUs data was used for the R script (.csv) ([here]())
  
  - Appendix S4: Read counts of 12S OTUs data was used for the R script (.csv) ([here]())

##Instructions on how to set up all dependencies for data processing/analyses
 
To facilitate full reproducibility of our analyses we provide Jupyter notebooks illustrating our workflow and all necessary supplementary data in this repository.

Illumina data was processed (from raw reads to taxonomic assignments) using the [metaBEAT](https://github.com/HullUni-bioinformatics/metaBEAT) pipeline. The pipeline relies on a range of open bioinformatics tools, which we have wrapped up in a self contained docker image which includes all necessary dependencies [here](https://hub.docker.com/r/chrishah/metabeat/).

##Setting up the environment

In order to retrieve supplementary data (reference sequences etc.) start by cloning this repository to your current directory:
```
git clone --recursive https://github.com/HullUni-bioinformatics/Li_et_al_2018_eDNA_fish_monitoring.git
```



In order to make use of our self contained analysis environment you will have to install [Docker](https://www.docker.com/) on your computer. Docker is compatible with all major operating systems. See the [Docker documenation](https://docs.docker.com/) for details. On Ubuntu installing Docker should be as easy as:

```
sudo apt-get install docker.io
```

Once Docker is installed you can enter the environment by typing, e.g.:
```
sudo docker run -i -t --net=host --name metaBEAT -v $(pwd):/home/working chrishah/metabeat /bin/bash
```

This will download the metaBEAT image (if it's not yet present on your computer) and enter the 'container', i.e. the self contained environment (Note that `sudo` may be necessary in some cases). With the above command the container's directory `/home/working` will be mounted to your current working directory (as instructed by `$(pwd)`), in other words, anything you do in the container's `/home/working` directory will be synced with your current working directory on your local machine. 

##Data processing workflow as Jupyter notebooks

  - __12S__
 
Raw illumina data has been deposited with Genbank (BioProject: PRJNA454866; BioSample accession: SAMN09058600-SAMN09058613; Sequence Read Archive accessions: SRR7106552-SRR7106565) - see sample specific accessions [here](https://github.com/HullUni-bioinformatics/Li_et_al_2018_eDNA_fish_monitoring/blob/master/12S/supplementary_data/Sample_accessions.tsv). Before following the workflow below, you'll need to download the raw reads from SRA. To __download the 12S raw read data__ you can follow the steps in [this Jupyter notebook](https://github.com/HullUni-bioinformatics/Li_et_al_2018_eDNA_fish_monitoring/blob/master/12S/How_to_download_12S_Rawdata_from_SRA.ipynb).


With the data in place you should be able to __fully rerun/reproduce our analyses__ by following the steps outlined in the [this Jupyter notebook](https://github.com/HullUni-bioinformatics/Li_et_al_2018_eDNA_fish_monitoring/blob/master/12S/LP_12S.ipynb).

The workflow illustrated in the notebooks assumes that the raw Illumina data is present in a directory `raw_data` at the base of the repository structure and that the files are named according to the following convention:
'sampleID-marker', followed by '_R1' or '_R2' to identify the forward/reverse read file respectively. sampleID must corresponds to the first column in the file `Sample_accessions.tsv` [here](https://github.com/HullUni-bioinformatics/Li_et_al_2018_eDNA_fish_monitoring/blob/master/12S/supplementary_data/Sample_accessions.tsv).

The __Querymap for demultiplex and trimming__ can be found [here](https://github.com/HullUni-bioinformatics/Li_et_al_2018_eDNA_fish_monitoring/blob/master/12S/Querymap_demultiplex_trimming_12S.txt)

The __12S reference sequences__ (curated reference databases) used in analyses in Genbank format ([here](https://github.com/HullUni-bioinformatics/Li_et_al_2018_eDNA_fish_monitoring/tree/master/12S/supplementary_data/reference_DBs))
 

  - __Cytb__

Raw illumina data has been deposited with Genbank (BioProject: PRJNA454866; BioSample accession: SAMN09059971-SAMN09060202; Sequence Read Archive accessions: SRR7108634-SRR7108885) - see sample specific accessions [here](https://github.com/HullUni-bioinformatics/Li_et_al_2018_eDNA_fish_monitoring/blob/master/Cytb/supplementary_data/Sample_accessions.tsv). Before following the workflow below, you'll need to download the raw reads from SRA. To __download the raw read data__ you can follow the steps in [this Jupyter notebook](https://github.com/HullUni-bioinformatics/Li_et_al_2018_eDNA_fish_monitoring/blob/master/Cytb/How_to_download_Cytb_Rawdata_from_SRA.ipynb).


With the data in place you should be able to __fully rerun/reproduce our analyses__ by following the steps outlined in the [this Jupyter notebook](https://github.com/HullUni-bioinformatics/Li_et_al_2018_eDNA_fish_monitoring/blob/master/Cytb/LP_Cytb.ipynb).

The workflow illustrated in the notebooks assumes that the raw Illumina data is present in a directory `raw_reads` at the base of the repository structure and that the files are named according to the following convention:
'sampleID-marker', followed by '_R1' or '_R2' to identify the forward/reverse read file respectively. sampleID must corresponds to the first column in the file `Sample_accessions.tsv` [here](https://github.com/HullUni-bioinformatics/Li_et_al_2018_eDNA_fish_monitoring/blob/master/Cytb/supplementary_data/Sample_accessions.tsv).


The __Cytb reference sequences__ (curated reference databases) used in analyses in Genbank format ([here](https://github.com/HullUni-bioinformatics/Li_et_al_2018_eDNA_fish_monitoring/tree/master/Cytb/supplementary_data/reference_DBs))

