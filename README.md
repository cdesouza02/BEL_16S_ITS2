# Project Goals
To understand how 6 coral host species interact with their holobionts under stress (disease and heatwave).

## Questions
### BEL_16S_outputs
1. [Explore read counts after dada2 filtering](https://github.com/cdesouza02/BEL_16S_ITS2/blob/main/BEL_16S_outputs/1_distribution_of_reads_post_dada2.ipynb)
2. [Troubleshooting how to combine 5 dada2 runs](https://github.com/cdesouza02/BEL_16S_ITS2/blob/main/BEL_16S_outputs/2_combine_dada2_output.ipynb)
3. [Creating 16S phyloseq object](https://github.com/cdesouza02/BEL_16S_ITS2/blob/main/BEL_16S_outputs/3_physeq_taxa16S.ipynb)
4. [Processing 16S phyloseq object](https://github.com/cdesouza02/BEL_16S_ITS2/blob/main/BEL_16S_outputs/4_ps_processing.ipynb)
  a. rarefation curves
  b. filtering thresholds
  c. 16S alpha diversity plots
5. Analysis of 16S phyloseq obejct (plots)
   a. How are the microbiomes of [individual coral colonies changing overtime and/or in response to stress?](https://github.com/cdesouza02/BEL_16S_ITS2/blob/main/BEL_16S_outputs/5_colonies_overtime.ipynb)
   b. What [variables are driving clustering](https://github.com/cdesouza02/BEL_16S_ITS2/blob/main/BEL_16S_outputs/5_ordination_stats.ipynb) of microbiome populations?
   c. How are different [coral species changing their microbiome interactions over time and in response to stress?](https://github.com/cdesouza02/BEL_16S_ITS2/blob/main/BEL_16S_outputs/5_species.ipynb)
6. [Focused P. asteroides analysis](https://github.com/cdesouza02/BEL_16S_ITS2/blob/main/BEL_16S_outputs/6_p_asteroides.ipynb.ipynb)
7. [Usable 16S phyloseq object](ps_16S.rds)

### BEL_ITS2_outputs
1. [Make submission sheet to send ITS2 sequences to symportal](https://github.com/cdesouza02/BEL_16S_ITS2/blob/main/BEL_ITS2_outputs/1_symportal.ipynb)
2. [Phyloseq object made from symportal run with all samples](https://github.com/cdesouza02/BEL_16S_ITS2/blob/main/BEL_ITS2_outputs/2_physeq_all_its2.ipynb)
3. [Processing ITS2 phyloseq object](https://github.com/cdesouza02/BEL_16S_ITS2/blob/main/BEL_ITS2_outputs/3_its2_ps_processing.ipynb)
4. Analysis of ITS2 phyloseq object
   a. Does exposure to stress change [the symbiont population (abundance or community)?](https://github.com/cdesouza02/BEL_16S_ITS2/blob/main/BEL_ITS2_outputs/4_Health_status.ipynb)
   b. [Assigning ITS2 color based on clade](https://github.com/cdesouza02/BEL_16S_ITS2/blob/main/BEL_ITS2_outputs/4_clades.ipynb)
   c. [Are the dominant symbionts in individual colonies changing overtime or in response to change?](https://github.com/cdesouza02/BEL_16S_ITS2/blob/main/BEL_ITS2_outputs/4_its2_colonies.ipynb)
   d. What [variables outside of species are driving clustering of symbiont communities?](https://github.com/cdesouza02/BEL_16S_ITS2/blob/main/BEL_ITS2_outputs/4_its2_ordination.ipynb)
5. [Usable ITS2 phyloseq object](ps_all_ITS2.rds)
6. Individual symportal runs
  a. physeq_individ_runs.ipynb
  b. ps_ITS2.rds

   
   
#### ITS2 questions and references
- Do coral species at the same reefs host different symbionts?
   - Siderastrea sidera has been desribed as stress tolerant when compared to other species, do we see anything in this data that corroborates that (Jones et.al 2025, Loya et al. 2001)?
   - Is there any evidence of "environmental memory" when face repeated heatwaves (Brown and Barott 2022; Hackerott et al. 2021)?

- Does bleaching impact the algal symbiont alpha diversity? 
    - Do bleached colonies recover after periods of stress?



## Labwork_plan
- **checking_files.ipynb** = checking number of reads from different dada2 runs
#### CBC_Extraction_Plan
-planning out sample extractions for coral holobiont project
#### CBC_PCR_Plan
- **.csv** = lab notes from each PCR preformed
- **16S_double_band.ipynb** = testing number of reads with 16S double banding
##### creating pcr metadata sheet
- **BEL_immune.ipynb**= samples from NSF Rapid
- **BEL_nonimmune.ipynb**= including historical coral samples

## qc_reports
- multiqc .html files for each sequencing run
- split by forward and reverse reads

## scripts
- bash and R scripts used for coral holobiont analysis
- named by the program being run in the script
  
## **sample_062026_PCR.csv** is the metadata for analysis
