# seda_DNA_Corsican_wreck_comparison
seda_DNA_Corsican_wreck_comparison

### Installing pipeline :

First, open your terminal. Then, run these two command lines :

    cd -place_in_your_local_computer
    git clone https://github.com/PLStenger/seda_DNA_Corsican_wreck_comparison.git

### Update the pipeline in local by :

    git pull

### References :

https://github.com/ZimmermannHH/BeringSea_shotgun_sequencing/


### Running :

    # bash scripts by (example)
    sbatch 10_mapdamage_cascade_species.sh
    # R scripts by (example)
    Rscript 11_compare_mapdamage_databases.R

    sbatch \
      --export=ALL,REF="/home/plstenge/seda_DNA_Corsican_wreck_comparison/10_mapdamage_cascade/runs/run_20261103_150905_9459089/references/Triticum_aestivum_IWGSC_CS_RefSeq_v2.1/reference.fa",SAMPLES_STR="sed6 sed8" \
      /home/plstenge/seda_DNA_Corsican_wreck_comparison/00_scripts/12_triticum_genome_ABD.sh \
      /home/plstenge/seda_DNA_Corsican_wreck_comparison/11_triticum_ITS_cascade/runs/run_20261112_115451_9952692
