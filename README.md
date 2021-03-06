# ChIPSeq-Analysis
I have used Snakemake based pipeline to analyse ChIP-seq data.

I have modified scripts from https://github.com/crazyhottommy/pyflow-ChIPseq to make snakefile run on the cluster at our facility.

Useful Resources:  
snakemake:  
https://hpc-carpentry.github.io/hpc-python/11-snakemake-intro/

snakemake + ChIP-seq:  
https://snakemake.readthedocs.io/en/stable/tutorial/basics.html
https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html
https://molb7621.github.io/workshop/Classes/snakemake-tutorial.html

Helpful github repositories:  
https://github.com/crazyhottommy/pyflow-ChIPseq
https://github.com/ENCODE-DCC
https://gist.github.com/davfre/8596159

ChIP-seq analysis:  
http://biocluster.ucr.edu/~rkaundal/workshops/R_feb2016/ChIPseq/ChIPseq.html
https://ngschool.eu/node/51
https://github.com/shenlab-sinai/chip-seq_preprocess
https://github.com/mforde84/ENCODE_TF_ChIP_pipeline
https://www.ebi.ac.uk/training/online/course/ebi-next-generation-sequencing-practical-course/gene-regulation/chip-seq-analysis
http://www.biologie.ens.fr/~mthomas/other/chip-seq-training/

Scripts from Marzi et al:   
https://epigenetics.essex.ac.uk/AD_H3K27ac/   
https://epigenetics.essex.ac.uk/AD_H3K27ac/code/index.html


Practical Guidelines for the Comprehensive Analysis of ChIP-seq Data  
http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003326

## Next Steps
1. phantom peaks
2. ENCODE metrics
3. macs2 peakcalling with different options eg: --nomodel, --broad
4. differential binding analyses: 
https://github.com/ying-w/chipseq-compare/tree/master/edgeR-DiffBind-DBChIP
https://github.com/crazyhottommy/ChIP-seq-analysis/blob/master/part3.1_Differential_binding_DiffBind_lib_size.md
DiffBind: http://bioconductor.org/packages/release/bioc/html/DiffBind.html    
chromswitch: https://bioconductor.org/packages/chromswitch    
https://academic.oup.com/bioinformatics/article/34/13/2286/4846890    
SCIDDO: https://github.com/ptrebert/sciddo    
https://www.biorxiv.org/content/early/2018/10/13/441766   

