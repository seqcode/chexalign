# To create Fig1, run the following scripts
# Note that you might need to move jar file and ChIP-exo datasets.
java -Xmx8G -jar chexalign.public.jar --geninfo sacCer3.info --cpoints saccharomyces_cerevisiae_tRNA_gene_noM_distsort.tss --exptRap1 Rap1.bam --exptHmo1 Hmo1.bam --exptSfp1 Sfp1.bam --exptIfh1 Ifh1.bam --exptFhl1 Fhl1.bam --ctrl Control.bam --format BAM --out rp-127rpgs --gap 100 --cwin 1400 > rp-127rpgs.out

for mat in rp-127rpgs/intermediate-results/rp-127rpgs_*.mat;do
python makecomposite.py $mat rc 
python makeheatmap_2color.py $mat 20 rc
done
