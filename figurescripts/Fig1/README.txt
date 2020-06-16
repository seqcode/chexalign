# To create Fig1, run the following scripts
java -Xmx8G org.seqcode.projects.chexalign.ChExAlign --geninfo sacCer3.info --cpoints saccharomyces_cerevisiae_tRNA_gene_noM_distsort.tss --design rp.design --out rp-127rpgs --gap 100 --cwin 1400 > rp-127rpgs.out

for mat in rp-127rpgs/intermediate-results/rp-127rpgs_*.mat;do
python makecomposite.py $mat rc 
python makeheatmap_2color.py $mat 20 rc
done
