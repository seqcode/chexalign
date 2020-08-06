# To create Fig1, run the following scripts. Jar, ChIP-exo data, and sacCer3.info can be found in example file.

java -Xmx8G -jar chexalign.public.jar --geninfo sacCer3.info --cpoints saccharomyces_cerevisiae_tRNA_gene_noM_distsort.tss --exptRap1 Rap1.bam --exptHmo1 Hmo1.bam --exptSfp1 Sfp1.bam --exptIfh1 Ifh1.bam --exptFhl1 Fhl1.bam --ctrl Control.bam --format BAM --out rp-127rpgs --gap 200 --cwin 1400 > rp-127rpgs.out

#Figure 1A and 1B
for mat in rp-127rpgs/intermediate-results/rp-127rpgs_*.mat;do
python makecomposite.py $mat rc 
python makeheatmap_2color.py $mat 20 rc
done

#Figure 1C
java -cp chexalign.public.jar org.seqcode.projects.chexalign.alignment.AlignmentPostAnalysis --geninfo sacCer3.info --mback sacCer3.back --mthres 0.1 --posf rp-127rpgs/intermediate-results/rp-127rpgs_original.regions --notstranded --out rp-127rpgs_original_rap1motif --seq sacCer3.fa --motif rap1.motif
python makemotifplot.py rp-127rpgs_original_rap1motif.out blue

#Figure 1D
java -cp chexalign.public.jar org.seqcode.projects.chexalign.alignment.AlignmentPostAnalysis --geninfo sacCer3.info --mback sacCer3.back --mthres 0.1 --posf rp-127rpgs/intermediate-results/rp-127rpgs_aligned.regions --notstranded --out rp-127rpgs_aligned_rap1motif --seq sacCer3.fa --motif rap1.motif

python makemotifplot.py rp-127rpgs_aligned_rap1motif.out blue

#Figure 1E (left and middle panel)
python plotStrandSeparateCompositeMultiExpt.py rp-127rpgs/rp-127rpgs_composite rp-127rpgs_composite 550 50 normalize
python plotStrandSeparateCompositeMultiExpt.py rp-127rpgs/rp-127rpgs_aligned_composite rp-127rpgs_aligned_composite 300 300 normalize

#To create Figure 1E (right panel)
java -Xmx8G -jar chexalign.public.jar --geninfo sacCer3.info --cpoints rap1-at-127RPGs.points --exptRap1 Rap1.bam --exptHmo1 Hmo1.bam --exptSfp1 Sfp1.bam --exptIfh1 Ifh1.bam --exptFhl1 Fhl1.bam --ctrl Control.bam --format BAM --out rp-127rpgs-rap1motif --gap 200 --cwin 1400 > rp-127rpgs-rap1motif.out

python plotStrandSeparateCompositeMultiExpt.py rp-127rpgs-rap1motif/rp-127rpgs-rap1motif_composite rp-127rpgs-rap1motif_composite 200 400 normalize
