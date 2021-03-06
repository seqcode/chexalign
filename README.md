[![Anaconda-Server Badge](https://anaconda.org/bioconda/chexalign/badges/installer/conda.svg)](https://anaconda.org/bioconda/chexalign) 
[![Anaconda-Server Badge](https://anaconda.org/bioconda/chexalign/badges/downloads.svg)](https://anaconda.org/bioconda/chexalign)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/chexalign/badges/platforms.svg
)](https://anaconda.org/bioconda/chexalign)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/chexalign/badges/license.svg)](https://anaconda.org/bioconda/chexalign)

# ChExAlign
ChExAlign: alignment and quantification of ChIP-exo crosslinking patterns

ChExAlign is a computational framework that aligns ChIP-exo crosslinking patterns from multiple proteins across a set of regulatory regions, and which detects and quantifies protein-DNA crosslinking events within the aligned profiles. The output of the alignment approach is a set of composite profiles that represent the crosslinking signatures of the complex across analyzed regulatory regions. We then use a probabilistic mixture model to deconvolve individual crosslinking events within the aligned ChIP-exo profiles, enabling consistent measurements of protein-DNA crosslinking strengths across multiple proteins.

Citation:
--------------
N Yamada, MJ Rossi, N Farrell, BF Pugh, S Mahony .“Alignment and quantification of ChIP-exo crosslinking patterns reveal the spatial organization of protein-DNA complexes”. bioRxiv. [doi:10.1101/868604](https://doi.org/10.1101/868604).

Downloading Executables
--------------
 * ChExAlign version 0.12 (2020-08-05): [JAR](http://lugh.bmb.psu.edu/software/chexalign/chexalign_v0.12.jar)
 * ChExAlign version 0.11 (2020-04-06): [JAR](http://lugh.bmb.psu.edu/software/chexalign/chexalign_v0.11.jar)
 * ChExAlign version 0.1 (2019-12-10): [JAR](http://lugh.bmb.psu.edu/software/chexalign/chexalign_v0.1.jar)

Building from Source
--------------
If you want to build the code yourself, you will need to first download and build the seqcode-core library (https://github.com/seqcode/seqcode-core) and add its build/classes and lib directories to your CLASSPATH.

Dependencies:
--------------
1. ChExAlign requires Java 8+. 
2. ChExAlign loads all data to memory, so you will need a lot of available memory if you are running analysis over many conditions or large datasets.

Running ChExAlign
--------------
Running from a jar file:

```{r, engine='sh', count_lines}
java -Xmx10G -jar chexalign.jar <options - see below>
```

In the above, the “-Xmx10G” argument tells java to use up to 10GB of memory. If you have installed source code from github, and if all classes are in your CLASSPATH, you can run ChExMix as follows:

```{r, engine='sh', count_lines}
java -Xmx10G org.seqcode.projects.chexalign.ChExAlign <options - see below>
```

Options (Required/important options are in __bold__)

__General__:

 * --__out__ \<prefix>: Output file prefix. All output will be put into a directory with the prefix name. 

__Specifying the Genome__:

 * --__geninfo__ \<genome info file\>: This file should list the lengths of all chromosomes on separate lines using the format chrName\<tab\>chrLength. You can generate a suitable file from UCSC 2bit format genomes using the UCSC utility “twoBitInfo”. The chromosome names should be exactly the same as those used in your input list of genomic regions. 
   
      The genome info files for some UCSC genome versions:  
      | [hg18](http://lugh.bmb.psu.edu/software/multigps/support/hg18.info) | [hg19](http://lugh.bmb.psu.edu/software/multigps/support/hg19.info) | [hg38](http://lugh.bmb.psu.edu/software/multigps/support/hg38.info) | [mm8](http://lugh.bmb.psu.edu/software/multigps/support/mm8.info) | [mm9](http://lugh.bmb.psu.edu/software/multigps/support/mm9.info) | [mm10](http://lugh.bmb.psu.edu/software/multigps/support/mm10.info) | [rn4](http://lugh.bmb.psu.edu/software/multigps/support/rn4.info) | [rn5](http://lugh.bmb.psu.edu/software/multigps/support/rn5.info) | [danRer6](http://lugh.bmb.psu.edu/software/multigps/support/danRer6.info) | [ce10](http://lugh.bmb.psu.edu/software/multigps/support/ce10.info) | [dm3](http://lugh.bmb.psu.edu/software/multigps/support/dm3.info) | [sacCer2](http://lugh.bmb.psu.edu/software/multigps/support/sacCer2.info) | [sacCer3](http://lugh.bmb.psu.edu/software/multigps/support/sacCer3.info) |

__Loading Data__:

 * --__exptCONDNAME-REPNAME__ \<file\>: Defines a file containing reads from a signal experiment. Replace CONDNAME and REPNAME with appropriate condition and replicate labels.
 * --__ctrlCONDNAME-REPNAME__ \<file\>: Optional arguments. Defines a file containing reads from a control experiment (such as mock ChIP, IgG or Input control). Replace CONDNAME and REPNAME with appropriate labels to match a signal experiment (i.e. to tell ChExAlign which condition/replicate this is a control for). If you leave out a REPNAME, this file will be used as a control for all replicates of CONDNAME.  
 * --__format__ \<SAM/BAM/BED/IDX\>: Format of data files. All files must be the same format if specifying experiments on the command line. Supported formats are SAM/BAM, BED, and IDX index files.
 
Instead of using the above options to specify each and every ChIP-exo data file on the command-line, you can instead use a design file:
 
 * --__design__ \<file\>: A file that specifies the data files and their condition/replicate relationships. See [here](http://lugh.bmb.psu.edu/software/multigps/example.design) for an example design file. The file should be formatted to contain the following pieces of information for each data file, in this order and tab-separated:
  
    * File name
    * Label stating if this experiment is “signal” or “control”
    * File format (SAM/BAM/BED/IDX) – mixtures of formats are allowed in design files
    * Condition name
    * Replicate name (optional for control experiments – if used, the control will only be used for the corresponding named signal replicate)
    
 Limits on how many reads can have their 5′ end at the same position in each replicate. Please use the options if your experiments contain many PCR duplicates:

 * --readfilter: Flag to turn on filtering reads, recommended for highly duplicated experiments. It estimates a global per-base limit from a Poisson distribution parameterized by the number of reads divided by the number of mappable bases in the genome. The per-base limit is set as the count corresponding to the 10^-7 probability level from the Poisson. Default = no read filter
 * --fixedpb \<value\>: Fixed per-base limit.
 * --poissongausspb \<value\>: Filter per base using a Poisson threshold parameterized by a local Gaussian sliding window (i.e. look at neighboring positions to decide what the per-base limit should be). 
 * --nonunique: Flag to use non-unique reads. 
 * --mappability \<value\>: Fraction of the genome that is mappable for these experiments. Default=0.8.
 * --nocache: Flag to turn off caching of the entire set of experiments (i.e. run slower with less memory)
 
 __Scaling Data__:
 
 If you do not provide a control experiment, data scaling will not be performed. 
 
 * --noscaling: Flag to turn off auto estimation of signal vs control scaling factor.
 * --medianscale: Flag to use scaling by median ratio of binned tag counts. Default = scaling by NCIS.
 * --regressionscale: Flag to use scaling by regression on binned tag counts. Default = scaling by NCIS.
 * --sesscale: Flag to use scaling by SES (Diaz, et al. Stat Appl Genet Mol Biol. 2012).
 * --fixedscaling \<scaling factor\>: Multiply control counts by total tag count ratio and then by this factor. Default: scaling by NCIS.
 * --scalewin \<window size\>: Window size for estimating scaling ratios. Default is 10Kbp. Use something much smaller if scaling via SES (e.g. 200bp).
 * --plotscaling: Flag to plot diagnostic information for the chosen scaling method.
 
__Running ChExAlign__:

 * --__bed__ \<file\>: Bed format file of genomic positions to perform alignment. This can be genomic positions from peak results, or genomic annotations if you know that target proteins bind there.

Instead of using the bed file to specify genomic positions, you can instead use cpoints or points file:

 * --cpoints \<file\>: File of genomic positions in point format.
 * --points \<file\>: File of genomic positions in stranded point format.
 
 * --cwin \<int\>: Window size for analyzing read profiles (default=400).
  
__Alining Crosslinking Patterns__:

 * --gap \<value\>: Gap open penalty (default=100). Use 1000 for non-gapped alignment.
 * --extscaling \<value\>: Gap extension scaling factor (default=0.1). Increasing this parameter results in greater gap extension penalty.
 * --sort: Flag to output per region alignment by the order of genomic position input file (default=off).

__Quantifying Crosslinking Events__:

 * --r \<int\>: Max. model update rounds (default=3).
 * --xlsigma \<value\>: Crosslinking component sigma (default=6)
 * --noposprior: Flag to turn off inter-experiment positional prior (default=on).

Example
--------------
This example runs ChExAlign v0.11 on ribosomal protein gene (RPG) datasets (NCBI Sequence Read Archive under accession number SRP041518) presented in Figure 1 from the paper. The bam files include ChIP-exo data for Rap1, Hmo1, Sfp1, Ifh1, Fhl1, and control experiment. The version of ChExAlign and all files except for python scripts required to run this analysis are in this file: [chexalign-yeast-example.tar.gz](http://lugh.bmb.psu.edu/software/chexalign/examples/chexalign-yeast-example.tar.gz).

Let’s demonstrate how ChExAlign works. Note that this example only uses the top 20 RPG sites to save time.

```{r, engine='sh', count_lines}
java -Xmx8G -jar chexalign.public.jar --geninfo sacCer3.info --cpoints rp-127rpgs.20.spoints --exptRap1 Rap1.bam --exptHmo1 Hmo1.bam --exptSfp1 Sfp1.bam --exptIfh1 Ifh1.bam --exptFhl1 Fhl1.bam --ctrl Control.bam --format BAM --out chexalign-test-results --gap 200 --nosort --cwin 1000 > chexalign-test-results.out 2>&1
```

Expected results are named chexalign-test and can be found in the same file.

To visualize your results, you need to download python scripts located [here](https://github.com/seqcode/chexalign/tree/master/pythonscripts) and put these scripts under `chexalign-yeast-example` directory. Run [this](https://github.com/seqcode/chexalign/blob/master/pythonscripts/plotStrandSeparateCompositeMultiExpt.py) to make an alignment figure.

```{r, engine='sh', count_lines}
python plotStrandSeparateCompositeMultiExpt.py flist.txt 500 100 normalize
```

`flist.txt` should contain a list of ChExAlign output files. Each line should indicate a single file path as following:

```{r, engine='sh', count_lines}
chexalign-test-results/chexalign-test-results_composite.Rap1.txt
chexalign-test-results/chexalign-test-results_composite.Sfp1.txt
chexalign-test-results/chexalign-test-results_composite.Ifh1.txt
```

Run [this](https://github.com/seqcode/chexalign/blob/master/pythonscripts/reduce.py) to make PCA and MDS plots.

```{r, engine='sh', count_lines}
python reduce.py inputflist.txt 
```

`inputflist.txt` should contain a list of ChExAlign output files similar to the following:

```{r, engine='sh', count_lines}
chexalign-test-results/chexalign-test-results_site-component-ML.Sfp1.txt
chexalign-test-results/chexalign-test-results_site-component-ML.Hmo1.txt
chexalign-test-results/chexalign-test-results_site-component-ML.Ifh1.txt
```

Output files
--------------


Contact
--------------
For queries, please contact Naomi (nuy11@psu.edu) or Shaun Mahony (mahony@psu.edu).

Major History:
-------------- 
Version 0.12 (2020-08-05): Fixed bugs on sequence generation using jar file. 
Version 0.11 (2020-04-06): Added classes to perform alignment using random regions and to output Pearson correlation scores from random input arrays. Support bed file as input. Re-organize output files.
Version 0.1 (2019-11-05): Initial release.
