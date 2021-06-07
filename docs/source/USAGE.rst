.. title:: Usage

Modules
-----------------------------------------------------------

1. **General command line**

  .. code-block:: shell

	$ python shuffle.py [-h] {shuffle,plot-dis,plot-comp,compare} ...
	Shuffles SNPs and outputs homozygosity runs along contigs.

	positional arguments:
	  {shuffle,plot-dis,plot-comp,compare}
		shuffle             Shuffles SNPs and output closest distance of observed and shuffled data.
		plot-dis            Plot shuffled vs observed distributions based on shuffle command results.
		plot-comp           Compare heterozygosity of regions in and out of a bed file.
		compare             Compare heterozygosity of regions in and out of a bed file.

	optional arguments:
	  -h, --help            show this help message and exit
	  
	  
2. **Shuffle**

The shuffler will shuffle SNPs along contigs based on a .BED file. The variants inside regions defined by the .BED file will be shuffled in the span of this region (i.e. a CDS). The variants **outside** of the regions will be shuffled along the whole length of the contig not spanned by the regions. The shuffled variants distribution can be compared with the observed distribution to check whether CDS or other regions are denser or less dense in variants.

  .. code-block:: shell
  
	usage: shuffle.py shuffle [-h] [-s SAMPLE] [-lc LOW_COVERAGE]
		                      [-hc HIGH_COVERAGE] [-mq MIN_QUAL] [-p [PLOT]]
		                      [-cs CONTIG_SIZE]
		                      VCF FASTA BED OUT

	positional arguments:
	  VCF                   <STRING>
		                           A path to a the short variants calls (vcf file).
	  FASTA                 <STRING>
		                           A path to the reference genome (fasta file).
	  BED                   <STRING>
		                           A path to the regions locations (sorted and merged bed file).
	  OUT                   <STRING>
		                           An output directory path for the out files.

	optional arguments:
	  -h, --help            show this help message and exit
	  -s SAMPLE, --sample SAMPLE <STRING>
		                         Sample to use for short variants positions. Default: ['ancestor']
	  -lc LOW_COVERAGE, --low-coverage LOW_COVERAGE <INT>
		                           Minimum coverage required to consider a good variant. Default: [None]
	  -hc HIGH_COVERAGE, --high-coverage HIGH_COVERAGE <INT>
		                           Maximum coverage threshold to consider a good variant. Default: [None]
	  -mq MIN_QUAL, --min-qual MIN_QUAL <INT>
		                           Minimum QUAL required to consider a good variant. Default: [None]
	  -p [PLOT], --plot [PLOT]
		                           Output distribution plots. Default: True
	  -cs CONTIG_SIZE, --contig-size CONTIG_SIZE <INT>
		                           Minimum size of a contig to plot. Default: [500000]bp

Example output: `<sample>_DataFrame.csv`

+--------+-------------+------------------------------------------------------------------------+
| Column | Name        | Explanation                                                            |
+========+=============+========================================================================+
| 1      |             | index                                                                  |
+--------+-------------+------------------------------------------------------------------------+
| 2      | POS         | position along the contig                                              |
+--------+-------------+------------------------------------------------------------------------+
| 3      | SNP         | python boolean determining if site is a SNP or not                     |
+--------+-------------+------------------------------------------------------------------------+
| 4      | REGIONS     | python boolean determining if site was in a .BED defined region or not |
+--------+-------------+------------------------------------------------------------------------+
| 5      | shuffledSNP | python boolean determining if site became a SNP after shuffling        |
+--------+-------------+------------------------------------------------------------------------+
| 6      | CHROM       | name of the contig                                                     |
+--------+-------------+------------------------------------------------------------------------+

  .. code-block:: shell
  
	,POS,SNP,REGIONS,shuffledSNP,CHROM
	0,0,False,False,False,Chrom_3
	1,1,False,False,False,Chrom_3
	2,2,False,False,False,Chrom_3
	3,3,False,False,False,Chrom_3
	4,4,False,False,False,Chrom_3
	5,5,False,False,False,Chrom_3
	6,6,False,False,False,Chrom_3
	7,7,False,False,False,Chrom_3
	8,8,False,False,False,Chrom_3


Example graphics output: simulated data with genes under selection and large homozygous regions at extremities of contigs. The expected curve (orange) shows that heterozygosity should indeed be higher at telomeres while accounting for lower overall %heterozygosity in genes.
  .. image:: ../../gitexample-1.png
  
2. **Plot-dis**

Takes .CSV output of `shuffle.py shuffle` and outputs graphs. See output of the `shuffle` module.

  .. code-block:: shell

	usage: shuffle.py plot-dis [-h] [-o OUTPUT] [-co COLOR_OBSERVED]
		                       [-cs COLOR_SHUFFLED] [-wd WIDTH] [-hg HEIGHT]
		                       [-s SIZE] [-b BINS]
		                       CSV FASTA
		                       
	positional arguments:
	  CSV                   <STRING> A path to a the output of the shuffle command (csv file).
	  FASTA                 <STRING> A path to the reference genome (fasta file).
	  
	optional arguments:
	  -h, --help            show this help message and exit
	  -o OUTPUT, --output OUTPUT <STRING>
		                           Directory name to output. Default: path/to/cur_dir/['plots']
	  -co COLOR_OBSERVED, --color-observed COLOR_OBSERVED <STRING>
		                           A matplotlib valid color for the observed distribution. Default: ['dodgerblue']
	  -cs COLOR_SHUFFLED, --color-shuffled COLOR_SHUFFLED <STRING>
		                           A matplotlib valid color for the shuffled distribution. Default: ['orange']
	  -wd WIDTH, --width WIDTH <INT>
		                           Width of plot to output. Default: [14]
	  -hg HEIGHT, --height HEIGHT <INT>
		                           Height of plot to output. Default: [7]
	  -s SIZE, --size SIZE  <INT>
		                           Minimum size of a contig to plot. Default: [500000]bp
	  -b BINS, --bins BINS  <INT>
		                           Number of bins (windows) to make on each chromosome. Default: [501]bp

3. **Compare**

Takes a .VCF, a .FASTA and a .BED file and compares heterozygosity distribution inside and outside defined regions. This helps verifying that some regions have lower or higher densities of variants.

  .. code-block:: shell

	usage: shuffle.py compare [-h] [-s SAMPLE] [-o OUTPUT] [-ci COLOR_IN]
		                      [-co COLOR_OUT] [-wd WIDTH] [-hg HEIGHT]
		                      [-cs CONTIG_SIZE] [-rs REGION_SIZE] [-MS MAX_SIZE]
		                      [-ms MIN_SIZE] [-bn BIN_NUMBER] [-MH MAX_HET]
		                      [-mh MIN_HET] [-bhn BIN_HET_NUMBER]
		                      VCF FASTA BED

	positional arguments:
	  VCF                   <STRING>
		                             A path to a the short variants calls (vcf file).
	  FASTA                 <STRING>
		                             A path to the reference genome (fasta file).
	  BED                   <STRING>
		                             A path to the regions locations (sorted and merged bed file).

	optional arguments:
	  -h, --help            show this help message and exit
	  -s SAMPLE, --sample SAMPLE <STRING>
		                             Sample name to read in VCF. Default: ['ancestor']
	  -o OUTPUT, --output OUTPUT <STRING>
		                             Directory name to output. Default: path/to/cur_dir/['plots']
	  -ci COLOR_IN, --color-in COLOR_IN <STRING>
		                             A matplotlib valid color for the inside distribution. Default: ['dodgerblue']
	  -co COLOR_OUT, --color-out COLOR_OUT <STRING>
		                             A matplotlib valid color for the outside distribution. Default: ['orange']
	  -wd WIDTH, --width WIDTH <INT>
		                             Width of plot to output. Default: [14]
	  -hg HEIGHT, --height HEIGHT <INT>
		                             Height of plot to output. Default: [7]
	  -cs CONTIG_SIZE, --contig-size CONTIG_SIZE <INT>
		                             Minimum size of a contig to plot. Default: [500000]bp
	  -rs REGION_SIZE, --region-size REGION_SIZE <INT>
		                             Minimum size of a region to be considered in output plots. Default: [200]bp
	  -MS MAX_SIZE, --max-size MAX_SIZE <INT>
		                             Maximum size of a region to plot the histogram (regions in range [ms, MS]). Default: [500000]
	  -ms MIN_SIZE, --min-size MIN_SIZE <INT>
		                             Minimum size of a region to plot the histogram (regions in range [ms, MS]). Default: [200]
	  -bn BIN_NUMBER, --bin-number BIN_NUMBER <INT>
		                             Number of bins in histogram in range [ms, MS]. Default: [100]
	  -MH MAX_HET, --max-het MAX_HET <INT>
		                             Maximum heterozygosity range of a region to plot the histogram (in range [mh, MH]). Default: [10]
	  -mh MIN_HET, --min-het MIN_HET <INT>
		                            Minimum heterozygosity range of a region to plot the histogram (in range [mh, MH]). Default: [0]
	  -bhn BIN_HET_NUMBER, --bin-het-number BIN_HET_NUMBER <INT>
		                             Number of bins in histogram in range [mh, MH]. Default: [100]

Example output: `<sample>_Comparison_DataFrame.csv`

+--------+-------------+------------------------------------------------------------------------+
| Column | Name        | Explanation                                                            |
+========+=============+========================================================================+
| 1      | REGIONS     | index                                                                  |
+--------+-------------+------------------------------------------------------------------------+
| 2      | POS         | length of the region                                                   |
+--------+-------------+------------------------------------------------------------------------+
| 3      | REGIONS     | python boolean determining if the region is inside or not              |
+--------+-------------+------------------------------------------------------------------------+
| 4      | SNP         | integer number of variants found in this region                        |
+--------+-------------+------------------------------------------------------------------------+
| 5      | CHROM       | name of the contig                                                     |
+--------+-------------+------------------------------------------------------------------------+

  .. code-block:: shell

	REGIONS,POS,REGIONS,SNP,CHROM
	1,4391,False,38,Chrom_3
	2,1012,True,18,Chrom_3
	3,2802,False,31,Chrom_3
	4,1081,True,43,Chrom_3
	5,3863,False,149,Chrom_3
	6,2861,True,56,Chrom_3
	7,3356,False,109,Chrom_3
	8,14533,True,0,Chrom_3
	9,11302,False,62,Chrom_3
	
Example graphics output: histogram of heterozygosity inside and outside defined regions

  .. image:: ../../gitexample-2.png
  
3. **Plot-comp**

  .. code-block:: shell

	usage: shuffle.py plot-comp [-h] [-s SAMPLE] [-o OUTPUT] [-ci COLOR_IN]
		                        [-co COLOR_OUT] [-wd WIDTH] [-hg HEIGHT]
		                        [-cs CONTIG_SIZE] [-rs REGION_SIZE] [-MS MAX_SIZE]
		                        [-ms MIN_SIZE] [-bn BIN_NUMBER] [-MH MAX_HET]
		                        [-mh MIN_HET] [-bhn BIN_HET_NUMBER]
		                        CSV FASTA

	positional arguments:
	  CSV                   <STRING> A path to a the short variants calls (vcf
		                    file).
	  FASTA                 <STRING> A path to the reference genome (fasta file).

	optional arguments:
	  -h, --help            show this help message and exit
	  -s SAMPLE, --sample SAMPLE <STRING>
		                               Sample name (only for output file names). Default: ['unknown']
	  -o OUTPUT, --output OUTPUT <STRING>
		                               Directory name to output. Default: path/to/cur_dir/['compare_plots']
	  -ci COLOR_IN, --color-in COLOR_IN <STRING>
		                               A matplotlib valid color for the inside distribution. Default: ['dodgerblue']
	  -co COLOR_OUT, --color-out COLOR_OUT <STRING>
		                               A matplotlib valid color for the outside distribution. Default: ['orange']
	  -wd WIDTH, --width WIDTH <INT>
		                               Width of plot to output. Default: [14]
	  -hg HEIGHT, --height HEIGHT <INT>
		                               Height of plot to output. Default: [7]
	  -cs CONTIG_SIZE, --contig-size CONTIG_SIZE <INT>
		                               Minimum size of a contig to plot. Default: [500000]bp
	  -rs REGION_SIZE, --region-size REGION_SIZE <INT>
		                               Minimum size of a region to be considered in output plots. Default: [200]bp
	  -MS MAX_SIZE, --max-size MAX_SIZE <INT>
		                               Maximum size of a region to plot the histogram (regions in range [ms, MS]). Default: [500000]
	  -ms MIN_SIZE, --min-size MIN_SIZE <INT>
		                               Minimum size of a region to plot the histogram (regions in range [ms, MS]). Default: [200]
	  -bn BIN_NUMBER, --bin-number BIN_NUMBER <INT>
		                               Number of bins in histogram in range [ms, MS]. Default: [100]
	  -MH MAX_HET, --max-het MAX_HET <INT>
		                               Maximum heterozygosity range of a region to plot the histogram (in range [mh, MH]). Default: [10]
	  -mh MIN_HET, --min-het MIN_HET <INT>
		                               Minimum heterozygosity range of a region to plot the histogram (in range [mh, MH]). Default: [0]
	  -bhn BIN_HET_NUMBER, --bin-het-number BIN_HET_NUMBER <INT>
		                               Number of bins in histogram in range [mh, MH]. Default: [100]


