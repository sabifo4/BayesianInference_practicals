# Data formatting

Before proceeding with timetree inference, we need to make sure that:

1. The alignment file is in PHYLIP format and easy to read (i.e., ideally one sequence per line).
2. The tree file is in Newick format.

You already have the raw data available in [directory  `00_raw_data`](00_raw_data). When you finish this session, you can [read the step-by-step tutorial in the `raw_data` directory](../raw_data/README.md) to see how the data have been collected and processed. But for now, let's focus on formatting this data for the subsequent analyses!

## Alignment file

If you open [the alignment file](00_raw_data/cytb_aln.fasta), you will see that each aligned sequence is already in a unique line, which makes it easier to convert into PHYLIP format. If that was not the case, remember that you can always use various tools such as [this in-house PERL script called `one_line_fasta.pl`](../src/one_line_fasta.pl) so that FASTA files in which the sequences for each taxon are split into more than one line can be converted into FASTA files in which sequences are written in one line.

Now, we will just run [another in-house PERL script called `FASTAtoPHYL.pl`](../src/FASTAtoPHYL.pl) to generate an alignment in PHYLIP format:

```sh
# Move to `00_data_formatting/00_raw_data`
# If you are not in this directory yet,
# please move to this directory and run the
# following commands
aln_name=`ls *aln.fasta`
a_noext=$( echo $aln_name | sed 's/\.fasta//' )
num=$( grep '>' $aln_name | wc -l )
len=$( sed -n '2,2p' $aln_name | sed 's/\r//' | sed 's/\n//' | wc -L )
perl ../../src/FASTAtoPHYL.pl $aln_name $num $len 
# Create a directory for input data for `MCMCtree`
mkdir ../01_inp_data
mv $a_noext.phy ../01_inp_data/aln.phy
```

You will now see a new directory called `01_inp_data` inside the directory [`00_data_formatting`](README.md). If you navigate to this newly created `01_inp_data` directory, you will find the alignment in PHYLIP format (i.e., the input file we need!). You will also find a log file called `log_lenseq.txt` inside the `00_raw_data` directory where you can read how many taxa were parsed and the length of each sequence.

Now, we can also partition our dataset according to a codon partitioning scheme. As it is widely known, first-codon positions and second-codon positions within a gene are expected to evolve in a more similar manner when compared to third-codon positions. The latter tend to evolve much faster than first- and second-codon positions, and thus have a much higher transition:transversion bias due to the fact that most of transversions are nonsynonymous for the third-codon positions ([Bofkin and Goldman, 2006](https://academic.oup.com/mbe/article/24/2/513/1150702)). In that way, we will evaluate the impact that codon partitioning scheme (CP scheme) can have on divergence times estimation by using two alignments: one with all the codon positions and another with only first- and second-codon positions (12CP).

We will use [the `fasta-phylip-partitions` pipeline](https://github.com/sabifo4/fasta-phylip-partitions), a tool I wrote some time ago, to generate CP-partitioned alignments. To make things easy, I have compressed this tool and saved it in the [`src`](../src/) directory. You just need to run the next commands to partition the dataset:

```sh
# Run from `src` to uncompress the file
tar -xvf fasta-phylip-partitions.tar.gz
chmod 775 fasta-phylip-partitions/src/*sh
chmod 775 fasta-phylip-partitions/src/Tools/*
# Now, change to `00_raw_data/01_inp_data` and run the following
# code to prepare the input data for the pipeline
cd ../00_data_formatting/01_inp_data/
mkdir part_aln
cd part_aln
grep '>' ../../00_raw_data/cytb_aln.fasta | sed 's/>//' > species_names.txt
cp ../../00_raw_data/cytb_aln.fasta .
# Now, run the pipeline. In essence, the first argument is the current 
# directory ("."), the second argument the tag ID for the job
# ("cytb"), and then a tag specifying if the alignment 
# needs to be partitioned into CPs, which we do want, and hence use 
# "partY".
# If you do not have this pipeline on your system, you can
# run the code below as it is using the source code that has already
# been added in the `src` directory. For more details about 
# how to use this pipeline, you can find this in the
# following link: https://github.com/sabifo4/fasta-phylip-partitions/blob/main/README.md
# NOTE: If you are running this code in a directory that you have synched to Dropbox or another 
# cloud storage and you have issues, just move the folder out of the synched directory and run the 
# command below again
../../../src/fasta-phylip-partitions/src/Run_tasks.sh . cytb partY
# Move alignment with only 1st+2nd CPs to main dir in `01_inp_data`
cp phylip_format/02_concatenated_alignments/part12/part12_cytb_concat.aln ../aln_12CP.phy
```

Now, we have generated the two alignments that we will use in our timetree inference analysis: one with all the codon positions and the other one with only the first and the second codon positions.

You are now ready to parse the tree file!

## Tree file

The calibrated input tree needs to be in PHYLIP format too, and we will also need to calibrate it! We will first remove the branch lengths and generate a PHYLIP formatted file with only the tree topology. Then, we will use a template file that has already been prepared to calibrate the tree topology with two different sets of calibrations.

First, we will generate the uncalibrated tree topology by removing the branch lengths:

```sh
# Move to directory `00_raw_data`
# Then, run the following commands
printf "8  1\n" > ../01_inp_data/tree_uncalib.tree
sed 's/\:[0-9]*\.[0-9]*//g' cytb_rooted_bl.tree  >> ../01_inp_data/tree_uncalib.tree
```

Now, we will calibrate the tree topology with two different sets of calibrations:

* **Fossil record (set #1)**: we will use two calibrations based on the fossil record:
  * **Caniformia**: given that we do not have a specific calibration for arctoids, we will calibrate the root of our phylogeny with the divergence of canids (dogs) and arctoids (musteloids, ursids, pinnipeds). The minimum age is 37.71 Ma and the maximum age is 66.09 Ma. The fossil taxon and specimen is *Hesperocyon gregarius* (SMNH P1899.6; 54) from the Cypress Hills Formation, Duchesnian NALMA, Lac Pelletier local fauna, Saskatchewan. In `MCMCtree` format and using a time unit of 100 Ma: `B(0.3771,0.6609)`.
  * **Crown group Ursinae**: the oldest first appearance of a crown group bear in our dataset is *Ursus americanus* around 1.84 Ma. We will use this age as a minimum constrain for that clade. In `MCMCtree` format and using a time unit of 100 Ma: `L(0.0184)`.
* **Secondary calibrations (set #2)**: instead of relying on the fossil record, we can also define our time constrains for various nodes in our tree topology based on the results obtained by previous molecular clock-dating analyses, an approach known as "secondary calibrations". In order to define this type of calibrations, we will use the results obtained by [Álvarez-Carretero et al., 2022](https://doi.org/10.1038/s41586-021-04341-1) when dating the "Laurasiatheria the rest" subtree given that it includes the taxa we have in our dataset. Specifically, we will use the divergence times estimated in the aforementioned study as secondary calibrations to constrain the ages of the two nodes calibrated using the fossil record. The resources we will be using are the following:
  * `out.txt`: this file is output by `MCMCtree` and contains various tree topologies in Newick format before printing out the estimated divergence times. As part of the materials I used in [Álvarez-Carretero et al., 2022](https://doi.org/10.1038/s41586-021-04341-1), I host a [publicly available GitHub repository](https://github.com/sabifo4/mammals_dating/) that can be used to reproduce our results. At the end of [the `README.md` file that I wrote to explain how to carry out the analyses for all the subtrees](https://github.com/sabifo4/mammals_dating/tree/main/02_SeqBayes_S2), users are provided with a link to [all the results for each subtree](https://www.dropbox.com/s/6gb746s7usitmcl/SeqBayesS2_MCMCtree.zip?dl=0). This compressed file is around ~11G, and the output file is too large to upload on GitHub. If you download the compressed file, you can find this specific output file if you follow the next path: `Laurasiatheria_therest/7/mcmc_files_filt/out.txt`. In a nutshell, we are keen on the following section in the output file `out.txt`:

    ```text
    # constant sites:      0 (0.00%)
    Species tree for FigTree.
    ((((((((((((((((((((((((1_aonyx_cinerea, 2_lutrogale_perspicillata) 683 , 3_aonyx_capensis) 682 , (4_lutra_lutra, 5_lutra_sumatrana) 684 ) 681 ,
    ((((((((((((((((((((((((aonyx_cinerea: 0.027567, lutrogale_perspicillata: 0.027567) [&95%={0.0172704, 0.039762}]: 0.026800, aonyx_capensis: 0.05
    ```

    The first Newick tree includes the node labels that `MCMCtree` has used during the analysis. I have saved this tree topology as [`BSS_L_therest_nodes.tree`](calibs/BSS_L_therest_nodes.tree). If you open this file with `FigTree`, you will see that the two nodes we are keen on calibrating are labelled as 782 and 787, which are also displayed in [figure `BSS_nodes_arctoidea.png`](calibs/BSS_L_therest_node_nums.png).

  * [Summary tsv file with time estimates per node](https://github.com/sabifo4/mammals_dating/blob/main/02_SeqBayes_S2/02_MCMCtree_updcrh/00_ESS_and_chain_convergence/out_data/Laurasiatheria_therest_filt/Laurasiatheria_therest_filt_all_mean_est.tsv): this file is output by an R script I wrote as part of the pipeline I used in [Álvarez-Carretero et al., 2022](https://doi.org/10.1038/s41586-021-04341-1). The first column can be used to identify specific nodes which, if calibrated, include the name of the calibration. We can look for nodes 782 and 784, which will be identified as `t_n782` and `t_n784`. We can see that the mean divergence time and the mean CIs for these nodes are the following:

    ```text
                Mean_time   Mean_qlow   Mean_qup
    t_n782    23.642      20.16       26.6
    t_n784     8.955       6.99       10.86
    ```

  Considering the results published in [Álvarez-Carretero et al., 2022](https://doi.org/10.1038/s41586-021-04341-1), the secondary calibrations will be the following:
  * **Root (Arctoidea)**: minimum age of 20.16 Ma and maximum age of 26.6 Ma. In `MCMCtree` format and using a time unit of 100 Ma: `B(0.2016,0.266)`.
  * **Internal node (Ursinae)**: minimum age of 6.99 Ma and maximum age of 10.86 Ma. In `MCMCtree` format and using a time unit of 100 Ma: `B(0.0699,0.1086)`.

  You can also find this tsv file saved in the `calibs` directory as [`BSS_L_therest_filt_all_mean_est.tsv`](calibs/BSS_L_therest_filt_all_mean_est.tsv).

Now, we can run the R script [`Include_calibrations.R`](scripts/Include_calibrations.R) to generate the calibrated tree files, which will be saved in the `01_inp_data` directory. Note that this script requires some input files that we have already prepared to save some time:

* **Tree topology with labels**: the [`cytb_calibnames.tree` file](calibs/cytb_calibnames.tree) contains the tree topology with two labels: `[ROOT]` and `[INT]`. These labels are replaced with the corresponding calibrations in `MCMCtree` format by the R script (see below). It is important to note that the labels need to be written within square brackets and that the last line is a blank line!
* **Matching files**: the [`Calib_converter_1.txt`](calibs/Calib_converter_1.txt) and [`Calib_converter_2.txt`](calibs/Calib_converter_2.txt) files include the name of the labels used to flag the nodes to be calibrated in the tree topology above (i.e., `ROOT` and `INT`), followed by a pipe delimiter `|`, and then the calibration in `MCMCtree` format within single quotation marks. The R script will use this matching file to replace the labels in the tree topology with the corresponding calibrations. It is important that the last line of these files is a blank line!

Once you run the R script [`Include_calibrations.R`](scripts/Include_calibrations.R), you will see two new files in the `01_inp_data` directory, each with the corresponding calibrations: `tree_fosscal.tree` (i.e., fossil calibrations) and `tree_seccal.tree` (i.e., secondary calibrations).

Now, we can move onto the next step: [calculating the Hessian and the gradient with `PAML` software!](../01_PAML/00_Hessian/README.md)
