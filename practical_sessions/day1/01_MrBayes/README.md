# 1. Introduction 
For this tutorial, we assume that you have already installed `MrBayes` on your PC
or that, alternatively, you will be using the one installed on the cluster. 
If you want to use it on your PC but have had issues with its installation, please
follow the instructions on their [GitHub repository](https://github.com/NBISweden/MrBayes).

We will be diving this practical session into the following sections:   

   * Getting familiar with the commands in MrBayes   
   * Running MrBayes   
   * EXERCISES: Evaluating the effect that using different models can have on posterior parameter estimates   
   * Learning how to use Tracer for MCMC diagnostics   
   
Let's get started! 

# 2. Getting familiar with the nexus format 

## Input data 
First, we will be looking at the different parts the input file has so we can understand 
which commands we need to use to define our input data, i.e., our alignment.

In [`data`](data), you will find a nexus file: [`apes.nex`](data/apes.nex). Open this file with your preferred 
text editor as we will go through the main commands you see here. As a summary, we mention them below:   

   * `begin data`: this command specifies the beginning of the data block. More details [here](http://mrbayes.sourceforge.net/Help/begin.html).   
   * `dimensions`: after this command, you define the number of taxa (`ntax`) and the number of characters in your alignment (`nchar`).
   More details [here](http://mrbayes.sourceforge.net/Help/dimensions.html).
   * `format`: this command specifies the type of data of your alignment (`datatype`), e.g., dna, protein,rna, etc. You can also indicate 
   if it is interleaved or not. More details [here](http://mrbayes.sourceforge.net/Help/format.html).
   * `matrix`: after this command (normally next line), you need to include the alignment. Interleave format is also accepted.
   More details [here](http://mrbayes.sourceforge.net/Help/matrix.html).   
   * `end`: this command terminates the data block.   
   
>> **NOTES:**   
>> Add a `;` to terminate a command (i.e., `begin data;`). If you do not add a semicolon, `MrBayes` 
>> will understand that all the information/arguments you added correspond to the last command you 
>> wrote in the data block even though it is in a new line! Also, whatever information you add within 
>> square brackets is considered a comment, and hence `MrBayes` does not execute it. Use this to 
>> add useful comments in your nexus file so you remember what each block corresponds to:
>> E.g., `[ this is a comment and it is not run by MrBayes ]`.   

## `MrBayes` commands 
Once we have defined our input data, we can generate another nexus file where the commands will define 
the analysis that will be carried out by `MrBayes`. Open a new 
text file in your preferred text editor as we will create this nexus file together. Below, you will find 
a summary of the commands that we will go through:   

**BLOCK 1: Start `MrBayes` and read input data**:   
   * `being mrbayes`: this command initiates the block with instructions to run `MrBayes`.   
   * `log`: this command outputs to a file the screen messages. You can specify the file name (`filename`) and 
   if you want append or replace the file (`append/replace`). You start or stop this command by including the
   option `start` or `stop`, respectively. More details [here](http://mrbayes.sourceforge.net/Help/log.html).   
   * `execute`: this commands is used to read the input file with your data alignment, which you have previously generated.   
   * `outgroup`: use this command to write down the name of the taxon/taxa that are to be used as an outgroup.   

**BLOCK 2: Define analysis in `MrBayes`**:   
   * `charset`: this command is used to specify different character sets you have defined in your data block.
   More details [here](http://mrbayes.sourceforge.net/Help/charset.html).   
   * `partition`: this command uses the information you have passed to `charset` and the names you gave to each 
   set to specific a character partition. It follows the next format: `partition <name> = <num_parts>:<name_charset1>, ...,<name_chraset_n`.
      E.g.:   
	  ```
      # Example 1   
      charset gene1=1-500;   
	  charset gene2=501-1000;
	  partition by-gene=2: gene1, gene2;   
	  ```   
	  ```
      # Example 2   
      charset gene1_cp1=1-500\3;   
      charset gene1_cp2=2-500\3;   
      charset gene1_cp3=3-500\3;   
	  charset gene2=501-1000;
	  partition by-codpos=2: gene1_cp1, gene1_cp2, gene1_cp3, gene2;   
	  ```   
      More details [here](http://mrbayes.sourceforge.net/Help/partition.html).
   * `set`: this command is to be used if you have used the `partition` command as it "sets up" what you previously defined.   
   * `lset`: this command sets the parameters of the likelihood model. There are different options this 
   command can take, but we will focus on `nst`, `applyto`, and `rates`. More details [here](http://mrbayes.sourceforge.net/Help/lset.html).   
   * `prset`: use this command to set the priors for the phylogenetic model you want to use. There are many options that you 
   can pass to this command, but we will focus on `statefreqpr`, `shapepr`, and `revmatpr`. 
   More details [here](http://mrbayes.sourceforge.net/Help/prset.html).   
   * `unlink`: this command is to be used if you have used the `partition` command. As the name of this command says, it will 
   "unlink" model parameters across the data partitions you have defined. You can add type `all` or the specific name 
   of the partition/s for which you want to unlink the model parameters. By default, if the same parameter applies to different 
   partitions and if this parameter has the same prior, `MrBayes` will use a single value for this parameter. If you want to use different parameter 
   values for each partition you have established, then you need to use this command: it will "unlink" the model parameters so a specific 
   value will be inferred for each partition. If you use the command `link` instead of `unlink`, the opposite will occur.
   More details [here](http://mrbayes.sourceforge.net/Help/unlink.html).   
   * `mcmc`: this command is used to set up and start the MCMC analysis. There are different options that you can pass to this 
   command, but we will focus on `seed`, `ngen`, `nruns`, `nchains` (default is 4, 1 cold chain and 3 heated chains),
   `printfreq`, `samplefreq`, `diagnfreq`, `diagnstat`, `savebrlens`, and `filename`.
   More details [here](http://mrbayes.sourceforge.net/Help/mcmc.html).   

**BLOCK 3: Summarise trees and other model parameters (as many blocks as analysis you want to perform!)**:   
   * `sumt`: this command produces summary statistics for the trees that have been sampled during the MCMC. You can 
   specify the file name (`filename`) where you want the output to be written. By default, the burnin is established 
   to be 25% of the samples collected (you could modify this if required).
   More details [here](http://mrbayes.sourceforge.net/Help/sumt.html).
   * `sump`: this command prints the values that have been sampled for the model parameters during the MCMC. You can 
   specify the file name (`filename`) where you want the output to be written. By default, the burnin is established 
   to be 25% of the samples collected (you could modify this if required).
   More details [here](http://mrbayes.sourceforge.net/Help/sump.html).
   
**BLOCK 4: Stop `MrBayes` and end of nexus file**:   
   * `quit`: this command exits `MrBayes`.   
   * `end`: this command is used to indicate that the `MrBayes` block has come to an end.   
   
# 3. Run `MrBayes` 
Now, we will run `MrBayes` using the two nexus files that we have generated in the previous section. 

If you are a Windows user and want to use the executable file, we suggest you copy this file
to the directory where you have saved the nexus files and execute `MrBayes` from there. If you are running
`MrBayes` from the command line and you have properly installed the software, you can open a terminal 
from the directory where you have your nexus files and type `mb` to start `MrBayes`. In both cases, you 
should see the following screen: 

<p align="center">
  <img width="700" height="400" src="../../../figs/figs_MrBayes/fig1.png">
</p>

Now, we will run a live demonstration of how you can use the commands. 
Then, we will let `MrBayes` with the example data so you have time to go through the exercises!

## Exercise
You will split into different groups so you can come up with 2 different models under which you
could run `MrBayes` and analyse this example data set. Once you have decided which two models you
are going to run in `MrBayes`, set up your nexus file and then start the Bayesian inference analysis.

We will then gather together and each group will 
present the models chosen and explain what they expect to obtain once inference finishes. 

# 4. Analysing the `MrBayes` MCMC output

To analyse the MCMC output, we are going to use the program `Tracer`.
Run `Tracer` by double clicking on the `Tracer` icon. If you want to run `Tracer` from 
the command line, you can either set an alias in your `~/.bashrc` or `~/.bash_profile` 
file, run it from
the directory where the `jar` file can be found, or execute this `jar` file 
using an absolute or a relative path to the directory where the terminal is running.
Below, you can find a summary of how you would proceed for each of these three options 
when using the terminal:

```sh 
# Option 1: Run from the directory where 
# `tracer.jar` can be found. For instance, 
# if you open a terminal from
# `Tracer_v1.7.1/lib`, then:
java -jar tracer.jar


# Option 2: Use an absolute or a relative 
# path to execute the file. An example of 
# an absolute path is shown below:
java -jar Applications/Tracer_v1.7.1/lib/tracer.jar

# Option 3: Add an alias in your 
# `~/.bashrc` or `~/.bash_profile.
# Please replace "<path_to_Tracer>" with the absolute path
# to the location where you have saved `Tracer`
# Modify the name of the directory that you have unzipped
# if needed too (e.g., `Tracer_v1.7.1`)
alias tracer1.7.1='java -jar <path_to_Tracer>/Tracer_v1.7.1/lib/tracer.jar'
# Now, open a terminal from any location on 
# your PC and type the following command to 
# execute `Tracer`. Before you do this, please 
# make sure that you have X11 installed or, if 
# you are on Windows, that you have the Xming 
# Server running on your PC.
# If you are on Windows, you might need to run 
# the command `export DISPLAY=:0.0` before 
# you can execute `Tracer`
tracer1.7.1
```

By now, we expect the first analysis has finished so you can use the 
corresponding output file (i.e., the file which extension is `.p`) to analyse 
the sampled values gathered during the MCMC for each model parameter. 
If your runs have not finished, we will share with you an output file we obtained 
when running this analysis so you can proceed with the next tasks. 

Once you have started `Tracer`, you can select the `Import Trace File...` option from the
`File` menu. Then, select the file with extension `.p` that was output by `MrBayes` to load it 
onto `Tracer`. 

We will go through the most important features of `Tracer` together but, in general, 
we will focus on the effective sample sizes (ESSs) calculated 
for each of the model parameters, the frequency plot of the samples, and the trace plots. 


>> **QUESTIONS**:   
>>  - Do you think we need to run the chain longer?   
>>  - Is the ESS enough for the model parameters?   
>>  - How efficient is the chain?   
>>  - Do we need to specify a specific value for the burn-in phase?   
>>  - Which is the mean root age for this timetree?   
>>  - If another MCMC analysis has finished for another model you have run, compare 
>>    the results you have obtained and discuss with your group.


# 5. Viewing the annotated tree
`FigTree` is a user-friendly, graphical program for viewing trees. If you want to run `FigTree` from 
the command line, you can either set an alias in your `~/.bashrc` or `~/.bash_profile` 
file, run it from
the directory where the `jar` file can be found, or execute this `jar` file 
using an absolute or a relative path to the directory where the terminal is running.
Below, you can find a summary of how you would proceed for each of these three options 
when using the terminal:

```sh 
# Option 1: Run from the directory where 
# `figtree.jar` can be found. For instance, 
# if you open a terminal from
# `FigTree_v1.4.4/lib`, then:
java -jar figtree.jar


# Option 2: Use an absolute or a relative 
# path to execute the file. An example of 
# an absolute path is shown below:
java -jar Applications/FigTree_v1.4.4/lib/figtree.jar

# Option 3: Add an alias in your 
# `~/.bashrc` or `~/.bash_profile.
# Please replace "<path_to_FigTree>" with the absolute path
# to the location where you have saved `FigTree`
# Modify the name of the directory that you have unzipped
# if needed too (e.g., `FigTree_v1.4.4`)
alias figtree1.4.4='java -jar <path_to_Tracer>/FigTree_v1.4.4/lib/figtree.jar'
# Now, open a terminal from any location on 
# your PC and type the following command to 
# execute `FigTree`. Before you do this, please 
# make sure that you have X11 installed or, if 
# you are on Windows, that you have the Xming 
# Server running on your PC.
# If you are on Windows, you might need to run 
# the command `export DISPLAY=:0.0` before 
# you can execute `FigTree`
figtree1.4.4
```

Now, you can open in `FigTree` the file with the consensus tree that `MrBayes` has output,
i.e., file name that ends with `con.tre`. 
The tree will be displayed in the `FigTree` window. On the left hand side of the
window you can find the options and settings which control how the tree is displayed. We will 
see together the main options you can use to display the tree.

