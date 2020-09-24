# Ogrooter
Ogrooter calculates the posterior probability distribution of the root position under the outgroup criteria as described in Huelsenbeck <i>et al</i>. (2002). We used an outgroup rooting method where each branch length has an independent exponential prior. To calculate the posterior probability distribution of the root under the molecular clock see <a href="https://github.com/NielsenBerkeleyLab/rooter">Rooter</a>.

Huelsenbeck, J. P., Bollback, J. P., and Levine, A. M. 2002. Inferring  the  root  of  a  phylogenetic  tree. Systematic biology, 51(1): 32–43.
# Prerequisites
To compile this program, you will need to download the Eigen library, which you can find at http://eigen.tuxfamily.org/. Make certain to change the appropriate line in the Makefile (modify the INC definition), pointing your compilation to the Eigen install on your computer.
# Installation

```S
git clone https://github.com/NielsenBerkeleyLab/ogrooter.git
cd ogrooter
make
make install
```
# Command-line usage
Ogrooter can be used from the command line with arguments specifying input data and parameters. Only PHYLIP alignment format is supported. Only newick format is supported for tree files. To run Ogrooter with multiple outgroups just use multiple arguments of <code>-og <name_of_outgroup1> -og <name_of_outgroup2></code>, etc. Outgroups are forced to be monophyletic. To run Ogrooter use:

```S
ogrooter -i <input.phylip> -t <input.nwk> -og <name_of_outgroup> -l <number of MCMC iterations> -s <MCMC sampling frequency> -p <MCMC printing frequency> -o <output_prefix>
```
# Outputs
Ogrooter prints a log file <code><output_prefix>.log</code>with a summary of each root position of the credible set along with its posterior probability and cumulative probability of the root branch. It also prints out the trees <code><output_prefix>.t</code>for each bipartition in Nexus format.
