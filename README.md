# dingIIofficial

The python 3 Version of DING. This is the version where new features will still be added and refined. As of now, this version is still to be regarded as experimental as it has not been tested thoroughly. Feel free to report any bugs or ideas for new features!


An overview of the scripts dealing directly with the DING ILP:


|script  | purpose | input | output | dependencies |
| ------ | ------ | ------ | ------ | ------ |
|  dingII.py | Generate the DING ILP or create a custom model file in order to fine tune how many genes per family are to be matched | [UniMoG file](https://bibiserv.cebitec.uni-bielefeld.de/dcj?id=dcj_manual) containing a single genome pair, (custom model file)  | [gurobi lp file](https://www.gurobi.com/documentation/9.1/refman/lp_format.html) or custom model file | python3, networkx, dingII\_util |
| dingII\_parsesol.py | Calculate the distance and matching | UniMoG file with the original, unmatched genomes, [gurobi solution file](https://www.gurobi.com/documentation/9.1/refman/sol_format.html) containing a single solution with objective value | UniMoG file with relabeled genes according to matching | python3, networkx, dingII\_util|
| dingII\_util.py | Provide utility functions for generation and interpretation| - | - | python3,networkx, ilp\_util\_adj |

A Typical workflow using the maximal matching model would look like this:
1.  Generate the ILP: `./dingII.py  {unimog-file} -mm --writeilp {ilp-file}`
2.  Apply gurobi
3.  Get the matching and distance (and number of indels as well as a summary of runs ("indel-blocks")): `%s/dingII_parsesol.py {unimog-file} --solgur {gurobi-sol} --matching {unimog-matching} (--numindels --runs {run-file})`

where ` {unimog-file}` is the original unmatched genome pair in UniMoG-Format and `{gurobi-sol}` the gurobi solution file.


<details><summary>Run file format</summary>

Run output file format (`{run-file}`):
- Each line encompasses all runs within the same cycle
- Runs within a cycle are separated by TAB-characters
- If there is no A-run the line begins with a TAB
- Runs are the concatenated string of the oriented markers to be deleted

------------------------------------
More formally:

`{cycle-1}`


`{cycle-2}`


`...`


with `cycle-n = {A-run}\tab{B-run}...` or `cycle-n = \tab {B-run}`

with `{X-run}=(+/-)indel1(+/-)indel2...` 

------------------------------------



</details>


The work of this repository is described in the following paper:
* Bohnenkämper L., Braga M.D.V., Doerr D., Stoye J.: Computing the Rearrangement Distance of Natural Genomes. In: Schwartz R. (eds) Proc. of RECOMB 2020. LNCS 12074, 3-18. Springer Verlag, 2020. [DOI](https://doi.org/10.1007/978-3-030-45257-5_1)

For a more detailed description see the following arxiv preprint:
*  Bohnenkämper, L., Braga, M.D.V., Doerr, D., Stoye, J.: Computing the rearrangement distance of natural genomes. [arXiv:2001.02139](http://arxiv.org/abs/2001.02139) (2020)
