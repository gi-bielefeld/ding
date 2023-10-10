# dingiiofficial_wrapper

This is just a wrapper on [dingII] to make it easier to `pip`- and `conda`- install. All rights belong to the [dingII] developers.

## Installation
```bash
pip install git+https://github.com/leoisl/dingiiofficial_wrapper
```

## Usage
```
$ dingII -h

usage: dingII [-h] [-V] {generate,parsesol} ...

ding II - an algorithm solving the genomic distance problem for natural genomes, in which any marker may occur an arbitrary number of times

positional arguments:
  {generate,parsesol}  Available subcommands
    generate           Generate the DING ILP or create a custom model file in order to fine tune how many genes per family are to be matched
    parsesol           Parse a gurobi solution into a distance and optionally give a matching

options:
  -h, --help           show this help message and exit
  -V, --version        show program's version number and exit
```

### Generate subcommand
```
$ dingII generate -h
usage: dingII generate [-h] [-mm | -em | -im | -r RANGE RANGE] [-c CUSTOM] [-p PAIR PAIR | -pn PAIRNUMBER PAIRNUMBER] (--writemodel WRITEMODEL | --writeilp WRITEILP) unimog

positional arguments:
  unimog                The genomes provided in UniMOG format.

options:
  -h, --help            show this help message and exit
  -mm, --maximal        Set matching model to maximal matching.
  -em, --exemplary      Set matching model to exemplary matching.
  -im, --intermediate   Set matching model to intermediate matching.
  -r RANGE RANGE, --range RANGE RANGE
                        Provide upper and lower percentiles to be matched per marker in range [0,1]. Actual discrete bounds will always be rounded up.
  -c CUSTOM, --custom CUSTOM
                        Provide a custom matching file.
  -p PAIR PAIR, --pair PAIR PAIR
                        Give the two names of the genomes you want to compare, as specified in the unimog header.
  -pn PAIRNUMBER PAIRNUMBER, --pairnumber PAIRNUMBER PAIRNUMBER
                        Chose the two genomes via their position in the file (starting at 0). Default: 0,1
  --writemodel WRITEMODEL
                        Write the matching model to a file in order to customize it.
  --writeilp WRITEILP   Write the resulting ILP to the specified file.
```

### parsesol subcommand
```
$ dingII parsesol -h
usage: dingII parsesol [-h] [-p PAIR PAIR | -pn PAIRNUMBER PAIRNUMBER] [-m MATCHING] --solgur SOLGUR [--runs RUNS] [--numindels] unimog

positional arguments:
  unimog                The genomes provided in UniMOG format.

options:
  -h, --help            show this help message and exit
  -p PAIR PAIR, --pair PAIR PAIR
                        Give the two names of the genomes you want to compare, as specified in the unimog header.
  -pn PAIRNUMBER PAIRNUMBER, --pairnumber PAIRNUMBER PAIRNUMBER
                        Chose the two genomes via their position in the file (starting at 0). Default: 0,1
  -m MATCHING, --matching MATCHING
                        Give the matching as a pair of indexed genomes.
  --solgur SOLGUR       Gurobi solution file with a single solution.
  --runs RUNS           Write runs of indels to the specified file. Format: Each line represents a cycle, each consecutive sequence of oriented markers is a run. Runs in the same cycle are separated by tab characters an begin
                        with an A-run or a tab character if no A-run exists.
  --numindels           Give a possible number of indels in the sorting scenario. Note that this number is NOT the same for all optimal scenarios.
```


Original [dingII] README follows.

# dingIIofficial

The python 3 Version of DING. This is the version where new features will still be added and refined. As of now, this version is still to be regarded as experimental as it has not been tested thoroughly. Feel free to report any bugs or ideas for new features!


An overview of the scripts dealing directly with the DING ILP:


|script  | purpose | input | output | dependencies |
| ------ | ------ | ------ | ------ | ------ |
|  dingII.py | Generate the DING ILP or create a custom model file in order to fine tune how many genes per family are to be matched | [UniMoG file](https://bibiserv.cebitec.uni-bielefeld.de/dcj?id=dcj_manual) containing a single genome pair, (custom model file)  | [gurobi lp file](https://www.gurobi.com/documentation/9.1/refman/lp_format.html) or custom model file | python3, networkx, dingII\_util, ilp\_util\_adj |
| dingII\_parsesol.py | Calculate the distance and matching | UniMoG file with the original, unmatched genomes, [gurobi solution file](https://www.gurobi.com/documentation/9.1/refman/sol_format.html) containing a single solution with objective value | UniMoG file with relabeled genes according to matching | python3, networkx, dingII\_util, ilp\_util\_adj|
| dingII\_util.py | Provide utility functions for MRD generation and solution parsing | - | - | python3,networkx, ilp\_util\_adj |
ilp\_util\_adj.py | Provide further utility functions | - | - | python3 |

A Typical workflow using the maximal matching model would look like this:
1.  Generate the ILP: `./dingII.py  {unimog-file} -mm --writeilp {ilp-file}`
2.  Use a solver to obtain a gurobi solution `{gurobi-sol}`.
3.  Get the matching and distance (and number of indels as well as a summary of runs ("indel-blocks")): `./dingII_parsesol.py {unimog-file} --solgur {gurobi-sol} --matching {unimog-matching} ` (`--numindels --runs {run-file}`)

where ` {unimog-file}` is the original unmatched genome pair in UniMoG-Format and `{gurobi-sol}` the gurobi solution file with objective value (e.g. `# Objective value = 12`).


<details><summary>Run file format</summary>

Run output file format (`{run-file}`):
- Each line encompasses all runs within the same cycle of the decomposition
- Runs within a cycle are separated by TAB-characters
- If there is no A-run the line begins with a TAB
- Runs are the concatenated string of the oriented markers to be deleted

------------------------------------
More formally:

`{cycle-1}`


`{cycle-2}`


`...`


with `{cycle-n}` = `{A-run}\tab{B-run}...` or `cycle-n` = `\tab {B-run}`

and `{X-run}`=`(+/-)indel1(+/-)indel2...` 

------------------------------------



</details>

The work of this repository is described in the following paper:
* Bohnenkämper L., Braga M.D.V., Doerr D., Stoye J.: Computing the Rearrangement Distance of Natural Genomes. In: Schwartz R. (eds) Proc. of RECOMB 2020. LNCS 12074, 3-18. Springer Verlag, 2020. [DOI](https://doi.org/10.1007/978-3-030-45257-5_1)

For a more detailed description see the following arxiv preprint:
*  Bohnenkämper, L., Braga, M.D.V., Doerr, D., Stoye, J.: Computing the rearrangement distance of natural genomes. [arXiv:2001.02139](http://arxiv.org/abs/2001.02139) (2020)


[dingII]: https://gitlab.ub.uni-bielefeld.de/gi/dingiiofficial