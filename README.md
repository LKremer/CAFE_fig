CAFE_fig, a visualization tool for CAFE.
=========

CAFE_fig takes a .cafe output file and produces:
- a summary tree that shows the average expansion/contraction of families across the phylogeny
- a tree that denotes which branches evolve under which lambda (if available)
- a tree for each family of interest, i.e. families that the user specified by ID or families that showed significant change at a user-specified clade of interest



Requirements
------------


CAFE_fig requires Python3.4+ and ETE3:

`pip3 install 'ete3==3.0.0b35'`


Usage
------------

```
usage: CAFE_fig.py [-h] [-f FAMILIES [FAMILIES ...]] [-c CLADES [CLADES ...]]
                   [-a ALPHA_ERROR] [-d DUMP] [-g GFX_OUTPUT_FORMAT]
                   report_cafe

Parses a CAFE output file (.cafe) and plots a summary tree that shows the
average expansion/contractin across the phylogeny; a tree that shows which
clades evolved under the same lambda (if available); and a gene family
evolution tree for each user-specified gene family.

positional arguments:
  report_cafe           the file report.cafe (or similar name)

optional arguments:
  -h, --help            show this help message and exit
  -f FAMILIES [FAMILIES ...], --families FAMILIES [FAMILIES ...]
                        only show families with these IDs
  -c CLADES [CLADES ...], --clades CLADES [CLADES ...]
                        only show families that are expanded/contracted at
                        this clade. Format: [clade]=[leaf],[leaf] where clade
                        is the name of the last common ancestor of the two
                        leaves, e.g.: Isoptera=zne,mna
  -a ALPHA_ERROR, --alpha_error ALPHA_ERROR
                        p-value cutoff (default: 0.05)
  -d DUMP, --dump DUMP  don't open trees in a window, write them to files in
                        the specified directory instead (default: off)
  -g GFX_OUTPUT_FORMAT, --gfx_output_format GFX_OUTPUT_FORMAT
                        output format for the tree figures when using --dump
                        [svg|pdf|png] (default: pdf)
```

Example usage
------------

`./CAFE_fig.py example_result.cafe -c Isoptera=zne,mna -a 0.05 --dump test/ -g .pdf`

Reads "example_result.cafe" and dumps all figures in PDF format to the directory "test/". Trees will only be shown for families that showed a significant (p<=0.05) expansion/contraction at the node "Isoptera", which is the last common ancestor of "zne" and "mna".

Significant contractions are marked in magenta, significant expansions are marked in green (p<=0.001 = ******, p<=0.01 = ****, p<=0.05 = *).

![example_tree](http://i.imgur.com/221ra0l.png)
