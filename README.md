# BreakID-XNF

A modification of Jo Devriendt's [BreakID code](https://bitbucket.org/krr/breakid/src/master/) to handle the XNF (XOR-CNF) format encountered in cryptographic instances.

## Compiling

Simply run "make" in the "src/" folder. The executable "BreakID" will be compiled into that folder. To clean the build, run "make clean" in "src/".

To indicate that the instance is in XNF format, the option "-xnf" has to be specified inside the command line.

## Folders

This repository is composed of the following folders:

- src: contains the source code of our tool

- instances: contains the original and augmented instances considered during our experiments

- experiments: contains the results of our experiments. The folder "BreakID" contains the results returned by our tool when detecting symmetries. The folders "CMS-Base" and "CMS-Sym" contain the results obtained with the solver [CryptoMiniSat](https://github.com/msoos/cryptominisat) (version 5.11.19) respectively on the original and the augmented instances

# Reference

A paper is to appear in ICTAI 2026:

```
@inproceedings{BlommeCherifSaisICTAI2026,
  author       = {Anthony Blomme and
                  Sami Cherif and
                  Lakhdar Sais},
  title        = {Symmetry Detection and Breaking in Cryptographic XNF Instances},
  booktitle    = {38th {IEEE} International Conference on Tools with Artificial Intelligence,
                  {ICTAI} 2026, Boca Raton FL, USA, November 2-4, 2026},
  publisher    = {{IEEE}},
  year         = {2026},
}
```
