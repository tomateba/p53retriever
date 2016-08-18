# p53retriever

R package to locate and display p53 putative response elements on a DNA sequence
------------------------------------------------------------------------

### Overview

p53retriever is an R package implementing a pattern search algorithm that maps p53 response elements (REs) and ranks them according to predicted transactivation potentials in five classes. Besides canonical, full site REs, p53retriever contains specific pattern searches for non-canonical half sites and 3/4 sites. REs can be displayed along the input sequence in a graphical layout.

The rules formalized in the model are based on observations obtained during several years of experiments using the yeast-based assays, including the results presented in (Inga et al MCB2002; Tomso et al, PNAS 2005; Jegga et al PNAS 2008; Jordan et al PloS Genetics 2008; Menendez et PNAS 2010).

------------------------------------------------------------------------

### Reference

Please cite the following article when using `p53retriever`:

T Tebaldi, S Zaccara, F alessandrini, A Bisio, Y Ciribilli, A Inga. ***Whole-genome cartography of p53 response elements ranked on transactivation potential.*** *BMC Genomics. 2015 Jun 17;16:464*

[![doi](https://img.shields.io/badge/DOI-10.1186%2Fs12864--015--1643--9-green.svg?style=flat)](http://dx.doi.org/10.1186/s12864-015-1643-9)

------------------------------------------------------------------------

### License

[![Packagist](https://img.shields.io/packagist/l/doctrine/orm.svg?maxAge=2592000?style=flat)](https://opensource.org/licenses/MIT)

------------------------------------------------------------------------

### Output example

p53 REs identified in the promoter region of CDKN1A (-10k, +10k from TSS)

![p53 REs identified in the promoter of CDKN1A](https://cloud.githubusercontent.com/assets/9716233/6002270/c460ac18-aae9-11e4-8822-7f5272396634.png)

------------------------------------------------------------------------

### Contact
t.tebaldi@unitn.it
