# p53retriever
R package to locate and display p53 putative response elements on a DNA sequence

### Overview

p53retriever is an R package implementing a pattern search algorithm that maps p53 response elements (REs) and ranks them according to predicted transactivation potentials in five classes. Besides canonical, full site REs, p53tracker contains specific pattern searches for non-canonical half sites and 3/4 sites. REs can be displayed along the input sequence in a graphical layout.

The rules formalized in the model are based on observations obtained during several years of experiments using the yeast-based assays, including the results presented in (Inga et al MCB2002; Tomso et al, PNAS 2005; Jegga et al PNAS 2008; Jordan et al PloS Genetics 2008; Menendez et PNAS 2010).

### Authors

Toma Tebaldi, Sara Zaccara, Federica Alessandrini, Alessandra Bisio, Yari Ciribilli and Alberto Inga.

### Output example

p53 REs identified in the promoter region of CDKN1A (-10k, +10k from TSS)

![p53 REs identified in the promoter of CDKN1A](https://cloud.githubusercontent.com/assets/9716233/6002270/c460ac18-aae9-11e4-8822-7f5272396634.png)

### Contact
t.tebaldi@unitn.it
