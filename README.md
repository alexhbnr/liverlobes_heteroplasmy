# Workflow for extraction of mitochondrial DNA heteroplasmies from Illumina sequencing data

Scripts for extraction of mitochondrial DNA heteroplasmies from Illumina sequencing data and
detection of possible contamination as described in:

**Selection and random genetic drift shape the sharing of mtDNA heteroplasmies between human liver lobes**

Alexander Hübner, Manja Wachsmuth, Roland Schröder, Mingkun Li, Anna Maria
Eis-Hübinger, Burkhard Madea, Mark Stoneking

[bioRxiv DOI 10.1101/155796](https://doi.org/10.1101/155796)

The scripts uses as input data the output file produced by the [DREEP pipeline](http://dmcrop.sourceforge.net/index.html).

The Python scripts require the following Python modules:

* NumPy (>= v0.11)
* Pandas (>= v0.17)
