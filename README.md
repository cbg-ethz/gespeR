gespeR: Gene-Specific Phenotype EstimatoR
======

gespeR is a novel model to estimate gene-specific phenotypes from off-target confounded RNAi screens. The observed phenotype for each siRNA is modeled as a the weighted linear combination of gene-specific phenotypes from the on- and all off-target genes. This deconvolution approach yields highly reproducible phenotypes, essential for unbiased analyses of siRNA screening data.

###### Reference: 
Fabian Schmich, Ewa Szczurek, Saskia Kreibich, Sabrina Dilling, Daniel Andritschke, Alain Casanova, Shyan Huey Low, Simone Eicher, Simone Muntwiler, Mario Emmenlauer, Pauli Ramo, Raquel Conde-Alvarez, Christian von Mering, Wolf-Dietrich Hardt, Christoph Dehio and Niko Beerenwinkel.  
<b>[gespeR: a statistical model for deconvoluting off-target-confounded RNA interference screens](http://www.genomebiology.com)</b>  
<i>Genome Biology</i>, 2015.
 
### Installation
gespeR is hosted on GitHub and available through [Bioconductor](http://bioconductor.org/packages/gespeR/). The package is released under the GNU General Public License (GPL) version 3 and includes examples and a [vignette](http://bioconductor.org/packages/devel/bioc/vignettes/gespeR/inst/doc/gespeR.pdf).

### Data
In addition to the phenotypic readout, gespeR requires siRNA-to-gene target relation matrices, quantifying how strongly each siRNA downregulates transcript genes via on- and off-targeting. These matrices can be computed with miRNA target prediction tools, such as for instance [TargetScan](http://www.targetscan.org). 

###### Target Relation Matrices
Below, we provide pre-computed matrices for all libraries used in Schmich et al., 2015. Wrapper scripts to run TargetScan in batch mode for the prediction of siRNA-to-gene target relation matrices are also available on [GitHub](http://github.com/fschmich/tscan_wrapper), including a README with step-by-step instructions. All pre-computed siRNA-to-gene target relation matrices are stored in .rds files using R's serialization interface for single objects. Load the data into R by using the function readRDS(). Note that loading target relation matrices can require up to 5GB of RAM.
- [Ambion](http://n.ethz.ch/~fschmich/gespeR/IFX_AMBION.rds)
- [Dharmacon](http://n.ethz.ch/~fschmich/gespeR/IFX_DHARMACON.rds)
- [Dharmacon deconvoluted](http://n.ethz.ch/~fschmich/gespeR/IFX_DHARMACON_DECON.rds)
- [Qiagen](http://n.ethz.ch/~fschmich/gespeR/IFX_QIAGEN.rds)
- [Qiagen kinome](http://n.ethz.ch/~fschmich/gespeR/IFX_QIAGEN_KINOME.rds)
- [Validation](http://n.ethz.ch/~fschmich/gespeR/IFX_VALIDATION.rds)
- [Schultz et al., 2011](http://n.ethz.ch/~fschmich/gespeR/SCHULTZ.rds)

###### Pathogen Infection Screen Phenotypes
High-content, image based phenotypes for pathogen infection RNAi screens from the [InfectX consortium](http://www.infectx.ch) are hosted at [PubChem](https://pubchem.ncbi.nlm.nih.gov/assay/assay.cgi?aid=1117357).

### Usage
Step-by-step instructions demonstrating how to download, pre-process and deconvolute pathogen infection screen phnotypes is available in form of an [R/Vignette](http://edit.ethz.ch/bsse/cbg/software/gespeR/Deconvolute_Pathogen_Screens).

### Contributions
- [Fabian Schmich](http://www.bsse.ethz.ch/cbg/group/people/person-detail.html?persid=183865)
- [Ewa Szczurek](mailto:ewa.szczurek@bsse.ethz.ch)
- [Niko Beerenwinkel](http://www.bsse.ethz.ch/cbg/group/people/person-detail.html?persid=149417)
 
###Contact
```
Fabian Schmich
fabian.schmich (at) bsse.ethz.ch
```
