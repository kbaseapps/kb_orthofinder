# kb_orthofinder release notes
=========================================

2.1.2
-----
* OrthoFinder will replace a few special characters within the protein
  ids with the underscore character, so this fix will make sure that 
  the original protein ids can be retrieved.

2.1.1
-----
* Major bug fix to handle possible missing links between features,
  mrnas, and cdss in genome object.
* Expansion of testing configuration to make it easier to iterate on a
  single large genome.

2.0.1
-----
* Fix in table output where predicted orthologs were listed in multiple rows

2.0.0
-----
* Major release of kb_orthofinder

    * The reference genomes used to generate the protein families were
      updated to include the latest gene model versions released as
      part of Phytozome V12

    * "Annotate Plant Enzymes with OrthoFinder" has a new output table
      which captures the homologous relationships within the protein
      families and shows the user any putative orthologs that did not
      make the "cut"

    * Updated to use the Python 3 and minor bug fixes

1.0.0
-----
* Major release of kb_orthofinder
* Includes bugfix for missing CDS fields

0.0.3
-----
* Improved docs

0.0.2
-----
* Added citations in PLOS format

0.0.0
-----
* Module created by kb-sdk init
