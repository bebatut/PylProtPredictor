Welcome to PylProtPredictor's documentation!
============================================

Pyrrolysine is an amino acid that is used in the biosynthesis of proteins in some methanogenic archaea and bacterium. It is encoded in mRNA by the UAG codon, which in most organisms is the 'amber' stop codon.

Some methanogenic archaea and bacterium have the pylT gene, which encodes an unusual transfer RNA (tRNA) with a CUA anticodon, and the pylS gene, which encodes a class II aminoacyl-tRNA synthetase that charges the pylT-derived tRNA with pyrrolysine. In some proteins, the UAG codon can then code for pyrrolysine, and no more for a STOP codon.

These proteins are difficult to identify. Indeed, in CDS prediction, UAG codons are seen as STOP codons. The predicted CDS are then cut when the first UAG codon is found.

Here, we propose a solution to detect proteins using Pyrrolisine amino acid.

.. toctree::
  :caption: PylProtPredictor documentation
  :maxdepth: 4

  pylprotpredictor
  use_case
  installation
  usage
  code  
  contributing
  Source code on GitHub <https://github.com/bebatut/PylProtPredictor>