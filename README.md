
1. embedding_distance.py - this script calculates the euclidean distance between protein sequences. I wrote this script to calculate the distances for a set of proteins with unknown products to those of the known products. Only the binding site residues were used for calculating the binding site distances. Files containing the embeddings, the list of binding site residues are used as inputs. The distances were calculated between residues that were at the identical positions in the binding site.

2. interface_label.py - this script was written to label residues as interface or non-interface based on a distance cut-off. We already have a file that has the distances between the residues within a particular cut off distance. PDB files and the atom distance files are used as inputs.
