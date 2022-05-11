# alpha_viewer
Class to view Alphafold models

alpha_viewer is based on the alphafold2 colab notebook visualization. alpha_viewer includes functions to plot PAE (`plot_pae`) and pLDDT (`plot_pLDDT`). It also allows for coloring of the py3Dmol view based on pLDDT with `show_confidence`. alpha_viewer contains `.obs`, a pandas data-frame, that can be used for custom annotations based on the aa postion within the chain(s). You can use a key of `.obs` to color the py3Dmol view of your protein(s) based that annotation with `show_annotation`.

alpha-viewer works has been test to work with `monomer`, `monomer_ptm` and `multimer` models.
So far it's sadly not possible to save the py3Dmol view into a plot for exporting.

If you use alpha_viewer in your publication please cite: 
