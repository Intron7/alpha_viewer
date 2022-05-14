# alpha_viewer
Class to view Alphafold models

alpha_viewer is based on the alphafold2 colab notebook visualization. It automatically chooses the best available prediction. alpha_viewer includes functions to plot PAE (`plot_pae`) and pLDDT (`plot_pLDDT`). It also allows for coloring of the py3Dmol view based on pLDDT with `show_confidence`. alpha_viewer contains `.obs`, a pandas data-frame, that can be used for custom annotations based on aa postion within the chain(s). You can use a key of `.obs` to color the py3Dmol view of your protein(s) based that annotation with `show_annotation`.

alpha-viewer has been tested to work with `monomer`, `monomer_ptm` and `multimer` models.
So far it's sadly not possible to save the py3Dmol view for exporting.

If you use alpha_viewer in your publication please cite: <a href="https://zenodo.org/badge/latestdoi/475046878"><img src="https://zenodo.org/badge/475046878.svg" alt="DOI"></a>



-------
### Install

you can simply install this repository from github with:\
`pip install git+https://github.com/Intron7/alpha_viewer.git`


It's recommended to use alpha_viewer within jupyterlab. Please enable jupyterlab extensions and install `jupyter-widgets/jupyterlab-manager` and `jupyterlab_3dmol`
