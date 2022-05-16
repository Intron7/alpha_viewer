import pickle
import py3Dmol
import matplotlib.pyplot as plt
import numpy as np
import json
import pandas as pd
import os
from IPython import display
from ipywidgets import GridspecLayout
from ipywidgets import Output




aa_dict = {0:["Alanine","Ala","A"],
           1:["Arginine","Arg","R"],
           2:["Asparagine","Asn","N"],
           3:["Aspartic acid","Asp","D"],
           4:["Cysteine","Cys","C"],
           5:["Glutamic acid","Glu","E"],
           6:["Glutamine","Gln","Q"],
           7:["Glycine","Gly","G"],
           8:["Histidine","His","H"],
           9:["Isoleucine","Ile","I"],
           10:["Leucine","Leu","L"],
           11:["Lysine","Lys","K"],
           12:["Methionine","Met","M"],
           13:["Phenylalanine","Phe","F"],
           14:["Proline","Pro","P"],
           15:["Serine","Ser","S"],
           16:["Threonine","Thr","T"],
           17:["Tryptophan","Trp","W"],
           18:["Tyrosine","Tyr","Y"],
           19:["Valine","Val","V"]}

cmap = ['#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#e377c2', '#b5bd61', '#17becf', '#aec7e8',
        '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31']

styles = ["cartoon","stick","sphere"]


def _translate_numbers_to_one_letter(sequence):
    """
    Translates numeric aa code to singleletter aa 
    """
    return [aa_dict[x][2] for x in sequence]
              
def _plot_plddt_legend(thresh, colors, title = None):
    """
    Creates Anotation for py3DMOL plots
    """
    plt.figure(figsize=(2, 2))
    for c in colors:
        plt.bar(0, 0, color=c)
    plt.legend(thresh, frameon=False, loc='center', fontsize=20)
    plt.xticks([])
    plt.yticks([])
    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    if title:
        plt.title(title, fontsize=20, pad=20)
    return plt

class alpha_viewer:
    """
    alpha_viewer creates 3D views of Alphafold Predictions. It's based on the alphafold colab notebooks visualistations.
    It automatically chooses the best model of `ranking_debug.json`.
    Contains `.obs` a pandas dataframe for easy custom annotations and visualisation.
        
    Parameters
    ----------
    path
        path to alphafold output folder
        
    """
    def __init__(self,path):
        self.path = path
        self._get_best_model()
        self._get_meta_data()
        self._get_pdb()
        self._get_feature_data()
        self._figures = False

    def _create_figure_folder(self)-> None:
        """
        checks for existing figure folder within alphafold folder
        """
        if self._figures == False:
            os.makedirs(f"{self.path}/figures", exist_ok=True)
            self._figures = True

    def _get_best_model(self)->None:
        """
        checks for best alphafold model
        """
        with open(f'{self.path}/ranking_debug.json', 'r') as j:
            ranking_data = json.load(j)
        self.best_model = ranking_data['order'][0]
        if "multimer" in self.best_model:
              self.type = "multimer"
        else:
              self.type = "monomer"
                
    def _get_meta_data(self):
        """
        loads `features.pkl` file
        """
        with open(f'{self.path}/result_{self.best_model}.pkl', 'rb') as f:
            self.meta_data = pickle.load(f)
            print(f"using model: {self.best_model} with pLDDT {self.meta_data['ranking_confidence']}")

    def _get_feature_data(self):
        """
        loads best `result_model_*.pkl` file
        """
        with open(f'{self.path}/features.pkl', 'rb') as f:
            self.feature_data = pickle.load(f)
        if self.type == "monomer":    
            seq = self.feature_data["sequence"][0].decode("utf8")
            self.obs = pd.DataFrame({"aa":[x for x in seq]})
        else:
            seq = self.feature_data["aatype"]
            self.obs = pd.DataFrame({"aa":_translate_numbers_to_one_letter(seq),
                                     "chain": self.feature_data["asym_id"].astype(int),
                                     "chain_postion":self.feature_data["residue_index"].astype(int)})
            
    def _get_pdb(self):
        """
        loads best relaxed `.pdb` file
        """
        with open(f"{self.path}/relaxed_{self.best_model}.pdb") as ifile:
            self.pdb = "".join([x for x in ifile])

    def plot_pae(self, save = False):
        """
        Plots Predicted Aligned Error (PAE)
        
        Parameters
        ----------
        save: bool (default: False)
            if you want to save PAE plot in alphafold_prediction_folder/figures
            
        """
        if "predicted_aligned_error" in self.meta_data.keys():
            plt.imshow(self.meta_data["predicted_aligned_error"], vmin=0., vmax=self.meta_data["max_predicted_aligned_error"], cmap='Greens_r')
            plt.colorbar(fraction=0.046, pad=0.04)
            if save:
                self._create_figure_folder()
                plt.savefig(f'{self.path}/figures/pae_{self.best_model}.png',dpi = 600)
            plt.show()
            plt.close()
        else:
            print(f"Please check your used Alphafold Model. You might want to run `monomer_ptm`")

    def plot_pLDDT(self, save = False):
        """
        Plots per-residue confidence measure (pLDDT) for each aa
        
        Parameters
        ----------
        save: bool (default: False)
            if you want to save pLDDT plot in alphafold_prediction_folder/figures
            
        """
        y = self.meta_data["plddt"]
        x = np.arange(len(y))
        plt.plot(x, y)
        plt.title("Predicted LDDT")
        plt.xlabel("Residue")
        plt.ylabel("pLDDT")
        plt.ylim(0,100)
        if save:
            self._create_figure_folder()
            plt.savefig(f'{self.path}/figures/pLDDT_{self.best_model}.png',dpi = 600)
        plt.show()
        plt.close()

    def show_confidence(self, style= "cartoon", show_sidechains = False):
        """
        Creates a py3Dmol view in the style of the alphafold colab-notebook colored by pLDDT
        
        Parameters
        ----------
        style: str (default:'cartoon')
            The style to display the py3Dmol view in
        
        show_sidechains: bool (default: False)
            Weather or not to use stick style in addition to cartoon style.
            Will be ignored if style is not `cartoon`
            
        """
        view = py3Dmol.view(width=800, height=800)
        view.addModelsAsFrames(self.pdb)
        plddt_bands = ['#FF7D45','#FFDB13','#65CBF3','#0053D6']
        if style not in styles:
            raise TypeError("style must be `cartoon`, `stick` or `sphere`")
        for i,line in enumerate(self.pdb.split("\n")):
            split = line.split()
            if len(split) == 0 or split[0] != "ATOM":
                continue
            if float(split[-2]) >90:
                color = plddt_bands[3]
            elif float(split[-2]) >70:
                color = plddt_bands[2]
            elif float(split[-2]) >50:
                color = plddt_bands[1]
            else:
                color = plddt_bands[0]
            
            if style == "cartoon" and show_sidechains == True:
                view.setStyle({'model': -1, 'serial': i+1}, {style: {'color': color},"stick":{}})
            else:
                view.setStyle({'model': -1, 'serial': i+1}, {style: {'color': color}})
        view.zoomTo()
        grid = GridspecLayout(1, 2)
        out = Output()
        with out:
            view.show()
        grid[0, 0] = out

        out = Output()
        with out:
            thresh = ['Very low (pLDDT < 50)',
                'Low (70 > pLDDT > 50)',
                'Confident (90 > pLDDT > 70)',
                'Very high (pLDDT > 90)']
            _plot_plddt_legend(thresh,plddt_bands,'Model Confidence').show()
        grid[0, 1] = out
        display.display(grid)


    def show_annotation(self, style= "cartoon", key = "aa"):
        """
        Creates a py3Dmol view colored by a custom annotaion column in `.obs`
        In multimere model use `key = 'chain'` to color each protein chain separately.
        
        Parameters
        ----------
        style: str (default:'catroon')
            The style to display the py3Dmol view in
        
        key: str (defaut: aa)
            The column key in `.obs` to use for coloring.
            By default it colors based on aa type.
            
        """
        view = py3Dmol.view(width=800, height=800)
        view.addModelsAsFrames(self.pdb)
        if style not in styles:
            raise TypeError("style must be `cartoon`, `stick` or `sphere`")
        self.obs[key] = self.obs[key].astype("category")
        thresh = self.obs[key].cat.categories.to_list()
        colordict = {}
        for idx, i in enumerate(thresh):
            colordict[i] = idx
        plddt_bands = cmap[:len(thresh)]
        multi_offset = 0 
        for i,line in enumerate(self.pdb.split("\n")):
            split = line.split()
            if len(split) == 0: 
                continue
            elif split[0] != "ATOM":
                if split[0] == "TER":
                    multi_offset += int(split[4])
                continue
            idx = int(split[5])-1 + multi_offset
            cdx = colordict[self.obs.loc[idx,key]]
            view.setStyle({'model': -1, 'serial': i+1}, {style: {'color': cmap[cdx]}})
        
        view.zoomTo()
        grid = GridspecLayout(1, 2)
        out = Output()
        with out:
            view.show()
        grid[0, 0] = out

        out = Output()
        with out:
            _plot_plddt_legend(thresh,plddt_bands, title=key).show()
        grid[0, 1] = out
        display.display(grid)


   