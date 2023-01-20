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
from typing import  Union



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

chain_dict = {
    1:"A",
    2:"B",
    3:"C",
    4:"D",
    5:"E",
    6:"F",
    7:"G",
    8:"H",
    9:"I",
}

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
    
    model: str (default:None)
        if you want to use model thats not the best rated one
    
    use_relaxed: bool (default: True)
        wheather or not to use the amberrelaxed model
        
    """
    def __init__(self,path,model=None,use_relaxed=True):
        self.path = path
        if model is None:
            self._get_best_model()
        else:
            self._check_model(model=model)
        self._get_meta_data()
        self._get_pdb(use_relaxed=use_relaxed)
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
    
    def _check_model(self,model=None)->None:
        """
        checks for best alphafold model
        """
        with open(f'{self.path}/ranking_debug.json', 'r') as j:
            ranking_data = json.load(j)
        if model in ranking_data['order']:
            self.best_model = model
            if "multimer" in self.best_model:
                self.type = "multimer"
            else:
                self.type = "monomer"
        else:
            raise NameError (f"The model not in use please used one of the following models {*ranking_data['order'],}")
                
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
            self.obs = pd.DataFrame({"aa":[x for x in seq],
                                    "chain": ["A" for x in seq],
                                    "chain_postion":[x for x in range(len(seq))],
                                    "pLDDT": self.meta_data["plddt"]})
        else:
            seq = self.feature_data["aatype"]
            self.obs = pd.DataFrame({"aa":_translate_numbers_to_one_letter(seq),
                                     "chain": [chain_dict[x] for x in self.feature_data["asym_id"].astype(int)],
                                     "chain_postion":self.feature_data["residue_index"].astype(int),
                                     "pLDDT": self.meta_data["plddt"]})
            
    def _get_pdb(self,use_relaxed=True):
        """
        loads best relaxed `.pdb` file
        """
        if use_relaxed:
            with open(f"{self.path}/relaxed_{self.best_model}.pdb") as ifile:
                self.pdb = "".join([x for x in ifile])
        else:
            with open(f"{self.path}/unrelaxed_{self.best_model}.pdb") as ifile:
                self.pdb = "".join([x for x in ifile])

    def plot_pae(self,color="r", save = False):
        """
        Plots Predicted Aligned Error (PAE)
        
        Parameters
        ----------
        color:str (default:"r")
            color to show multimer protein boundaries. Has no effect on monomer proteins.
        save: bool (default: False)
            if you want to save PAE plot in alphafold_prediction_folder/figures
            
        """
        if "predicted_aligned_error" in self.meta_data.keys():
            plt.imshow(self.meta_data["predicted_aligned_error"], vmin=0., vmax=self.meta_data["max_predicted_aligned_error"], cmap='Greens_r')
            if self.type=="multimer":
                counter = 0
                sizes = self.obs.groupby("chain").size()
                for idx,length in enumerate(sizes):
                    if idx+1<len(sizes):
                        plt.hlines(length+counter-1,xmin=0,xmax= self.obs.shape[0],
                        color=color, linestyle='dashed')
                        plt.vlines(length+counter-1,ymin=0,ymax= self.obs.shape[0],
                        color= color, linestyle='dashed')
                        counter+=length
            plt.colorbar(fraction=0.046, pad=0.04)
            plt.ylim(0,self.obs.shape[0]-1)
            plt.xlim(0,self.obs.shape[0]-1)

            if save:
                self._create_figure_folder()
                plt.savefig(f'{self.path}/figures/pae_{self.best_model}.png',dpi = 600)
            plt.show()
            plt.close()
        else:
            print(f"Please check your used Alphafold Model. You might want to run `monomer_ptm`")

    def plot_pLDDT(self, add_pLDDT_background= False, save = False):
        """
        Plots per-residue confidence measure (pLDDT) for each aa
        
        Parameters
        ----------
        add_pLDDT_background: bool (default: False)
            adds a colored background based on pLDDT
        save: bool (default: False)
            if you want to save pLDDT plot in alphafold_prediction_folder/figures
            
        """
        y = self.meta_data["plddt"]
        x = np.arange(len(y))
        if add_pLDDT_background:
            plt.plot(x,y,c="black")
            plt.axhspan(0, 50, facecolor='#ff7d45')
            plt.axhspan(50, 70, facecolor='#ffdb13')
            plt.axhspan(70, 90, facecolor='#65cbf3')
            plt.axhspan(90, 100, facecolor='#0053d6')
        else:
            plt.plot(x,y)
        if self.type=="multimer":
            counter = 0
            sizes = self.obs.groupby("chain").size()
            for idx,length in enumerate(sizes):
                if idx+1<len(sizes):

                    plt.vlines(length+counter-1,ymin=0,ymax= 100,
                    color= "k", linestyle='dashed')
                    counter+=length
        plt.title("Predicted LDDT")
        plt.xlabel("Residue") 
        plt.ylabel("pLDDT")
        plt.ylim(0,100)
        plt.xlim(0,len(y))
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


    def show_annotation(self, style= "cartoon", annotation_key = "aa"):
        """
        Creates a py3Dmol view colored by a custom annotaion column in `.obs`
        In multimere model use `key = 'chain'` to color each protein chain separately.
        
        Parameters
        ----------
        style: str (default:'catroon')
            The style to display the py3Dmol view in
        
        annotation_key: str (defaut: aa)
            The column key in `.obs` to use for coloring.
            By default it colors based on aa type.
            
        """
        if annotation_key == "pLDDT":
            self.show_confidence(style=style)
        else:
            view = py3Dmol.view(width=800, height=800)
            view.addModelsAsFrames(self.pdb)
            if style not in styles:
                raise TypeError("style must be `cartoon`, `stick` or `sphere`")
            self.obs[annotation_key] = self.obs[annotation_key].astype("category")
            thresh = self.obs[annotation_key].cat.categories.to_list()
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
                cdx = colordict[self.obs.loc[idx,annotation_key]]
                view.setStyle({'model': -1, 'serial': i+1}, {style: {'color': cmap[cdx]}})
            
            view.zoomTo()
            grid = GridspecLayout(1, 2)
            out = Output()
            with out:
                view.show()
            grid[0, 0] = out

            out = Output()
            with out:
                _plot_plddt_legend(thresh,plddt_bands, title=annotation_key).show()
            grid[0, 1] = out
            display.display(grid)

    def show_substructure(self,
                        split_key: str,
                        use_cat: Union[str, list],
                        annotation_key: str= "pLDDT",
                        style = "cartoon"):
        """
        Creates a py3Dmol view colored by a custom annotaion column in `.obs`
        In multimere model use `key = 'chain'` to color each protein chain separately.
        
        Parameters
        ----------
        split_key: str
            The column key in `.obs` to choose the 
            substructure(s) from.
        
        use_cat: Union[str, list],
            String or List of substructure(s) to view.
            This can even be a single

        annotation_key:
            The column key in `.obs` to use for coloring.
            By default it colors based on pLDDT.

        style: str (default:'catroon')
            The style to display the py3Dmol view in
        """
        pdb_list = self.pdb.split("\n")
        subset_pdb = []
        if type(use_cat) is str:
            use_cat = [use_cat]
        subset_obs = self.obs.loc[self.obs[split_key].isin(use_cat)].copy()
        subset_obs["idx"] =subset_obs.chain+(subset_obs.chain_postion+1).astype(str)
        valid_idx= set(subset_obs["idx"])
        for line in pdb_list:
            split = line.split()
            if len(split) == 0: 
                pass
            elif split[0] == "TER":
                if split[3] in subset_obs.chain:
                    subset_pdb.append(line)
            elif split[0] == "ATOM":
                if split[4]+split[5] in valid_idx:
                    subset_pdb.append(line)
        
        idx = 1 
        indexed_subset_pdb = []
        for line in subset_pdb:
            new_line = line[:5]+f"              {str(idx)}"[-6:]+line[11:]
            idx+=1
            indexed_subset_pdb.append(new_line)
        indexed_subset_pdb = "\n".join(indexed_subset_pdb)
        
        subset_obs["adjusted"] = None
        view = py3Dmol.view(width=800, height=800)
        view.addModelsAsFrames(indexed_subset_pdb)
        if style not in styles:
            raise TypeError("style must be `cartoon`, `stick` or `sphere`")
        if annotation_key == "pLDDT":
            plddt_bands = ['#FF7D45','#FFDB13','#65CBF3','#0053D6']
            for i,line in enumerate(indexed_subset_pdb.split("\n")):
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
        else:
            subset_obs[annotation_key] = subset_obs[annotation_key].astype("category")
            thresh = subset_obs[annotation_key].cat.categories.to_list()
            colordict = {}
            for idx, i in enumerate(thresh):
                colordict[i] = idx
            plddt_bands = cmap[:len(thresh)]
            for i,line in enumerate(indexed_subset_pdb.split("\n")):
                split = line.split()
                if len(split) == 0: 
                    continue
                elif split[0] != "ATOM":
                    continue
                idx_c = split[4]+split[5]
                cdx = colordict[subset_obs.loc[subset_obs["idx"]==idx_c,annotation_key].values[0]]
                view.setStyle({'model': -1, 'serial': i+1}, {style: {'color': cmap[cdx]}})
            
            view.zoomTo()
            grid = GridspecLayout(1, 2)
            out = Output()
            with out:
                view.show()
            grid[0, 0] = out

            out = Output()
            with out:
                _plot_plddt_legend(thresh,plddt_bands, title=annotation_key).show()
            grid[0, 1] = out
            display.display(grid)
