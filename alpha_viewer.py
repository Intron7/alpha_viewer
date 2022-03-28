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


def plot_plddt_legend(thresh, colors, title = None):
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
    def __init__(self,folder) -> None:
        self.path = folder
        self.get_best_model()
        self.get_meta_data()
        self.get_pdb()
        self.get_feature_data()
        self.figures = False
        self.cmap = ['#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#e377c2', '#b5bd61', '#17becf', '#aec7e8',
        '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31']

    def create_figure_folder(self)-> None:
        if self.figures == False:
            os.makedirs(f"{self.path}/figures", exist_ok=True)
            self.figures = True

    def get_best_model(self)->None:
        with open(f'{self.path}/ranking_debug.json', 'r') as j:
            ranking_data = json.load(j)
        self.best_model = ranking_data['order'][0]
    
    def get_meta_data(self):
        with open(f'{self.path}/result_{self.best_model}.pkl', 'rb') as f:
            self.meta_data = pickle.load(f)
            print(f"using model: {self.best_model} with pLDDT {self.meta_data['ranking_confidence']}")

    def get_feature_data(self):
        with open(f'{self.path}/features.pkl', 'rb') as f:
            self.feature_data = pickle.load(f)
        seq = self.feature_data["sequence"][0].decode("utf8")
        self.obs = pd.DataFrame({"aa":[x for x in seq]})

    def get_pdb(self):
        with open(f"{self.path}/relaxed_{self.best_model}.pdb") as ifile:
            self.pdb = "".join([x for x in ifile])


    def plot_pae(self, save = False):
        if "predicted_aligned_error" in self.meta_data.keys():
            plt.imshow(self.meta_data["predicted_aligned_error"], vmin=0., vmax=self.meta_data["max_predicted_aligned_error"], cmap='Greens_r')
            plt.colorbar(fraction=0.046, pad=0.04)
            if save:
                self.create_figure_folder()
                plt.savefig(f'{self.path}/figures/pae_{self.best_model}.png',dpi = 600)
            plt.show()
            plt.close()
        else:
            print(f"Please check your used Alphafold Model. You might want to run `monomer_ptm`")

    def plot_pLDDT(self, save = False):
        y = self.meta_data["plddt"]
        x = np.arange(len(y))
        plt.plot(x, y)
        plt.title("Predicted LDDT")
        plt.xlabel("Residue")
        plt.ylabel("pLDDT")
        plt.ylim(0,100)
        if save:
            self.create_figure_folder()
            plt.savefig(f'{self.path}/figures/pLDDT_{self.best_model}.png',dpi = 600)
        plt.show()
        plt.close()

    def show_confidence(self, style= "cartoon", show_sidechains = False):
        view = py3Dmol.view(width=800, height=800)
        view.addModelsAsFrames(self.pdb)
        plddt_bands = ['#FF7D45','#FFDB13','#65CBF3','#0053D6']
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
            idx = int(split[1])
            
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
            plot_plddt_legend(thresh,plddt_bands,'Model Confidence').show()
        grid[0, 1] = out
        display.display(grid)
        self.view = view

    def show_annotation(self, style= "cartoon", key = "aa"):
        view = py3Dmol.view(width=800, height=800)
        view.addModelsAsFrames(self.pdb)
        self.obs[key] = self.obs[key].astype("category")
        thresh = self.obs[key].cat.categories.to_list()
        colordict = {}
        for idx, i in enumerate(thresh):
            colordict[i] = idx
        plddt_bands = self.cmap[:len(thresh)]
        for i,line in enumerate(self.pdb.split("\n")):
            split = line.split()
            if len(split) == 0 or split[0] != "ATOM":
                continue
            idx = int(split[5])-1
            cdx = colordict[self.obs.loc[idx,key]]
            view.setStyle({'model': -1, 'serial': i+1}, {style: {'color': self.cmap[cdx]}})

        view.zoomTo()
        grid = GridspecLayout(1, 2)
        out = Output()
        with out:
            view.show()
        grid[0, 0] = out

        out = Output()
        with out:
            plot_plddt_legend(thresh,plddt_bands, title=key).show()
        grid[0, 1] = out
        display.display(grid)


    def show_glycosylation(self, style= "cartoon"):
        try:
            with open(f"{self.path}/glyc.txt") as glyc_raw:
                glyc_raw = glyc_raw.read().splitlines()
                glyc_raw = [x.split()[0] for x in glyc_raw]
                glyc_raw = "".join(glyc_raw)
                glyc_raw.replace(" ","")
                glyc_sites = glyc_raw
            view = py3Dmol.view(width=800, height=800)
            view.addModelsAsFrames(self.pdb)
            plddt_bands = ['red','white']
            for i,line in enumerate(self.pdb.split("\n")):
                split = line.split()
                if len(split) == 0 or split[0] != "ATOM":
                    continue
                idx = int(split[5])
                if glyc_sites[idx-1] == "N":
                    color = plddt_bands[0]
                else:
                    color = plddt_bands[1]
                view.setStyle({'model': -1, 'serial': i+1}, {style: {'color': color}})

            view.zoomTo()
            grid = GridspecLayout(1, 2)
            out = Output()
            with out:
                view.show()
            grid[0, 0] = out

            out = Output()
            with out:
                thresh = ['Possible Glycosylation',
                    'No Glycosylation']
                plot_plddt_legend(thresh,plddt_bands).show()
            grid[0, 1] = out
            display.display(grid)
        except IOError:
            print("Error: File does not appear to exist. Please create a `glyc.txt` in the alphafold output folder.")


