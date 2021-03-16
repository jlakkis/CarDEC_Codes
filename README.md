# CarDEC Evaluations

CarDEC (**C**ount **a**dapted **r**egularized **D**eep **E**mbedded **C**lustering) is a joint deep learning computational tool that is useful for analyses of single-cell RNA-seq data. The CarDEC method's repository can be found [here](https://github.com/jlakkis/CarDEC).

This repository is dedicated to providing the code used to perform all evaluations in the [CarDEC paper](https://www.biorxiv.org/content/10.1101/2020.09.23.310003v1). It includes code used to generate results for CarDEC, and for every competing method:

1. scVI
2. DCA + Combat
3. MNN
4. Scanorama
5. scDeepCluster

# General Flow

It is recommended the user proceeds as follows.

1. Clone this repository to their local machine
2. Download the data from Box.
3. Install all necessary packages.
4. Run all evaluations.
5. Run Rscripts to generate final plots.

## Clone this repository to your local machine

Clone this repository to your local machine using [the standard procedure](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository).

## Download the data from Box

Download the [data from Box](https://upenn.box.com/s/xlwg9e0vtj8a0xq6l87f2knwquclpjaw), and place them into the [currently empty data folder](https://github.com/jlakkis/CarDEC_Codes/tree/main/Data).

## Install all necessary packages

The user will need to install multiple packages: anaconda, two conda environments containing many dependencies, and a version of R >= 4.0

### Install Anaconda

First, install [Anaconda](https://www.anaconda.com/products/individual) if you do not already have it, so that you can access conda commands in terminal.

### Set up conda environments

Next, use [cardec.yml](https://github.com/jlakkis/CarDEC_Codes/blob/main/cardec.yml) and [cardec_alternatives.yml](https://github.com/jlakkis/CarDEC_Codes/blob/main/cardec_alternatives.yml) to set up the "cardec" and "cardec_alternatives" environments respectively.

To do this, simply cd in the cloned "CarDEC_Codes" repository. Once in this directory, run the following two commands.

```
$ conda env create -f cardec.yml
$ conda env create -f cardec_alternatives.yml
```

### Install R

Lastly, install a [version of R](https://www.r-project.org/). It is highly recommended that the user installs R version >= 4.0. [Rstudio](https://rstudio.com/products/rstudio/download/) is also reccomended for installation, but not required.

## Run all evaluations

Next, it is recommended that the user run all of the evaluation notebooks. The user should activate either the cardec or cardec_alternatives environment before opening jupyter to run the python notebooks. This is necessary because these two environments have "nb_conda_kernels" installed, which will allow the user to switch anaconda environments in the jupyter app. The following command will activate the cardec environment.

```
$ conda activate cardec
```

Then, open jupyter. The user can use either jupyter notebook or jupyter lab. The following command will open jupyter lab.

```
$ jupyter lab
```

### Run CarDEC Notebooks

It is recommended that the user first run the CarDEC notebook. Simply, open each of the following notebooks in jupyter. Make sure to set the activate conda kernel in jupyter to "cardec" and then run all cells. Repeat this for every notebook listed below.

1. [CarDEC Macaque.ipynb](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20CarDEC/CarDEC%20Macaque.ipynb)
2. [CarDEC Mouse Cortex.ipynb](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20CarDEC/CarDEC%20Mouse%20Cortex.ipynb)
3. [CarDEC Mouse Retina.ipynb](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20CarDEC/CarDEC%20Mouse%20Retina.ipynb)
4. [CarDEC PBMC.ipynb](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20CarDEC/CarDEC%20PBMC.ipynb)
5. [CarDEC Pancreas.ipynb](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20CarDEC/CarDEC%20Pancreas.ipynb)
6. [CarDEC Liver Runtime.ipynb](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20CarDEC/CarDEC%20Liver%20Runtime.ipynb)
7. [CarDEC Mouse Cortex-SCT.ipynb](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20CarDEC/CarDEC%20Mouse%20Cortex-SCT.ipynb)
8. [CarDEC Pancreas Revisions.ipynb](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20CarDEC/CarDEC%20Pancreas%20Revisions.ipynb)
9. [CarDEC_monocyte.ipynb](Evaluations of Competing Methods/script_for_monocyte/CarDEC_monocyte.ipynb)

### Run MNN R Scripts

Next, it is recommended that the user run all scripts to evaluate MNN. For each file in the list below, the user should open R (or Rstudio), and execute the script.

1. [MNN_Cortex.R](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20Competing%20Methods/MNN_Cortex.R)
2. [MNN_Cortex_HVG.R](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20Competing%20Methods/MNN_Cortex_HVG.R)
3. [MNN_Liver_Runtime.R](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20Competing%20Methods/MNN_Liver_Runtime.R)
4. [MNN_PBMC.R](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20Competing%20Methods/MNN_PBMC.R)
5. [MNN_PBMC_HVG.R](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20Competing%20Methods/MNN_PBMC_HVG.R)
6. [MNN_Pancreas.R](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20Competing%20Methods/MNN_Pancreas.R)
7. [MNN_Pancreas_HVG.R](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20Competing%20Methods/MNN_Pancreas_HVG.R)
8. [MNN_Retina.R](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20Competing%20Methods/MNN_Retina.R)
9. [MNN_Retina_HVG.R](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20Competing%20Methods/MNN_Retina_HVG.R)
10. [MNN_sampleMacaque.R](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20Competing%20Methods/MNN_sampleMacaque.R)
11. [MNN_sampleMacaque_HVG.R](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20Competing%20Methods/MNN_sampleMacaque_HVG.R)
12. [MNN_monocyte.R](Evaluations of Competing Methods/script_for_monocyte/MNN_monocyte.R)


### Run Other Method Python Notebooks

In the next step, the user should run the python notebooks to evaluate all methods other than CarDEC and MNN. Simply, open each of the following notebooks in jupyter. Make sure to set the activate conda kernel in jupyter to "cardec_alternatives" and then run all cells. Repeat this for every notebook listed below.

1. [Competing Methods Macaque.ipynb](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20Competing%20Methods/Competing%20Methods%20Macaque.ipynb)
2. [Competing Methods Mouse Cortex.ipynb](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20Competing%20Methods/Competing%20Methods%20Mouse%20Cortex.ipynb)
3. [Competing Methods Mouse Retina.ipynb](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20Competing%20Methods/Competing%20Methods%20Mouse%20Retina.ipynb)
4. [Competing Methods PBMC.ipynb](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20Competing%20Methods/Competing%20Methods%20PBMC.ipynb)
5. [Competing Methods Pancreas.ipynb](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20Competing%20Methods/Competing%20Methods%20Pancreas.ipynb)
6. [Competing Methods for monocyte](Evaluations of Competing Methods/script_for_monocyte/README.md) 
7. [DCA Liver Runtime.ipynb](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20Competing%20Methods/DCA%20Liver%20Runtime.ipynb)
8. [Scanorama Liver Runtime.ipynb](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20Competing%20Methods/Scanorama%20Liver%20Runtime.ipynb)
8. [scVI Liver Runtime.ipynb](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20Competing%20Methods/scVI%20Liver%20Runtime.ipynb)

**Remark:** The `Competing Methods for monocyte` is a folder.  All reproducing codes related to monocyte dataset can be found in this folder. 

### Run CV Score Notebooks

Lastly, the user should run the python notebooks used to generate the coefficient of variation plots demonstrated in many of the CarDEC paper's figures. Simply, open each of the following notebooks in jupyter. Make sure to set the activate conda kernel in jupyter to "cardec" and then run all cells. Repeat this for every notebook listed below.

1. [Batch Calibration Tests Cortex.ipynb](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20CarDEC/Batch%20Calibration%20Tests%20Cortex.ipynb)
2. [Batch Calibration Tests Macaque.ipynb](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20CarDEC/Batch%20Calibration%20Tests%20Macaque.ipynb)
3. [Batch Calibration Tests Mouse Retina.ipynb](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20CarDEC/Batch%20Calibration%20Tests%20Mouse%20Retina.ipynb)
4. [Batch Calibration Tests PBMC.ipynb](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20CarDEC/Batch%20Calibration%20Tests%20PBMC.ipynb)
5. [Batch Calibration Tests Pancreas.ipynb](https://github.com/jlakkis/CarDEC_Codes/blob/main/Evaluations%20of%20CarDEC/Batch%20Calibration%20Tests%20Pancreas.ipynb)

## Run Rscripts to Generate Final Plots

This last step is purely optional. In the previous steps, all analysis was completed. This final step involves using Rscripts to generate final figures. These Rscripts do not perform any actual analysis, they are simply used in order to generate prettier plots than Python for the paper. For example, all UMAP plots in the paper were generated by running all analysis in Python, exporting the computed UMAP coordinates to a csv file, and then reading this csv into R to build a prettier UMAP plot using ggplot2.

If the user wishes to generate the final plots, they just need to open each folder and run any Rscripts they find. These Rscripts should run in under 30 seconds each since they just read in small csv files and generate UMAP plots. The scripts have names like "figure_make.R", "figure_make_HVGo.R", "figure_make_bybatch.R", etcetra. A few figure folders will not contain Rscripts, which means that no R postprocessing was done to generate final figures.
