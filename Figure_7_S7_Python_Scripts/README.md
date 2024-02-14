# Overview

To create the resulting `all_data.tsv` file in the `Analysis` directory, analysis from the different TIFF and CZI files has to be performed and integrated. Please first install the Python package requirements in `requirements.txt` before proceeding with any analysis - i.e. `pip install -U requirements.txt`. Note that this process necessitates an [Anaconda](https://docs.anaconda.com/free/anaconda/install/index.html) environment or a [virtual Python](https://docs.python.org/3/library/venv.html) environment. Please choose the method that you're most comfortable with.

The general analysis pipeline is a 2-step process:

1. Run `./calculate_pH.py` with the corresponding command line arguments (see next section), which generates several files for each Image File and Mask dataset pair. For multiple File and Mask pairs, this step must be performed for each set before proceeding to the next step.
2. Run `./aggregate_violin_plot.py` (no command line arguments are required), and all the data from the different runs will be merged into `all_data.tsv`.


# Example Runs of `calculate_pH.py` using the current data

Note that the following patch coordinates are illustrative and specific to the dataset in question. For a breakdown of options, please consult the output from `./calculate_pH.py -h`. These examples can be run from the command line (Unix, Linux, Mac, WSL) in the current directory. Navigate to that location, and run the commands below:

## NCL_Pre

```
# Field 1
./calculate_pH.py \
-c ../Figure_7_S7_Python_Data/snarf4.txt \
-w ../Figure_7_S7_Python_Data/all_wavelengths.txt \
-z ../Figure_7_S7_Python_Data/231106_NCL_Pre_Field1.czi \
-m ../Figure_7_S7_Python_Data/231106_NCL_Pre_Field1_Mask.tif \
-pc 250 350 150 250 \
-pc 25 125 425 525 \
-pc 400 475 450 525

# Field 2
./calculate_pH.py \
-c ../Figure_7_S7_Python_Data/snarf4.txt \
-w ../Figure_7_S7_Python_Data/all_wavelengths.txt \
-z ../Figure_7_S7_Python_Data/231106_NCL_Pre_Field2.czi \
-m ../Figure_7_S7_Python_Data/231106_NCL_Pre_Field2_Mask.tif \
-pc 75 175 100 200 \
-pc 425 500 225 300 \
-pc 450 525 425 500
```

## NCL_Pre_LYAR

```
# Field 1
./calculate_pH.py \
-c ../Figure_7_S7_Python_Data/snarf4.txt \
-w ../Figure_7_S7_Python_Data/all_wavelengths.txt \
-z ../Figure_7_S7_Python_Data/231106_NCL_Pre_LYAR_Field1.czi \
-m ../Figure_7_S7_Python_Data/231106_NCL_Pre_LYAR_Field1_Mask.tif \
-pc 275 375 375 475 \
-pc 275 375 375 475 \
-pc 325 425 375 475

# Field 2
./calculate_pH.py \
-c ../Figure_7_S7_Python_Data/snarf4.txt \
-w ../Figure_7_S7_Python_Data/all_wavelengths.txt \
-z ../Figure_7_S7_Python_Data/231106_NCL_Pre_LYAR_Field2.czi \
-m ../Figure_7_S7_Python_Data/231106_NCL_Pre_LYAR_Field2_Mask.tif \
-pc 250 350 50 150 \
-pc 275 375 75 175 \
-pc 25 100 200 275
```


## NPM1_Mat

```
# Field 1
./calculate_pH.py \
-c ../Figure_7_S7_Python_Data/snarf4.txt \
-w ../Figure_7_S7_Python_Data/all_wavelengths.txt \
-z ../Figure_7_S7_Python_Data/231106_NPM1_Mat_Field1.czi \
-m ../Figure_7_S7_Python_Data/231106_NPM1_Mat_Field1_Mask.tif \
-pc 150 250 125 225 \
-pc 50 150 50 150 \
-pc 250 350 225 325

# Field 2
./calculate_pH.py \
-c ../Figure_7_S7_Python_Data/snarf4.txt \
-w ../Figure_7_S7_Python_Data/all_wavelengths.txt \
-z ../Figure_7_S7_Python_Data/231106_NPM1_Mat_Field2.czi \
-m ../Figure_7_S7_Python_Data/231106_NPM1_Mat_Field2_Mask.tif \
-pc 150 250 125 225 \
-pc 50 150 50 150 \
-pc 250 350 225 325
```


## NdelLPre_L0p2uM

```
Field 1
./calculate_pH.py \
-c ../Figure_7_S7_Python_Data/snarf4.txt \
-w ../Figure_7_S7_Python_Data/all_wavelengths.txt \
-z ../Figure_7_S7_Python_Data/230615_NdelLPre_L0p2uM_Field1.czi \
-m ../Figure_7_S7_Python_Data/230615_NdelLPre_L0p2uM_Field1_Mask.tif \
-pc 60 110 240 290 \
-pc 700 750 320 370 \
-pc 530 580 400 450

# Field 2
./calculate_pH.py \
-c ../Figure_7_S7_Python_Data/snarf4.txt \
-w ../Figure_7_S7_Python_Data/all_wavelengths.txt \
-z ../Figure_7_S7_Python_Data/230615_NdelLPre_L0p2uM_Field2.czi \
-m ../Figure_7_S7_Python_Data/230615_NdelLPre_L0p2uM_Field2_Mask.tif \
-pc 140 190 480 530 \
-pc 660 710 580 630 \
-pc 350 400 490 540
```