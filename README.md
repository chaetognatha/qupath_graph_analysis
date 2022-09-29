# graph analysis with network x

dependencies: pandas numpy matplotlib libpysal networkx seaborn scipy
new dependencies: alive-progress
optional dependencies: geopandas

# how to run:

python nwx_analyze.py -h

## example
### NOTE: you now need to specify a path to where you keep your images!!!
python nwx_analyze.py -i ~/Documents/Ovarial_22/ML_TMA1/obj_class_TMA1.csv -p CD45:PANCK --tiff_dir ~/Documents/Ovarial_22/imj_out_1

# Notes

## Regarding class1:class2
while developing the script, class1 was assumed to be immune cells and class2 cancer cells, this means for example that group degree centrality is calculated for whatever is on the left side!

# QUICKSTART
* Make sure you have conda!
* Output the results of your cell detections in QuPath to the default format
* Open the output tsv to determine what decimal point QuPath was using (if it is not . you need to specify)
* Create new conda environment with all dependencies:

```bash
conda create -n graph_analysis -c conda-forge -c bioconda pandas numpy matplotlib libpysal networkx seaborn scipy geopandas alive-progress
```

* Activate your environment (start here if you already set the environment up)
```conda activate graph_analysis ```

* Clone this repo and navigate to the local directory with the conda prompt

* Get the help text for constructing your CLI command

```python nwx_analyze.py --help ```

# Troubleshooting:
If you get an error that seems to indicate a package is missing, first check that your conda environment is properly set up and that you have activated it - then if the problem persists write an issue on github!
