# chemoinformaticsTools
Scripts to analyse chemoinformatics software's output, collect data from the web...etc.

## :red_circle: diversity_clustering

A Python script to select a diverse subset of compounds from a library (.CSV or .xlsx file). It is based on MACCS keys fingerprints or Morgan circular fingerprints, and uses a MaxMin or hierachical clustering algorithm to select the most diverse compounds.

The following non-standard python libraries are needed for the script to function correctly:
* pandas
* numpy
* rdkit

## :red_circle: jatoonSuccessRate

[JATOON](http://joao.airesdesousa.com/jatoon/) (Java Tools for Neural Networks) was started in 2001 as a system of Java applets for training and applying neural networks.

jatoonSuccessRate.sh is a Bash script designed to :
* Read the Tools/Predict output from jatoonSOM for Kohonen self-organizing maps or counterpropagation neural networks
* Output the number of classes correctly predicted.

## :red_circle: webScraper

webScraper is a Python script designed to collect data from several webpages and output them in a CSV format, with Data-Mining in mind. It is capable of multi-threading for faster data collection.

The present version is tailored to target [The Good Scents Company](http://www.thegoodscentscompany.com/) data, extracting the following properties : 'CAS Number','Name','SMILES','InChIKey','Molecular Weight','Odor Type','Odor Strength','Odor Description', and 'Taste Description', when available.  
 It is then up to the user to use data curation techniques to clean the data fetched from the website.
 
The following non-standard python libraries are needed for the script to function correctly:
* lxml
* requests
* pandas
* progressbar
