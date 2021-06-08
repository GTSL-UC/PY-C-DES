# PY-C-DES
PY-C-DES VERSION 2.1 BUILDING GUIDE

1. PREREQUISITES
Obtaining and running PY-C-DES depends on the following

git Python {Tested with: Python 3.7+} For easy editing consider downloading and using Anaconda Navigator and Spyder (Python 3.7)

2. OBTAINING PY-C-DES
The best way to obtain PY-C-DES is to clone the git repository from GitHub as follows:

git clone {{{ THIS DIRECTORY }}}

This will create the directory PY-C-DES in the current working directory.

3. DIRECTORY CONTENTS
In the top level of the PY-C-DES directory are the following files:

COMPRESSOR FILES -'8st-des' .tci, .tcs, .igv details an 8-stage compressor. PY-C-DES takes the information from these files to calculate the compressor and fluid properties.

            To change the desired file to run select the 'PY-C-DES V.2.1.py'
            and adjust the input file names in the source code. 
            
            Later versions will allow the user to type the input file name into
            the terminal to choose which input files to run. 
AIR_XXX.csv - Tables that detail the fluid properties of air throughout a wide spectrum of conditions.

CH4_XXX.csv - Tables that detail the fluid properties of methane throughout a wide spectrum of conditions.

CO2_XXX.csv - Tables that detail the fluid properties of carbon dioxide throughout a wide spectrum of conditions.

README - This file.

Copyright.txt - Contains information pertaining to the use and distribution of PY-C-DES

License.txt - PY-C-DES license information
