# PY-C-DES
PY-C-DES VERSION 2.5, 2.6 BUILDING GUIDE

INTRO: PY-C-DES

   PY-C-DES takes input properties (inlet conditions, and stage by stage info such as temperature rise, certain type of fluid, alpha angles, mach #'s, ect...) and simulates a compressor design to find the compressor and fluid properties using mean-line calculations. The outputs of the code provide a wide array of properties and flow values at each rotor and stator; as well as a overall compressor design properties such as adiabatic efficiency and total pressure ratio along with 2D geometric coordinates of the blades, tip, and hub. It can be run both with the perfect gas assumption and with the real gas properties -- hence the addition of fluid folders and property tables in this directory. It is meant to be used for the starting stage of preliminary compressor design; but could also be used as an educational tool. 

1. PREREQUISITES
            Obtaining and running PY-C-DES depends on the following

            git 
            Python 
            {Tested with: Python 3.7+} 
            For easy editing consider downloading and using Anaconda Navigator and Spyder (Python 3.7)
            Numpy
            Matplotlib
            Scipy

2. OBTAINING PY-C-DES
            The best way to obtain PY-C-DES is to clone the git repository from GitHub as follows:

            git clone https://github.com/GTSL-UC/PY-C-DES.git

     This will create the directory PY-C-DES in the current working directory.
     
     Downloading the zipped folder and extracting it into the desired directory, 
            then running the code works as well.

3. DIRECTORY CONTENTS
            In the top level of the PY-C-DES directory are the following files:

COMPRESSOR INPUT FILES -'8st-des' .tci, .tcs, .igv details an 8-stage compressor. PY-C-DES takes the information from these files to calculate the compressor and fluid properties.
            
   8st-des.tci -- details overall compressor properties like # of stages, radius option, mass flow rate, rotation speed, inlet conditions, and if the code will us the perfect gas assumption or utilize real gas properties. 
   
   8st-des.tcs -- details compressor properties for each stage such as inlet flow angles, mach numbers, enthalpy rise or temperature rise (depending on the setting) and much more. 
   
   8st-des.igv -- optional -- PY-C-DES will automatically pick up if the IGV file is in the directory and will utilize it if it is there. 
   It details the Inlet Guide Vein properties
            
AIR - Folder with tables that detail the fluid properties of air throughout a wide spectrum of conditions.

CH4 - Folder with tables that detail the fluid properties of methane throughout a wide spectrum of conditions.

CO2 - Folder with tables that detail the fluid properties of carbon dioxide throughout a wide spectrum of conditions.

Other full table properties such as O2 and H2 follow the same format. 

README - This file.

Copyright.txt - Contains information pertaining to the use and distribution of PY-C-DES

License.txt - PY-C-DES license information

4. RUNNING PY-C-DES
To run PY-C-DES, make sure you're in the current directory where these files were downloaded. 

In the command prompt, type:
            
            python PyCDes_V_2_6.py 8st-des
  
  Where 8st-des specifies which compressor inputs to use. PY-C-DES will find the .tci, .tcs, and igv if they are in the directory. 
  
 If using Spyder, then specify the input file prefix with Run > configuration per file >  then select 'command line option' and type '8st-des'
 Hitting the green run button runs the code. 

5. OUTPUTS OF PY-C-DES

After successfully running PY-C-DES, new output files will be created in the current directory. The files will be named after the fluid type and the date the code was run. For instance, if the fluid is 'AIR' and it's run on 8/29/2020; then the output files will be named 'AIR_08_29_20....' This name will be used for this ReadMe example. 

The types of files exported are: 
                
   AIR_08_29_20.txt -- The main export, details the fluid properties at each rotor and stator throughout the compressor. Such as enthalpy, entropy, viscosities, temperatures, mach #s, pressures, and turbomachinery coefficients. (And Many More!)
   
   AIR_08_29_20_Geometry.txt -- details geometry of the compressor with flow angles, radii lengths, chord lengths, number of blades for each station, aspect ratios, and more. It outlines the geometry values on a Per Stage Basis. 
   
   AIR_08_29_20_GeometryVectors.txt -- details the geometry in vector form, for easier reading by a computer and further design/analysis. 
   
   AIR_08_29_20.twal -- output in line with the T-AXI optimization suite. Exports the radius and axial coordinates from a 2d cross section of the compressor. 
      
   AIR_08_29_20.tstk -- output in line with the T-AXI optimization suite. Exports blade leading edge and trailing edge coordinates at both the hub and tip; as well as other properties of the compressor like loss and blockage coefficients. 
   
   AIR_08_29_20.tinf -- output in line with the T-AXI optimization suite. Exports the Number of Blades per each blade row, as well as the corresponding RV-Theta velocity. 
   
   AIR_08_29_20.pdf -- The fun output. Details many compressor and fluid properties in easy to read figures. Over 50 figures output for easy analysis of the design. 
   
   NOTE --> The specific heat and gamma graphs use the real gas calculations for the values, regardless of which fluid is selected ('real' or 'perfect gas' setting being selected.) Hence, why they still vary even if the 'perfect gas' is selected. This is done to demonstrate how the specific heat or gamma would vary in the current compressor conditions. However, all the equations for the fluid properties use the gamma or specific heat that corresponds with the fluid setting, so all other properties are in line with the fluid type selected. 

6. EDITING THE COMPRESSOR in PY-C-DES
   
   Changing the compressor properties in the '8st-des' files is very easy. Open the '8st-des' input files with any editor (notepad works just fine) and change any value or setting your heart desires for your current work or project; make sure to save the file before running. If there's going to be >10 stages, make sure there is data for >10 stages in the TCS file first before running, otherwise there will be an error.
   
   One of the main parameters to change is the 'radius type' in the .tci file. If the constant hub is selected then it keeps the hub constant throughout the entire compressor; the tip and pitch selection work that way as well. The setting with most variable outputs is the '4-Arbitrary' setting. This makes the radii of the compressor for each stage depend on the .tcs file radii values - so if the arbitrary setting is desired, make sure to look over the radii values in the .tci file, so that they make sense and the compressor isn't some monstrosity. 
   
   
   
   
