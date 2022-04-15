# Overview
This is the python portion of our project. These programs serve as a data processing pipeline. 

## Prerequisites 
Python, numpy, pandas, matplotlib, xarray, and netCDF4 

## File overview
1) `csv_extractor.py`: Reads .NC file, outputs several CSV files. In order to get testing data, we can use .NC files from CDIP. This program extracts the raw acceleration data and other information to simulate what we would get from the gliders SD card. <br />

2) `preprocessor.py`: Finds CSV files in the ncFiles folder, parses it, and stores it in nc file. This would normally be the first step after getting the CSV files from the gliders SD card. <br />

3) `Driver.py`: This is the driver program for all the calculations. This takes in a .NC file, preforms all the relevant calculations, then outputs a new filled out .NC file. <br />

`getdata.py` and `precessdata.py`: These two programs are helper functions for the Driver.py program. They have various functions reguarding extracting the data from the source .NC file and calculating all our coefficients<br />

`tests\test_non_directional_values.py` and `tests\test_directional_values.py`: These two programs are testing functions for Driver.py's results. Given a .NC file from CDIP, and the resulting .NC file from Driver.py, we can compare our calculations with CDIP. <br />


## Tutorial:
1) Extract CSV files from CDIP NC file with csv_extractor.py (optional if we are working with cdip data)<br /> 
&nbsp;-This creates several CSV files which simulates what we would get from the glider<br /> 
&nbsp;-Note: The --cdip flag means we need to take the step in converting cdips displacement data to acceleration data<br /> 
> python .\csv_extractor.py --cdip .\ncFiles\067.20201225_1200.20201225_1600.nc

2) Run preprocessor.py on the csv files that were created from the previous step or collected from the gliders SD card. This result in a file named filename_output.nc<br /> 
> python .\preprocessor.py --group=Meta --group=Wave --group=XYZ ".\ncFiles\067.20201225_1200.20201225_1600.nc" 

3) Run the driver.py on this _output.nc file. This will now result in a file named filename_output_proccessed.nc <br /> 
> python .\Driver.py ".\ncFiles\067.20201225_1200.20201225_1600_output.nc"

4) Testing (Optional): <br />
-Compare non-directional data of the fully processed file with CDIPs original file<br /> 
> python .\tests\test_non_directional_values.py --graph ".\ncFiles\067.20201225_1200.20201225_1600_output_processed.nc" ".\ncFiles\067.20201225_1200.20201225_1600.nc"

-Compare directional spectra of the fully processed file with CDIPs original file<br /> 
> python .\tests\test_directional_values.py --graph ".\ncFiles\067.20201225_1200.20201225_1600_output_processed.nc" ".\ncFiles\067.20201225_1200.20201225_1600.nc"