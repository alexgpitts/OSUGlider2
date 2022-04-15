# Ocean wave measurements (CDIP Testing environment)

This is legacy code our team used for research on how to calculate things from CDIP files. THe program takes in a CDIP .NC file and process it with various calculation methods.

This program needs numpy, matplotlib, netCDF4, and xarray installed to run. To display different graphs, change the plotting perimeter boolean values towards the top of the driver file.

This file takes command lines<br /> 
ie: cdip_driver.py [--welch | --banding] [--hann | --boxcar] [--raw] [--ds] [--graph] -o [--output] nc [nc ...] <br />
where: <br />
    1) --welch, --banding chooses the calculation method <br />
    2) --hann or --boxcar selects the windowing method <br />
    3) --norm shows the normal fft using built in numpy fft method <br />
    4) --raw plots the raw acceleration data <br />
    5) --ds plots the directional spectra coefficients compared with CDIPs calculations <br />
    6) --graph chooses if you want to plot the resulting PSD calculated from the chosen calculation method <br />
    7) -o ".\outputFolder" specifies the output folder for outputted .nc file to go. By default, the outputted data is placed in the same folder as the .nc file that was read from <br />
    8) nc is the .nc file you wish to process. Takes a sting as input <br />


Current Status: <br />
As of now, we have most of the calculations performed on some test data from CDIP. 

Some of our next steps: 
1) Fix normalization for the welch method using a boxcar window
2) Create a testing suite for running large amounts of CDIP data through the program. 


Example Output: `python .\cdip_driver.py --banding --graph --raw --norm --ds ".\ncFiles\067.20201225_1200.20201225_1600.nc"`
![builds](../ProjectImages/python_output.png?raw=true)
Alternativly to run welch method: `python .\cdip_driver.py --welch --hann --raw --norm --ds --graph ".\ncFiles\067.20201225_1200.20201225_1600.nc"`  