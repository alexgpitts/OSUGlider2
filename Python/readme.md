# Tutorial

These programs serve as a data processing pipeline
## For example:
1) Extract CSV files from CDIP NC file with csv_extractor.py (optional if we are working with cdip data)<br /> 
&nbsp;-This creates several CSV files which simulates what we would get from the glider<br /> 
&nbsp;-Note: The --cdip flag means we need to take the step in converting cdips displacement data to acceleration data<br /> 
> python .\csv_extractor.py --cdip ".\ncFiles\067.20201225_1200.20201225_1600.nc"

2) Run preprocessor.py on the csv files that were created from the previous step or collected from the gliders SD card. This result in a file named filename_processed.nc<br /> 
>python .\preprocessor.py --group=Meta --group=Wave --group=XYZ ".\ncFiles\067.20220102_0000.20220102_0600.nc"  