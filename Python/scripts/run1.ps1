# PI.ps1
python .\csv_extractor.py --cdip .\ncFiles\067.20201225_1200.20201227_0400.nc
python .\preprocessor.py --group=Meta --group=Wave --group=XYZ ".\ncFiles\067.20201225_1200.20201227_0400.nc"
python .\Driver.py ".\ncFiles\067.20201225_1200.20201227_0400_output.nc"
