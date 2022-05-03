# PI.ps1
python .\csv_extractor.py --cdip .\ncFiles\067.20220102_0000.20220102_0600.nc
python .\preprocessor.py --group=Meta --group=Wave --group=XYZ ".\ncFiles\067.20220102_0000.20220102_0600.nc"
python .\Driver.py ".\ncFiles\067.20220102_0000.20220102_0600_output.nc"
python .\tests\test_non_directional_values.py --graph ".\ncFiles\067.20220102_0000.20220102_0600_output_processed.nc" ".\ncFiles\067.20220102_0000.20220102_0600.nc"
python .\tests\test_directional_values.py --graph ".\ncFiles\067.20220102_0000.20220102_0600_output_processed.nc" ".\ncFiles\067.20220102_0000.20220102_0600.nc"
