# PI.ps1
git clean -fx
$file=".\ncFiles\067.20201225_1200.20201227_0400"
python .\csv_extractor.py --cdip $file".nc"
python .\preprocessor.py --group=Meta --group=Wave --group=XYZ $file".nc"
python .\Driver.py $file"_output.nc"
python .\tests\test_vals_v2.py $file"_output_processed.nc" $file".nc"

