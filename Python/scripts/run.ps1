# PI.ps1
git clean -fxq
$file=".\ncFiles\067.20201225_1200.20201227_0400"
# $file=".\ncFiles\067.20201225_1200.20201225_1600"
python .\csv_extractor.py --cdip $file".nc"
python .\preprocessor.py --group=Meta --group=Wave --group=XYZ $file".nc"
python .\Driver.py $file"_output.nc"
python .\tests\test_vals_v2.py $file"_output_processed.nc" $file".nc"
git clean -fxq

