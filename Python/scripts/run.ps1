# PI.ps1
git clean -fx
$file=".\ncFiles\067.20201225_1200.20201227_0400.nc"
python .\csv_extractor.py --cdip $file
git clean -fx
