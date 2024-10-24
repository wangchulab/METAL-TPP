#fitting curve, get scale parameters
python ../scripts/iter_fit_rawdata.py "37.0 41.0 45.0 48.0 51.0 54.0 58.0 62.0" control.csv experimental.csv

#parsing all data and compare curves
python ../scripts/parse_draw_CSVs.py "37.0 41.0 45.0 48.0 51.0 54.0 58.0 62.0" scale.txt > results.txt
