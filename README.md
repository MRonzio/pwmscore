# PWM matrix score calculator. 
Simple PWM score calculator based on [1].


## Usage
```bash
CalcScoreJaspar copy.py [-h] -M JASPAR_MATRIX -o TEST_OLIGO [-s STRAND]

JASPAR matrix score calculator

optional arguments:
  -h, --help            show this help message and exit
  -M JASPAR_MATRIX, --matrix JASPAR_MATRIX
                        jaspar PWM in jaspar format
  -o TEST_OLIGO, --oligo TEST_OLIGO
                        target name
  -s STRAND, --strand STRAND
                        oligo strand. Options: "pos" for positive, "neg" for negative. Default "pos".
```

## Example with test data
```CalcScoreJaspar.py -M test_data/MAA060_1_JASPAR2016 -o CTCAGCCAATCAGCGC -s pos```

## Requirements
python3
pandas


## Reference
[1] https://doi.org/10.1093/nar/gkt448
