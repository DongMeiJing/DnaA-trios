# DnaA-trios
DnaA-trio is an important element in bacterial replication origins. This project can be used to search DnaA-trios in bacterial replication origins.

## Implementation
This project is implemented Python 3.7.7.

## Usage
```sh
git clone https://github.com/DongMeiJing/DnaA-trios.git
cd ./DnaA-trios/Search_DnaA_trios
pip install -r requirements.txt
python Search_DnaA_trios.py
```

## Input Files
The input data for this project are bacterial replication origins sequences (DoriC_oric_sequences.csv) and DnaA-box motif (Species_specific_DnaA_Box_motif.csv and Epsilonproteobacteria_DnaA_Box_motif.csv) located in the folder Search_DnaA_trios. These data come from the DoriC database (http://tubic.org/doric).

## Output Files
This project finally outputs a file named 'DnaA_trios.csv'.
