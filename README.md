# BioBART
b7001699_csc8391

## Description
This is an expansion on protein annotation program that utilised the output of automatic annotation program and generate descriptive text about the annotations that were  completed.
Usage of VRAM was tested to be ~12gb on an A100 but might peak at around ~14gb.

## Commands and Usage
The application must and take a general Interpro output in tsv form and a FASTA file of the corresponding protein sequence.

The application has the following commands:
* `-i` `--input` - Path to the Interpro tsv file.
* `-f` `--fasta` - Path to the corresponding fasta file.
* `-l` `--long_text` - Default would only generate one liner, by incuding this command, output will be 2 to 3 line of texts.
* `-o`, `--output-file` - Path or name of the output file, default would be BioBART.fasta.

## Examples
There are included test file of test.tsv and test.fasta. Please ensure package in `requirements.txt` are fulfilled before run.
* `python ./bio_bart_test_copy.py -i ./test.tsv -f ./test.fasta` - Should output a `biobart.fasta` with descriptive test of the annotation along with protein sequence.