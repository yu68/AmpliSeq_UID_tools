# use Levenshtein distance algorithm for match 
python extract_UID.py -q test_input.fastq -I identifier.txt -M 1 -o test_output

# use pairwise-alignment for match
python extract_UID.py -q test_input.fastq -I identifier.txt -o test_output2
