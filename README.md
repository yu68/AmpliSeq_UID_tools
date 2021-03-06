# AmpliSeq\_UID\_tools

---

#### extract\_UID.py
Extract UID sequences from the read with structure "UU...UUII...IIRRRR...RRRR" or "RRRR...RRRRII...IIUU...UU". Trim the Identifier (II...II) and UID (UU...UU), keep the actual reads (RRRR...RRRR). Then the UID information are stored in the read name information.

__input:__
```
@IM_READ_NAME 1:N:0:12 
UU...UUII...IIRRRR...RRRR
+
&%^^&$#%^##&&^@I^%#^%(($!
```
  or
```
@IM_READ_NAME 1:N:0:12
RRRR...RRRRII...IIUU...UU
+
&%^^&$#%^##&&^@I^%#^%(($!
```

__output:__
```
@IM_READ_NAME_UU...UU 1:N:0:12
RRRR...RRRR
+
@I^%#^%(($!
```
  or 
```
@IM_READ_NAME_UU...UU 1:N:0:12
RRRR...RRRR
+
&%^^&$#%^##
```



__Usage:__
```
usage: extract_UID.py [-h] [-f | -q] [-I IDENTIFIER] [-l UIDLEN] [-H]
                      [--keep_identifier] [-m MAX_MIS] [-M METHOD] [-o OUTPUT]
                      input

Extract UID sequences and put it in fastq/a name

positional arguments:
  input                 input fastq/fasta file with UID and identifiers

optional arguments:
  -h, --help            show this help message and exit
  -f, --fasta           add this option for fasta input file
  -q, --fastq           add this option for fastq input file
  -I IDENTIFIER, --identifier IDENTIFIER
                        identifier sequence file
  -l UIDLEN, --uidLen UIDLEN
                        length of UID, will search for +-1 basepair,
                        default=10
  -H, --head            set this if uid is at the beginning of the reads,
                        otherwise, uid is at the end of reads
  --keep_identifier     set to keep the identifier sequences in the output
                        fastq file
  -m MAX_MIS, --max_mis MAX_MIS
                        max(mismatch+indel) allowed for identifier match,
                        otherwise move reads into 'unassigned' file. default:
                        2
  -M METHOD, --Method METHOD
                        Method for alignment(1-2): 1. Levenshtein distance
                        algorithm; 2. pairwise Alignment; Default: 2
  -o OUTPUT, --output OUTPUT
                        output fastq/fasta file name

require: Bio

```

__Example:__
```
# use Levenshtein distance algorithm for match
python extract_UID.py -q test_input.fastq -I identifier.txt -M 1 -o test_output

# use pairwise-alignment for match
python extract_UID.py -q test_input.fastq -I identifier.txt -o test_output2

```

__Explanations:__
* the identifier sequence file contains identifier sequences one per line. The program will check the alignment between read sequence and each identifier sequence and find the alignment with smallest mismatchs and indels. All the identifier sequences need to have same length
* If the location of best alignment with identifiers is close to the actual UID (with UID length specified) and identifier locations, UIDs are extracted.
* '-H' need to be specified if UID is at the begining of the reads (case: "UU...UUII...IIRRRR...RRRR"), otherwise (case:"RRRR...RRRRII...IIUU...UU") don't set '-H'
* '--keep_identifier' need to be specified if you want to keep identifier sequences (II...II) untrimmed. 
* two alignment methods are included: 1. Levenshtein distance algorithm; 2. pairwise Alignment;
* If the best alignment has mismatch and indel number larger than the number specified by '-m', we cannot locate the identifier and the original reads are stored in "unassign_<output>.fastq/a" file.

