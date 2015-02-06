# AmpliSeq\_UID\_tools

---

#### extract\_UID.py
Extract UID sequences from the read with structure "UU...UUII...IIRRRR...RRRR" or "RRRR...RRRRII...IIUU...UU". Trim the Identifier (II...II) and UID (UU...UU), keep the actual reads (RRRR...RRRR). Then ths UID information are stored in the read name information.

__input:__
```
@IM\_READ\_NAME 1:N:0:12 
UU...UUII...IIRRRR...RRRR
+
&%^^&$#%^##&&^@I^%#^%(($!
```

__output:__
```
@IM\_READ\_NAME\_UU...UU 1:N:0:12
RRRR...RRRR
+
@I^%#^%(($!
```

__Usage:__
```
usage: extract_UID.py [-h] [-f | -q] [-I IDENTIFIER] [-l UIDLEN] [-H]
                      [-m MAX_MIS] [-o OUTPUT]
                      input

Extract UID sequences and put it in fastq/a name

positional arguments:
  input                 input fastq/fasta file 1 with UID and identifiers

optional arguments:
  -h, --help            show this help message and exit
  -f, --fasta           add this option for fasta input file
  -q, --fastq           add this option for fastq input file
  -I IDENTIFIER, --identifier IDENTIFIER
                        identifier sequence file
  -l UIDLEN, --uidLen UIDLEN
                        length of UID, will search for +-1 basepair,
                        default=10
  -H, --head            set this if uid is at teh beginning of the reads,
                        otherwise, uid is at the end of reads
  -m MAX_MIS, --max_mis MAX_MIS
                        max(mismatch+indel) allowed for identifier match,
                        otherwise move reads into 'unassigned' file. default:
                        2
  -o OUTPUT, --output OUTPUT
                        output fastq/fasta file name

require: Bio

```
