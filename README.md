# Advnced-Python
this is my portfolio of code for advanced computational biology
## Sequence Objects
```python
from Bio.Seq import Seq
```


```python
my_seq = Seq("GATCG")
```


```python
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G



```python
# here we are printing the length of each sequence
print(len(my_seq))
```

    5



```python
print(my_seq[0])
```

    G



```python
print(my_seq[4])
```

    G



```python
print(my_seq[2])
```

    T



```python
Seq("AAAA").count("AA")
```




    2




```python
my_seq = Seq("GGGGGAAAATTTTCCCTTCTCATAT")
```


```python
len(my_seq)
```




    25




```python
my_seq.count("G")
```




    5




```python
100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
```




    40.0




```python
from Bio.SeqUtils import gc_fraction
```


```python
my_seq = Seq("GGGGGAAAATTTTCCCTTCTCATAT")
```


```python
gc_fraction(my_seq)
```




    0.4




```python
# here we will slice fractions with start stop and stride
my_seq[4:12]
```




    Seq('GAAAATTT')




```python
my_seq[0::3]
```




    Seq('GGATTCCAT')




```python
my_seq[1::3]
```




    Seq('GGATCTTT')




```python
my_seq[2:3]
```




    Seq('G')




```python
# you can also run the sequence backwards with negative numbers
my_seq[::-1]
```




    Seq('TATACTCTTCCCTTTTAAAAGGGGG')




```python
str(my_seq)
```




    'GGGGGAAAATTTTCCCTTCTCATAT'




```python
fasta_format_string = ">Name\n%s\n" % my_seq
```


```python
print(fasta_format_string)
```

    >Name
    GGGGGAAAATTTTCCCTTCTCATAT
    



```python
seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
```


```python
seq1 + seq2
```




    Seq('ACGTAACCGG')




```python
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTTGCA")]
```


```python
spacer = Seq("N" *10)
```


```python
spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTTGCA')




```python
dna_seq = Seq("acgtACGT")
```


```python
dna_seq
```




    Seq('acgtACGT')




```python
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
dna_seq.lower()
```




    Seq('acgtacgt')




```python
#this will be false because of case sensitivity
"gtac" in dna_seq
```




    False




```python
dna_seq = dna_seq.upper()
```


```python
# This one will be true
"GTAC" in dna_seq
```




    True




```python
my_seq = Seq("GATCTATATATATAATCTCTCT")
```


```python
# this gives us the complementn strand of dna
my_seq.complement()
```




    Seq('CTAGATATATATATTAGAGAGA')




```python
# this gives us the reverse complement strand
my_seq.reverse_complement()
```




    Seq('AGAGAGATTATATATATAGATC')




```python
# this displays the ambiguity code for nucleotides
protein_seq = Seq("EVRNAK")
protein_seq.complement()
```




    Seq('EBYNTM')




```python
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
```


```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
template_dna = coding_dna.reverse_complement()
```


```python
template_dna
```




    Seq('CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT')




```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
messenger_rna = coding_dna.transcribe()
```


```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
template_dna.reverse_complement().transcribe()
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
# this can turn dna to rna and back to dna via reverse transcription
messenger_rna.back_transcribe()
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
# this shows translation (rna turning into proteins)
messenger_rna.translate()
```




    Seq('MAIVMGR*KGAR*')




```python
# we can also code genes from codon tables
coding_dna.translate(table="Vertebrate Mitochondrial")
```




    Seq('MAIVMGRWKGAR*')




```python
# you can also translate the nucleotides up to the first stop codon
coding_dna.translate(to_stop = True)
```




    Seq('MAIVMGR')




```python
coding_dna.translate(table = 2, to_stop=True)
```




    Seq('MAIVMGRWKGAR')




```python
# you can also specify what the stop codon is
coding_dna.translate(table = 2, stop_symbol = "!")
```




    Seq('MAIVMGRWKGAR!')




```python
# for situations when there isnt a typical stop codon such as certain bacteria, you can import sequence from NCBI
```


```python
gene = Seq("GTGAAAAAAGATGCAAGCTATATACGCGTATATCGGGAAATATGTCTCGCTAGTCGATCGATCGATCGATCGCTACGTAGCTACGGTAGCTAGCATCGATCGATCGATCGGATCGATCGATCGATCGTAGCTAGCTACGATGCTAC")
```


```python
gene.translate(table = "Bacterial")
```




    Seq('VKKDASYIRVYREICLASRSIDRSLRSYGS*HRSIDRIDRSIVASYDA')




```python
gene.translate(table = "Bacterial", to_stop = True)
```




    Seq('VKKDASYIRVYREICLASRSIDRSLRSYGS')




```python
from Bio.Data import CodonTable
```


```python
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```


```python
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
# we can print the tables for vizualization
print(standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
mito_table.stop_codons
```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
mito_table.start_codons
```




    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']




```python
seq = Seq("ACGT")
```


```python
# you can also test equivalence of 2 codons
"ACGT" == seq1
```




    True




```python
seq1 == "ACGT"
```




    True




```python
# sequences can be identified by its configuration:
unknown_seq = Seq(None, 10)
```


```python
unknown_seq
```




    Seq(None, length=10)




```python
len(unknown_seq)
```




    10




```python
# we can pull sequences from biopython
seq = Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length = 159344973)
```


```python
seq[1000:1020]
```




    Seq(None, length=20)




```python
seq[117512690:117512700]
```




    Seq('CCTGAATGTG')




```python
# this shows that you can have partial information of a chromosome without the entire length
seq[117512670:]
```




    Seq({13: 'TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT'}, length=41832303)




```python
seq = Seq("ACGT")
```


```python
undefined_seq = Seq(None, length =10)
```


```python
seq + undefined_seq + seq
```




    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)




```python
my_seq = Seq("GCCATTGTAATGGGCCGCTAAAGGGTGCCCGA")
```


```python
# we can create mutable sequences
from Bio.Seq import MutableSeq
```


```python
mutable_seq = MutableSeq(my_seq)
```


```python
mutable_seq
```




    MutableSeq('GCCATTGTAATGGGCCGCTAAAGGGTGCCCGA')




```python
# here we will add a mutation
mutable_seq.remove("T")
```


```python
mutable_seq
```




    MutableSeq('GCCATGTAATGGGCCGCTAAAGGGTGCCCGA')




```python
# and the reverse mutated sequence
mutable_seq.reverse()
```


```python
mutable_seq
```




    MutableSeq('AGCCCGTGGGAAATCGCCGGGTAATGTACCG')




```python
# if we want to change our gene and return to immutable:
new_seq = Seq(mutable_seq)
```


```python
new_seq
```




    Seq('GCCATGTAATGGGCCGCTAAAGGGTGCCCGA')




```python
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
```


```python
my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG"
```


```python
reverse_complement(my_string)
```




    'CTAACCAGCAGCACGACCACCCTTCCAACGACCCATAACAGC'




```python
transcribe(my_string)
```




    'GCUGUUAUGGGUCGUUGGAAGGGUGGUCGUGCUGCUGGUUAG'




```python
back_transcribe(my_string)
```




    'GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG'




```python
translate(my_string)
```




    'AVMGRWKGGRAAG*'

## Sequence Annotation
```python
from Bio.SeqRecord import SeqRecord
```


```python
# this will pull up the help file
help(SeqRecord)
```

    Help on class SeqRecord in module Bio.SeqRecord:
    
    class SeqRecord(builtins.object)
     |  SeqRecord(seq, id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=None, features=None, annotations=None, letter_annotations=None)
     |  
     |  A SeqRecord object holds a sequence and information about it.
     |  
     |  Main attributes:
     |   - id          - Identifier such as a locus tag (string)
     |   - seq         - The sequence itself (Seq object or similar)
     |  
     |  Additional attributes:
     |   - name        - Sequence name, e.g. gene name (string)
     |   - description - Additional text (string)
     |   - dbxrefs     - List of database cross references (list of strings)
     |   - features    - Any (sub)features defined (list of SeqFeature objects)
     |   - annotations - Further information about the whole sequence (dictionary).
     |     Most entries are strings, or lists of strings.
     |   - letter_annotations - Per letter/symbol annotation (restricted
     |     dictionary). This holds Python sequences (lists, strings
     |     or tuples) whose length matches that of the sequence.
     |     A typical use would be to hold a list of integers
     |     representing sequencing quality scores, or a string
     |     representing the secondary structure.
     |  
     |  You will typically use Bio.SeqIO to read in sequences from files as
     |  SeqRecord objects.  However, you may want to create your own SeqRecord
     |  objects directly (see the __init__ method for further details):
     |  
     |  >>> from Bio.Seq import Seq
     |  >>> from Bio.SeqRecord import SeqRecord
     |  >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
     |  ...                    id="YP_025292.1", name="HokC",
     |  ...                    description="toxic membrane protein")
     |  >>> print(record)
     |  ID: YP_025292.1
     |  Name: HokC
     |  Description: toxic membrane protein
     |  Number of features: 0
     |  Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF')
     |  
     |  If you want to save SeqRecord objects to a sequence file, use Bio.SeqIO
     |  for this.  For the special case where you want the SeqRecord turned into
     |  a string in a particular file format there is a format method which uses
     |  Bio.SeqIO internally:
     |  
     |  >>> print(record.format("fasta"))
     |  >YP_025292.1 toxic membrane protein
     |  MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
     |  <BLANKLINE>
     |  
     |  You can also do things like slicing a SeqRecord, checking its length, etc
     |  
     |  >>> len(record)
     |  44
     |  >>> edited = record[:10] + record[11:]
     |  >>> print(edited.seq)
     |  MKQHKAMIVAIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
     |  >>> print(record.seq)
     |  MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
     |  
     |  Methods defined here:
     |  
     |  __add__(self, other)
     |      Add another sequence or string to this sequence.
     |      
     |      The other sequence can be a SeqRecord object, a Seq object (or
     |      similar, e.g. a MutableSeq) or a plain Python string. If you add
     |      a plain string or a Seq (like) object, the new SeqRecord will simply
     |      have this appended to the existing data. However, any per letter
     |      annotation will be lost:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
     |      >>> print("%s %s" % (record.id, record.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(record.letter_annotations))
     |      ['solexa_quality']
     |      
     |      >>> new = record + "ACT"
     |      >>> print("%s %s" % (new.id, new.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNNACT
     |      >>> print(list(new.letter_annotations))
     |      []
     |      
     |      The new record will attempt to combine the annotation, but for any
     |      ambiguities (e.g. different names) it defaults to omitting that
     |      annotation.
     |      
     |      >>> from Bio import SeqIO
     |      >>> with open("GenBank/pBAD30.gb") as handle:
     |      ...     plasmid = SeqIO.read(handle, "gb")
     |      >>> print("%s %i" % (plasmid.id, len(plasmid)))
     |      pBAD30 4923
     |      
     |      Now let's cut the plasmid into two pieces, and join them back up the
     |      other way round (i.e. shift the starting point on this plasmid, have
     |      a look at the annotated features in the original file to see why this
     |      particular split point might make sense):
     |      
     |      >>> left = plasmid[:3765]
     |      >>> right = plasmid[3765:]
     |      >>> new = right + left
     |      >>> print("%s %i" % (new.id, len(new)))
     |      pBAD30 4923
     |      >>> str(new.seq) == str(right.seq + left.seq)
     |      True
     |      >>> len(new.features) == len(left.features) + len(right.features)
     |      True
     |      
     |      When we add the left and right SeqRecord objects, their annotation
     |      is all consistent, so it is all conserved in the new SeqRecord:
     |      
     |      >>> new.id == left.id == right.id == plasmid.id
     |      True
     |      >>> new.name == left.name == right.name == plasmid.name
     |      True
     |      >>> new.description == plasmid.description
     |      True
     |      >>> new.annotations == left.annotations == right.annotations
     |      True
     |      >>> new.letter_annotations == plasmid.letter_annotations
     |      True
     |      >>> new.dbxrefs == left.dbxrefs == right.dbxrefs
     |      True
     |      
     |      However, we should point out that when we sliced the SeqRecord,
     |      any annotations dictionary or dbxrefs list entries were lost.
     |      You can explicitly copy them like this:
     |      
     |      >>> new.annotations = plasmid.annotations.copy()
     |      >>> new.dbxrefs = plasmid.dbxrefs[:]
     |  
     |  __bool__(self)
     |      Boolean value of an instance of this class (True).
     |      
     |      This behaviour is for backwards compatibility, since until the
     |      __len__ method was added, a SeqRecord always evaluated as True.
     |      
     |      Note that in comparison, a Seq object will evaluate to False if it
     |      has a zero length sequence.
     |      
     |      WARNING: The SeqRecord may in future evaluate to False when its
     |      sequence is of zero length (in order to better match the Seq
     |      object behaviour)!
     |  
     |  __bytes__(self)
     |  
     |  __contains__(self, char)
     |      Implement the 'in' keyword, searches the sequence.
     |      
     |      e.g.
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Fasta/sweetpea.nu", "fasta")
     |      >>> "GAATTC" in record
     |      False
     |      >>> "AAA" in record
     |      True
     |      
     |      This essentially acts as a proxy for using "in" on the sequence:
     |      
     |      >>> "GAATTC" in record.seq
     |      False
     |      >>> "AAA" in record.seq
     |      True
     |      
     |      Note that you can also use Seq objects as the query,
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> Seq("AAA") in record
     |      True
     |      
     |      See also the Seq object's __contains__ method.
     |  
     |  __eq__(self, other)
     |      Define the equal-to operand (not implemented).
     |  
     |  __format__(self, format_spec)
     |      Return the record as a string in the specified file format.
     |      
     |      This method supports the Python format() function and f-strings.
     |      The format_spec should be a lower case string supported by
     |      Bio.SeqIO as a text output file format. Requesting a binary file
     |      format raises a ValueError. e.g.
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
     |      ...                    id="YP_025292.1", name="HokC",
     |      ...                    description="toxic membrane protein")
     |      ...
     |      >>> format(record, "fasta")
     |      '>YP_025292.1 toxic membrane protein\nMKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF\n'
     |      >>> print(f"Here is {record.id} in FASTA format:\n{record:fasta}")
     |      Here is YP_025292.1 in FASTA format:
     |      >YP_025292.1 toxic membrane protein
     |      MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
     |      <BLANKLINE>
     |      
     |      See also the SeqRecord's format() method.
     |  
     |  __ge__(self, other)
     |      Define the greater-than-or-equal-to operand (not implemented).
     |  
     |  __getitem__(self, index)
     |      Return a sub-sequence or an individual letter.
     |      
     |      Slicing, e.g. my_record[5:10], returns a new SeqRecord for
     |      that sub-sequence with some annotation preserved as follows:
     |      
     |      * The name, id and description are kept as-is.
     |      * Any per-letter-annotations are sliced to match the requested
     |        sub-sequence.
     |      * Unless a stride is used, all those features which fall fully
     |        within the subsequence are included (with their locations
     |        adjusted accordingly). If you want to preserve any truncated
     |        features (e.g. GenBank/EMBL source features), you must
     |        explicitly add them to the new SeqRecord yourself.
     |      * With the exception of any molecule type, the annotations
     |        dictionary and the dbxrefs list are not used for the new
     |        SeqRecord, as in general they may not apply to the
     |        subsequence. If you want to preserve them, you must explicitly
     |        copy them to the new SeqRecord yourself.
     |      
     |      Using an integer index, e.g. my_record[5] is shorthand for
     |      extracting that letter from the sequence, my_record.seq[5].
     |      
     |      For example, consider this short protein and its secondary
     |      structure as encoded by the PDB (e.g. H for alpha helices),
     |      plus a simple feature for its histidine self phosphorylation
     |      site:
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> from Bio.SeqFeature import SeqFeature, SimpleLocation
     |      >>> rec = SeqRecord(Seq("MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLAT"
     |      ...                     "EMMSEQDGYLAESINKDIEECNAIIEQFIDYLR"),
     |      ...                 id="1JOY", name="EnvZ",
     |      ...                 description="Homodimeric domain of EnvZ from E. coli")
     |      >>> rec.letter_annotations["secondary_structure"] = "  S  SSSSSSHHHHHTTTHHHHHHHHHHHHHHHHHHHHHHTHHHHHHHHHHHHHHHHHHHHHTT  "
     |      >>> rec.features.append(SeqFeature(SimpleLocation(20, 21),
     |      ...                     type = "Site"))
     |      
     |      Now let's have a quick look at the full record,
     |      
     |      >>> print(rec)
     |      ID: 1JOY
     |      Name: EnvZ
     |      Description: Homodimeric domain of EnvZ from E. coli
     |      Number of features: 1
     |      Per letter annotation for: secondary_structure
     |      Seq('MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLATEMMSEQDGYLAESINKDIEE...YLR')
     |      >>> rec.letter_annotations["secondary_structure"]
     |      '  S  SSSSSSHHHHHTTTHHHHHHHHHHHHHHHHHHHHHHTHHHHHHHHHHHHHHHHHHHHHTT  '
     |      >>> print(rec.features[0].location)
     |      [20:21]
     |      
     |      Now let's take a sub sequence, here chosen as the first (fractured)
     |      alpha helix which includes the histidine phosphorylation site:
     |      
     |      >>> sub = rec[11:41]
     |      >>> print(sub)
     |      ID: 1JOY
     |      Name: EnvZ
     |      Description: Homodimeric domain of EnvZ from E. coli
     |      Number of features: 1
     |      Per letter annotation for: secondary_structure
     |      Seq('RTLLMAGVSHDLRTPLTRIRLATEMMSEQD')
     |      >>> sub.letter_annotations["secondary_structure"]
     |      'HHHHHTTTHHHHHHHHHHHHHHHHHHHHHH'
     |      >>> print(sub.features[0].location)
     |      [9:10]
     |      
     |      You can also of course omit the start or end values, for
     |      example to get the first ten letters only:
     |      
     |      >>> print(rec[:10])
     |      ID: 1JOY
     |      Name: EnvZ
     |      Description: Homodimeric domain of EnvZ from E. coli
     |      Number of features: 0
     |      Per letter annotation for: secondary_structure
     |      Seq('MAAGVKQLAD')
     |      
     |      Or for the last ten letters:
     |      
     |      >>> print(rec[-10:])
     |      ID: 1JOY
     |      Name: EnvZ
     |      Description: Homodimeric domain of EnvZ from E. coli
     |      Number of features: 0
     |      Per letter annotation for: secondary_structure
     |      Seq('IIEQFIDYLR')
     |      
     |      If you omit both, then you get a copy of the original record (although
     |      lacking the annotations and dbxrefs):
     |      
     |      >>> print(rec[:])
     |      ID: 1JOY
     |      Name: EnvZ
     |      Description: Homodimeric domain of EnvZ from E. coli
     |      Number of features: 1
     |      Per letter annotation for: secondary_structure
     |      Seq('MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLATEMMSEQDGYLAESINKDIEE...YLR')
     |      
     |      Finally, indexing with a simple integer is shorthand for pulling out
     |      that letter from the sequence directly:
     |      
     |      >>> rec[5]
     |      'K'
     |      >>> rec.seq[5]
     |      'K'
     |  
     |  __gt__(self, other)
     |      Define the greater-than operand (not implemented).
     |  
     |  __init__(self, seq, id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=None, features=None, annotations=None, letter_annotations=None)
     |      Create a SeqRecord.
     |      
     |      Arguments:
     |       - seq         - Sequence, required (Seq or MutableSeq)
     |       - id          - Sequence identifier, recommended (string)
     |       - name        - Sequence name, optional (string)
     |       - description - Sequence description, optional (string)
     |       - dbxrefs     - Database cross references, optional (list of strings)
     |       - features    - Any (sub)features, optional (list of SeqFeature objects)
     |       - annotations - Dictionary of annotations for the whole sequence
     |       - letter_annotations - Dictionary of per-letter-annotations, values
     |         should be strings, list or tuples of the same length as the full
     |         sequence.
     |      
     |      You will typically use Bio.SeqIO to read in sequences from files as
     |      SeqRecord objects.  However, you may want to create your own SeqRecord
     |      objects directly.
     |      
     |      Note that while an id is optional, we strongly recommend you supply a
     |      unique id string for each record.  This is especially important
     |      if you wish to write your sequences to a file.
     |      
     |      You can create a 'blank' SeqRecord object, and then populate the
     |      attributes later.
     |  
     |  __iter__(self)
     |      Iterate over the letters in the sequence.
     |      
     |      For example, using Bio.SeqIO to read in a protein FASTA file:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Fasta/loveliesbleeding.pro", "fasta")
     |      >>> for amino in record:
     |      ...     print(amino)
     |      ...     if amino == "L": break
     |      X
     |      A
     |      G
     |      L
     |      >>> print(record.seq[3])
     |      L
     |      
     |      This is just a shortcut for iterating over the sequence directly:
     |      
     |      >>> for amino in record.seq:
     |      ...     print(amino)
     |      ...     if amino == "L": break
     |      X
     |      A
     |      G
     |      L
     |      >>> print(record.seq[3])
     |      L
     |      
     |      Note that this does not facilitate iteration together with any
     |      per-letter-annotation.  However, you can achieve that using the
     |      python zip function on the record (or its sequence) and the relevant
     |      per-letter-annotation:
     |      
     |      >>> from Bio import SeqIO
     |      >>> rec = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
     |      >>> print("%s %s" % (rec.id, rec.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(rec.letter_annotations))
     |      ['solexa_quality']
     |      >>> for nuc, qual in zip(rec, rec.letter_annotations["solexa_quality"]):
     |      ...     if qual > 35:
     |      ...         print("%s %i" % (nuc, qual))
     |      A 40
     |      C 39
     |      G 38
     |      T 37
     |      A 36
     |      
     |      You may agree that using zip(rec.seq, ...) is more explicit than using
     |      zip(rec, ...) as shown above.
     |  
     |  __le__(self, other)
     |      Define the less-than-or-equal-to operand (not implemented).
     |  
     |  __len__(self)
     |      Return the length of the sequence.
     |      
     |      For example, using Bio.SeqIO to read in a FASTA nucleotide file:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Fasta/sweetpea.nu", "fasta")
     |      >>> len(record)
     |      309
     |      >>> len(record.seq)
     |      309
     |  
     |  __lt__(self, other)
     |      Define the less-than operand (not implemented).
     |  
     |  __ne__(self, other)
     |      Define the not-equal-to operand (not implemented).
     |  
     |  __radd__(self, other)
     |      Add another sequence or string to this sequence (from the left).
     |      
     |      This method handles adding a Seq object (or similar, e.g. MutableSeq)
     |      or a plain Python string (on the left) to a SeqRecord (on the right).
     |      See the __add__ method for more details, but for example:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
     |      >>> print("%s %s" % (record.id, record.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(record.letter_annotations))
     |      ['solexa_quality']
     |      
     |      >>> new = "ACT" + record
     |      >>> print("%s %s" % (new.id, new.seq))
     |      slxa_0001_1_0001_01 ACTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(new.letter_annotations))
     |      []
     |  
     |  __repr__(self)
     |      Return a concise summary of the record for debugging (string).
     |      
     |      The python built in function repr works by calling the object's __repr__
     |      method.  e.g.
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> rec = SeqRecord(Seq("MASRGVNKVILVGNLGQDPEVRYMPNGGAVANITLATSESWRDKAT"
     |      ...                     "GEMKEQTEWHRVVLFGKLAEVASEYLRKGSQVYIEGQLRTRKWTDQ"
     |      ...                     "SGQDRYTTEVVVNVGGTMQMLGGRQGGGAPAGGNIGGGQPQGGWGQ"
     |      ...                     "PQQPQGGNQFSGGAQSRPQQSAPAAPSNEPPMDFDDDIPF"),
     |      ...                 id="NP_418483.1", name="b4059",
     |      ...                 description="ssDNA-binding protein",
     |      ...                 dbxrefs=["ASAP:13298", "GI:16131885", "GeneID:948570"])
     |      >>> print(repr(rec))
     |      SeqRecord(seq=Seq('MASRGVNKVILVGNLGQDPEVRYMPNGGAVANITLATSESWRDKATGEMKEQTE...IPF'), id='NP_418483.1', name='b4059', description='ssDNA-binding protein', dbxrefs=['ASAP:13298', 'GI:16131885', 'GeneID:948570'])
     |      
     |      At the python prompt you can also use this shorthand:
     |      
     |      >>> rec
     |      SeqRecord(seq=Seq('MASRGVNKVILVGNLGQDPEVRYMPNGGAVANITLATSESWRDKATGEMKEQTE...IPF'), id='NP_418483.1', name='b4059', description='ssDNA-binding protein', dbxrefs=['ASAP:13298', 'GI:16131885', 'GeneID:948570'])
     |      
     |      Note that long sequences are shown truncated. Also note that any
     |      annotations, letter_annotations and features are not shown (as they
     |      would lead to a very long string).
     |  
     |  __str__(self)
     |      Return a human readable summary of the record and its annotation (string).
     |      
     |      The python built in function str works by calling the object's __str__
     |      method.  e.g.
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
     |      ...                    id="YP_025292.1", name="HokC",
     |      ...                    description="toxic membrane protein, small")
     |      >>> print(str(record))
     |      ID: YP_025292.1
     |      Name: HokC
     |      Description: toxic membrane protein, small
     |      Number of features: 0
     |      Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF')
     |      
     |      In this example you don't actually need to call str explicitly, as the
     |      print command does this automatically:
     |      
     |      >>> print(record)
     |      ID: YP_025292.1
     |      Name: HokC
     |      Description: toxic membrane protein, small
     |      Number of features: 0
     |      Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF')
     |      
     |      Note that long sequences are shown truncated.
     |  
     |  count(self, sub, start=None, end=None)
     |      Return the number of non-overlapping occurrences of sub in seq[start:end].
     |      
     |      Optional arguments start and end are interpreted as in slice notation.
     |      This method behaves as the count method of Python strings.
     |  
     |  format(self, format)
     |      Return the record as a string in the specified file format.
     |      
     |      The format should be a lower case string supported as an output
     |      format by Bio.SeqIO, which is used to turn the SeqRecord into a
     |      string.  e.g.
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
     |      ...                    id="YP_025292.1", name="HokC",
     |      ...                    description="toxic membrane protein")
     |      >>> record.format("fasta")
     |      '>YP_025292.1 toxic membrane protein\nMKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF\n'
     |      >>> print(record.format("fasta"))
     |      >YP_025292.1 toxic membrane protein
     |      MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
     |      <BLANKLINE>
     |      
     |      The Python print function automatically appends a new line, meaning
     |      in this example a blank line is shown.  If you look at the string
     |      representation you can see there is a trailing new line (shown as
     |      slash n) which is important when writing to a file or if
     |      concatenating multiple sequence strings together.
     |      
     |      Note that this method will NOT work on every possible file format
     |      supported by Bio.SeqIO (e.g. some are for multiple sequences only,
     |      and binary formats are not supported).
     |  
     |  islower(self)
     |      Return True if all ASCII characters in the record's sequence are lowercase.
     |      
     |      If there are no cased characters, the method returns False.
     |  
     |  isupper(self)
     |      Return True if all ASCII characters in the record's sequence are uppercase.
     |      
     |      If there are no cased characters, the method returns False.
     |  
     |  lower(self)
     |      Return a copy of the record with a lower case sequence.
     |      
     |      All the annotation is preserved unchanged. e.g.
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Fasta/aster.pro", "fasta")
     |      >>> print(record.format("fasta"))
     |      >gi|3298468|dbj|BAA31520.1| SAMIPF
     |      GGHVNPAVTFGAFVGGNITLLRGIVYIIAQLLGSTVACLLLKFVTNDMAVGVFSLSAGVG
     |      VTNALVFEIVMTFGLVYTVYATAIDPKKGSLGTIAPIAIGFIVGANI
     |      <BLANKLINE>
     |      >>> print(record.lower().format("fasta"))
     |      >gi|3298468|dbj|BAA31520.1| SAMIPF
     |      gghvnpavtfgafvggnitllrgivyiiaqllgstvaclllkfvtndmavgvfslsagvg
     |      vtnalvfeivmtfglvytvyataidpkkgslgtiapiaigfivgani
     |      <BLANKLINE>
     |      
     |      To take a more annotation rich example,
     |      
     |      >>> from Bio import SeqIO
     |      >>> old = SeqIO.read("EMBL/TRBG361.embl", "embl")
     |      >>> len(old.features)
     |      3
     |      >>> new = old.lower()
     |      >>> len(old.features) == len(new.features)
     |      True
     |      >>> old.annotations["organism"] == new.annotations["organism"]
     |      True
     |      >>> old.dbxrefs == new.dbxrefs
     |      True
     |  
     |  reverse_complement(self, id=False, name=False, description=False, features=True, annotations=False, letter_annotations=True, dbxrefs=False)
     |      Return new SeqRecord with reverse complement sequence.
     |      
     |      By default the new record does NOT preserve the sequence identifier,
     |      name, description, general annotation or database cross-references -
     |      these are unlikely to apply to the reversed sequence.
     |      
     |      You can specify the returned record's id, name and description as
     |      strings, or True to keep that of the parent, or False for a default.
     |      
     |      You can specify the returned record's features with a list of
     |      SeqFeature objects, or True to keep that of the parent, or False to
     |      omit them. The default is to keep the original features (with the
     |      strand and locations adjusted).
     |      
     |      You can also specify both the returned record's annotations and
     |      letter_annotations as dictionaries, True to keep that of the parent,
     |      or False to omit them. The default is to keep the original
     |      annotations (with the letter annotations reversed).
     |      
     |      To show what happens to the pre-letter annotations, consider an
     |      example Solexa variant FASTQ file with a single entry, which we'll
     |      read in as a SeqRecord:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
     |      >>> print("%s %s" % (record.id, record.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(record.letter_annotations))
     |      ['solexa_quality']
     |      >>> print(record.letter_annotations["solexa_quality"])
     |      [40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5]
     |      
     |      Now take the reverse complement, here we explicitly give a new
     |      identifier (the old identifier with a suffix):
     |      
     |      >>> rc_record = record.reverse_complement(id=record.id + "_rc")
     |      >>> print("%s %s" % (rc_record.id, rc_record.seq))
     |      slxa_0001_1_0001_01_rc NNNNNNACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
     |      
     |      Notice that the per-letter-annotations have also been reversed,
     |      although this may not be appropriate for all cases.
     |      
     |      >>> print(rc_record.letter_annotations["solexa_quality"])
     |      [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]
     |      
     |      Now for the features, we need a different example. Parsing a GenBank
     |      file is probably the easiest way to get an nice example with features
     |      in it...
     |      
     |      >>> from Bio import SeqIO
     |      >>> with open("GenBank/pBAD30.gb") as handle:
     |      ...     plasmid = SeqIO.read(handle, "gb")
     |      >>> print("%s %i" % (plasmid.id, len(plasmid)))
     |      pBAD30 4923
     |      >>> plasmid.seq
     |      Seq('GCTAGCGGAGTGTATACTGGCTTACTATGTTGGCACTGATGAGGGTGTCAGTGA...ATG')
     |      >>> len(plasmid.features)
     |      13
     |      
     |      Now, let's take the reverse complement of this whole plasmid:
     |      
     |      >>> rc_plasmid = plasmid.reverse_complement(id=plasmid.id+"_rc")
     |      >>> print("%s %i" % (rc_plasmid.id, len(rc_plasmid)))
     |      pBAD30_rc 4923
     |      >>> rc_plasmid.seq
     |      Seq('CATGGGCAAATATTATACGCAAGGCGACAAGGTGCTGATGCCGCTGGCGATTCA...AGC')
     |      >>> len(rc_plasmid.features)
     |      13
     |      
     |      Let's compare the first CDS feature - it has gone from being the
     |      second feature (index 1) to the second last feature (index -2), its
     |      strand has changed, and the location switched round.
     |      
     |      >>> print(plasmid.features[1])
     |      type: CDS
     |      location: [1081:1960](-)
     |      qualifiers:
     |          Key: label, Value: ['araC']
     |          Key: note, Value: ['araC regulator of the arabinose BAD promoter']
     |          Key: vntifkey, Value: ['4']
     |      <BLANKLINE>
     |      >>> print(rc_plasmid.features[-2])
     |      type: CDS
     |      location: [2963:3842](+)
     |      qualifiers:
     |          Key: label, Value: ['araC']
     |          Key: note, Value: ['araC regulator of the arabinose BAD promoter']
     |          Key: vntifkey, Value: ['4']
     |      <BLANKLINE>
     |      
     |      You can check this new location, based on the length of the plasmid:
     |      
     |      >>> len(plasmid) - 1081
     |      3842
     |      >>> len(plasmid) - 1960
     |      2963
     |      
     |      Note that if the SeqFeature annotation includes any strand specific
     |      information (e.g. base changes for a SNP), this information is not
     |      amended, and would need correction after the reverse complement.
     |      
     |      Note trying to reverse complement a protein SeqRecord raises an
     |      exception:
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> protein_rec = SeqRecord(Seq("MAIVMGR"), id="Test",
     |      ...                         annotations={"molecule_type": "protein"})
     |      >>> protein_rec.reverse_complement()
     |      Traceback (most recent call last):
     |         ...
     |      ValueError: Proteins do not have complements!
     |      
     |      If you have RNA without any U bases, it must be annotated as RNA
     |      otherwise it will be treated as DNA by default with A mapped to T:
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> rna1 = SeqRecord(Seq("ACG"), id="Test")
     |      >>> rna2 = SeqRecord(Seq("ACG"), id="Test", annotations={"molecule_type": "RNA"})
     |      >>> print(rna1.reverse_complement(id="RC", description="unk").format("fasta"))
     |      >RC unk
     |      CGT
     |      <BLANKLINE>
     |      >>> print(rna2.reverse_complement(id="RC", description="RNA").format("fasta"))
     |      >RC RNA
     |      CGU
     |      <BLANKLINE>
     |      
     |      Also note you can reverse complement a SeqRecord using a MutableSeq:
     |      
     |      >>> from Bio.Seq import MutableSeq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> rec = SeqRecord(MutableSeq("ACGT"), id="Test")
     |      >>> rec.seq[0] = "T"
     |      >>> print("%s %s" % (rec.id, rec.seq))
     |      Test TCGT
     |      >>> rc = rec.reverse_complement(id=True)
     |      >>> print("%s %s" % (rc.id, rc.seq))
     |      Test ACGA
     |  
     |  translate(self, table='Standard', stop_symbol='*', to_stop=False, cds=False, gap=None, id=False, name=False, description=False, features=False, annotations=False, letter_annotations=False, dbxrefs=False)
     |      Return new SeqRecord with translated sequence.
     |      
     |      This calls the record's .seq.translate() method (which describes
     |      the translation related arguments, like table for the genetic code),
     |      
     |      By default the new record does NOT preserve the sequence identifier,
     |      name, description, general annotation or database cross-references -
     |      these are unlikely to apply to the translated sequence.
     |      
     |      You can specify the returned record's id, name and description as
     |      strings, or True to keep that of the parent, or False for a default.
     |      
     |      You can specify the returned record's features with a list of
     |      SeqFeature objects, or False (default) to omit them.
     |      
     |      You can also specify both the returned record's annotations and
     |      letter_annotations as dictionaries, True to keep that of the parent
     |      (annotations only), or False (default) to omit them.
     |      
     |      e.g. Loading a FASTA gene and translating it,
     |      
     |      >>> from Bio import SeqIO
     |      >>> gene_record = SeqIO.read("Fasta/sweetpea.nu", "fasta")
     |      >>> print(gene_record.format("fasta"))
     |      >gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA) gene, partial cds
     |      CAGGCTGCGCGGTTTCTATTTATGAAGAACAAGGTCCGTATGATAGTTGATTGTCATGCA
     |      AAACATGTGAAGGTTCTTCAAGACGAAAAACTCCCATTTGATTTGACTCTGTGCGGTTCG
     |      ACCTTAAGAGCTCCACATAGTTGCCATTTGCAGTACATGGCTAACATGGATTCAATTGCT
     |      TCATTGGTTATGGCAGTGGTCGTCAATGACAGCGATGAAGATGGAGATAGCCGTGACGCA
     |      GTTCTACCACAAAAGAAAAAGAGACTTTGGGGTTTGGTAGTTTGTCATAACACTACTCCG
     |      AGGTTTGTT
     |      <BLANKLINE>
     |      
     |      And now translating the record, specifying the new ID and description:
     |      
     |      >>> protein_record = gene_record.translate(table=11,
     |      ...                                        id="phya",
     |      ...                                        description="translation")
     |      >>> print(protein_record.format("fasta"))
     |      >phya translation
     |      QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
     |      SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV
     |      <BLANKLINE>
     |  
     |  upper(self)
     |      Return a copy of the record with an upper case sequence.
     |      
     |      All the annotation is preserved unchanged. e.g.
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> record = SeqRecord(Seq("acgtACGT"), id="Test",
     |      ...                    description = "Made up for this example")
     |      >>> record.letter_annotations["phred_quality"] = [1, 2, 3, 4, 5, 6, 7, 8]
     |      >>> print(record.upper().format("fastq"))
     |      @Test Made up for this example
     |      ACGTACGT
     |      +
     |      "#$%&'()
     |      <BLANKLINE>
     |      
     |      Naturally, there is a matching lower method:
     |      
     |      >>> print(record.lower().format("fastq"))
     |      @Test Made up for this example
     |      acgtacgt
     |      +
     |      "#$%&'()
     |      <BLANKLINE>
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
     |  
     |  letter_annotations
     |      Dictionary of per-letter-annotation for the sequence.
     |      
     |      For example, this can hold quality scores used in FASTQ or QUAL files.
     |      Consider this example using Bio.SeqIO to read in an example Solexa
     |      variant FASTQ file as a SeqRecord:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
     |      >>> print("%s %s" % (record.id, record.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(record.letter_annotations))
     |      ['solexa_quality']
     |      >>> print(record.letter_annotations["solexa_quality"])
     |      [40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5]
     |      
     |      The letter_annotations get sliced automatically if you slice the
     |      parent SeqRecord, for example taking the last ten bases:
     |      
     |      >>> sub_record = record[-10:]
     |      >>> print("%s %s" % (sub_record.id, sub_record.seq))
     |      slxa_0001_1_0001_01 ACGTNNNNNN
     |      >>> print(sub_record.letter_annotations["solexa_quality"])
     |      [4, 3, 2, 1, 0, -1, -2, -3, -4, -5]
     |      
     |      Any python sequence (i.e. list, tuple or string) can be recorded in
     |      the SeqRecord's letter_annotations dictionary as long as the length
     |      matches that of the SeqRecord's sequence.  e.g.
     |      
     |      >>> len(sub_record.letter_annotations)
     |      1
     |      >>> sub_record.letter_annotations["dummy"] = "abcdefghij"
     |      >>> len(sub_record.letter_annotations)
     |      2
     |      
     |      You can delete entries from the letter_annotations dictionary as usual:
     |      
     |      >>> del sub_record.letter_annotations["solexa_quality"]
     |      >>> sub_record.letter_annotations
     |      {'dummy': 'abcdefghij'}
     |      
     |      You can completely clear the dictionary easily as follows:
     |      
     |      >>> sub_record.letter_annotations = {}
     |      >>> sub_record.letter_annotations
     |      {}
     |      
     |      Note that if replacing the record's sequence with a sequence of a
     |      different length you must first clear the letter_annotations dict.
     |  
     |  seq
     |      The sequence itself, as a Seq or MutableSeq object.
     |  
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |  
     |  __hash__ = None
    



```python
from Bio.Seq import Seq
```


```python
simple_seq = Seq("GATC")
```


```python
simple_seq_r = SeqRecord(simple_seq)
```


```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
# we can also pass an ID to it
simple_seq_r.id = "AC12345"
```


```python
simple_seq_r.description = "Made up sequence for the VDB computational biology class"
```


```python
print(simple_seq_r.description)
```

    Made up sequence for the VDB computational biology class



```python
simple_seq_r.seq
```




    Seq('GATC')




```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB computational biology class', dbxrefs=[])




```python
simple_seq_r.annotations["evidence"] = "None. this is just an example"
```


```python
print(simple_seq_r.annotations["evidence"])
```

    None. this is just an example



```python
# you can also make per letter annotations:
simple_seq_r.letter_annotations["phred_quality"] = [40, 40, 38, 30]
```


```python
print(simple_seq_r.letter_annotations)
```

    {'phred_quality': [40, 40, 38, 30]}



```python
from Bio import SeqIO
```


```python
#https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.gb
```


```python
record = SeqIO.read("NC_005816.fna.txt", "fasta")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='gi|45478711|ref|NC_005816.1|', name='gi|45478711|ref|NC_005816.1|', description='gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
# we can individually pull the ID
record.id
```




    'gi|45478711|ref|NC_005816.1|'




```python
# and the description
record.description
```




    'gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
record.dbxrefs
```




    []




```python
record.annotations
```




    {}




```python
record.features
```




    []




```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.id
```




    'NC_005816.1'




```python
record.name
```




    'NC_005816'




```python
record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
record.letter_annotations
```




    {}




```python
len(record.annotations)
```




    13




```python
record.annotations["source"]
```




    'Yersinia pestis biovar Microtus str. 91001'




```python
record.dbxrefs
```




    ['Project:58037']




```python
len(record.features)
```




    41




```python
from Bio import SeqFeature
```


```python
# we do this because we believe that the start position for this genetic feature starts after position 5
start_pos = SeqFeature.AfterPosition(5)
```


```python
end_pos = SeqFeature.BetweenPosition(9, left = 8, right = 9)
```


```python
my_location = SeqFeature.SimpleLocation(start_pos, end_pos)
```


```python
#now that we have our start and end location, we can print our location
print(my_location)
```

    [>5:(8^9)]



```python
my_location.start
```




    AfterPosition(5)




```python
my_location.end
```




    BetweenPosition(9, left=8, right=9)




```python
# we can also do this if we dont want a description for our location
int(my_location.start)
```




    5




```python
exact_location = SeqFeature.SimpleLocation(5,9)
```


```python
print(exact_location)
```

    [5:9]



```python
exact_location.start
```




    ExactPosition(5)

## Sequence Input/ Output
```python
from Bio import SeqIO
```


```python
# this will print the sequence as a fasta file
for seq_record in SeqIO.parse("orchid.fasta.txt", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```

    gi|2765658|emb|Z78533.1|CIZ78533
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    gi|2765657|emb|Z78532.1|CCZ78532
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')
    753
    gi|2765656|emb|Z78531.1|CFZ78531
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')
    748
    gi|2765655|emb|Z78530.1|CMZ78530
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')
    744
    gi|2765654|emb|Z78529.1|CLZ78529
    Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')
    733
    gi|2765652|emb|Z78527.1|CYZ78527
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')
    718
    gi|2765651|emb|Z78526.1|CGZ78526
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730
    gi|2765650|emb|Z78525.1|CAZ78525
    Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')
    704
    gi|2765649|emb|Z78524.1|CFZ78524
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')
    740
    gi|2765648|emb|Z78523.1|CHZ78523
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
    709
    gi|2765647|emb|Z78522.1|CMZ78522
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')
    700
    gi|2765646|emb|Z78521.1|CCZ78521
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')
    726
    gi|2765645|emb|Z78520.1|CSZ78520
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')
    753
    gi|2765644|emb|Z78519.1|CPZ78519
    Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')
    699
    gi|2765643|emb|Z78518.1|CRZ78518
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')
    658
    gi|2765642|emb|Z78517.1|CFZ78517
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')
    752
    gi|2765641|emb|Z78516.1|CPZ78516
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')
    726
    gi|2765640|emb|Z78515.1|MXZ78515
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')
    765
    gi|2765639|emb|Z78514.1|PSZ78514
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')
    755
    gi|2765638|emb|Z78513.1|PBZ78513
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')
    742
    gi|2765637|emb|Z78512.1|PWZ78512
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')
    762
    gi|2765636|emb|Z78511.1|PEZ78511
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')
    745
    gi|2765635|emb|Z78510.1|PCZ78510
    Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')
    750
    gi|2765634|emb|Z78509.1|PPZ78509
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')
    731
    gi|2765633|emb|Z78508.1|PLZ78508
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')
    741
    gi|2765632|emb|Z78507.1|PLZ78507
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')
    740
    gi|2765631|emb|Z78506.1|PLZ78506
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')
    727
    gi|2765630|emb|Z78505.1|PSZ78505
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')
    711
    gi|2765629|emb|Z78504.1|PKZ78504
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')
    743
    gi|2765628|emb|Z78503.1|PCZ78503
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')
    727
    gi|2765627|emb|Z78502.1|PBZ78502
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')
    757
    gi|2765626|emb|Z78501.1|PCZ78501
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')
    770
    gi|2765625|emb|Z78500.1|PWZ78500
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')
    767
    gi|2765624|emb|Z78499.1|PMZ78499
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')
    759
    gi|2765623|emb|Z78498.1|PMZ78498
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')
    750
    gi|2765622|emb|Z78497.1|PDZ78497
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    788
    gi|2765621|emb|Z78496.1|PAZ78496
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    774
    gi|2765620|emb|Z78495.1|PEZ78495
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')
    789
    gi|2765619|emb|Z78494.1|PNZ78494
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')
    688
    gi|2765618|emb|Z78493.1|PGZ78493
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')
    719
    gi|2765617|emb|Z78492.1|PBZ78492
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')
    743
    gi|2765616|emb|Z78491.1|PCZ78491
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    737
    gi|2765615|emb|Z78490.1|PFZ78490
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    728
    gi|2765614|emb|Z78489.1|PDZ78489
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')
    740
    gi|2765613|emb|Z78488.1|PTZ78488
    Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')
    696
    gi|2765612|emb|Z78487.1|PHZ78487
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    gi|2765611|emb|Z78486.1|PBZ78486
    Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')
    731
    gi|2765610|emb|Z78485.1|PHZ78485
    Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')
    735
    gi|2765609|emb|Z78484.1|PCZ78484
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')
    720
    gi|2765608|emb|Z78483.1|PVZ78483
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    740
    gi|2765607|emb|Z78482.1|PEZ78482
    Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')
    629
    gi|2765606|emb|Z78481.1|PIZ78481
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    572
    gi|2765605|emb|Z78480.1|PGZ78480
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    587
    gi|2765604|emb|Z78479.1|PPZ78479
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')
    700
    gi|2765603|emb|Z78478.1|PVZ78478
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')
    636
    gi|2765602|emb|Z78477.1|PVZ78477
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    716
    gi|2765601|emb|Z78476.1|PGZ78476
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    592
    gi|2765600|emb|Z78475.1|PSZ78475
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    716
    gi|2765599|emb|Z78474.1|PKZ78474
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')
    733
    gi|2765598|emb|Z78473.1|PSZ78473
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    626
    gi|2765597|emb|Z78472.1|PLZ78472
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    737
    gi|2765596|emb|Z78471.1|PDZ78471
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    gi|2765595|emb|Z78470.1|PPZ78470
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    574
    gi|2765594|emb|Z78469.1|PHZ78469
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    594
    gi|2765593|emb|Z78468.1|PAZ78468
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')
    610
    gi|2765592|emb|Z78467.1|PSZ78467
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    730
    gi|2765591|emb|Z78466.1|PPZ78466
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    641
    gi|2765590|emb|Z78465.1|PRZ78465
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    702
    gi|2765589|emb|Z78464.1|PGZ78464
    Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')
    733
    gi|2765588|emb|Z78463.1|PGZ78463
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')
    738
    gi|2765587|emb|Z78462.1|PSZ78462
    Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')
    736
    gi|2765586|emb|Z78461.1|PWZ78461
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    gi|2765585|emb|Z78460.1|PCZ78460
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')
    745
    gi|2765584|emb|Z78459.1|PDZ78459
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')
    744
    gi|2765583|emb|Z78458.1|PHZ78458
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')
    738
    gi|2765582|emb|Z78457.1|PCZ78457
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')
    739
    gi|2765581|emb|Z78456.1|PTZ78456
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    gi|2765580|emb|Z78455.1|PJZ78455
    Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')
    745
    gi|2765579|emb|Z78454.1|PFZ78454
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')
    695
    gi|2765578|emb|Z78453.1|PSZ78453
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    745
    gi|2765577|emb|Z78452.1|PBZ78452
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    743
    gi|2765576|emb|Z78451.1|PHZ78451
    Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')
    730
    gi|2765575|emb|Z78450.1|PPZ78450
    Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')
    706
    gi|2765574|emb|Z78449.1|PMZ78449
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    744
    gi|2765573|emb|Z78448.1|PAZ78448
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    742
    gi|2765572|emb|Z78447.1|PVZ78447
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')
    694
    gi|2765571|emb|Z78446.1|PAZ78446
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')
    712
    gi|2765570|emb|Z78445.1|PUZ78445
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')
    715
    gi|2765569|emb|Z78444.1|PAZ78444
    Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')
    688
    gi|2765568|emb|Z78443.1|PLZ78443
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    784
    gi|2765567|emb|Z78442.1|PBZ78442
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')
    721
    gi|2765566|emb|Z78441.1|PSZ78441
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703
    gi|2765565|emb|Z78440.1|PPZ78440
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')
    744
    gi|2765564|emb|Z78439.1|PBZ78439
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
# this will print a cleaner version than the fasta
for seq_record in SeqIO.parse("orchid.gbk.txt", "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```

    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    Z78532.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')
    753
    Z78531.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')
    748
    Z78530.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')
    744
    Z78529.1
    Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')
    733
    Z78527.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')
    718
    Z78526.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730
    Z78525.1
    Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')
    704
    Z78524.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')
    740
    Z78523.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
    709
    Z78522.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')
    700
    Z78521.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')
    726
    Z78520.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')
    753
    Z78519.1
    Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')
    699
    Z78518.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')
    658
    Z78517.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')
    752
    Z78516.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')
    726
    Z78515.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')
    765
    Z78514.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')
    755
    Z78513.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')
    742
    Z78512.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')
    762
    Z78511.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')
    745
    Z78510.1
    Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')
    750
    Z78509.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')
    731
    Z78508.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')
    741
    Z78507.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')
    740
    Z78506.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')
    727
    Z78505.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')
    711
    Z78504.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')
    743
    Z78503.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')
    727
    Z78502.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')
    757
    Z78501.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')
    770
    Z78500.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')
    767
    Z78499.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')
    759
    Z78498.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')
    750
    Z78497.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    788
    Z78496.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    774
    Z78495.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')
    789
    Z78494.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')
    688
    Z78493.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')
    719
    Z78492.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')
    743
    Z78491.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    737
    Z78490.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    728
    Z78489.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')
    740
    Z78488.1
    Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')
    696
    Z78487.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78486.1
    Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')
    731
    Z78485.1
    Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')
    735
    Z78484.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')
    720
    Z78483.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    740
    Z78482.1
    Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')
    629
    Z78481.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    572
    Z78480.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    587
    Z78479.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')
    700
    Z78478.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')
    636
    Z78477.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    716
    Z78476.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    592
    Z78475.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    716
    Z78474.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')
    733
    Z78473.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    626
    Z78472.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    737
    Z78471.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78470.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    574
    Z78469.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    594
    Z78468.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')
    610
    Z78467.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    730
    Z78466.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    641
    Z78465.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    702
    Z78464.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')
    733
    Z78463.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')
    738
    Z78462.1
    Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')
    736
    Z78461.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78460.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')
    745
    Z78459.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')
    744
    Z78458.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')
    738
    Z78457.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')
    739
    Z78456.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78455.1
    Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')
    745
    Z78454.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')
    695
    Z78453.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    745
    Z78452.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    743
    Z78451.1
    Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')
    730
    Z78450.1
    Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')
    706
    Z78449.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    744
    Z78448.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    742
    Z78447.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')
    694
    Z78446.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')
    712
    Z78445.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')
    715
    Z78444.1
    Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')
    688
    Z78443.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    784
    Z78442.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')
    721
    Z78441.1
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703
    Z78440.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')
    744
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
# this will show the identifiers
identifiers = [seq_record.id for seq_record in SeqIO.parse("orchid.gbk.txt", "genbank")]
```


```python
identifiers
```




    ['Z78533.1',
     'Z78532.1',
     'Z78531.1',
     'Z78530.1',
     'Z78529.1',
     'Z78527.1',
     'Z78526.1',
     'Z78525.1',
     'Z78524.1',
     'Z78523.1',
     'Z78522.1',
     'Z78521.1',
     'Z78520.1',
     'Z78519.1',
     'Z78518.1',
     'Z78517.1',
     'Z78516.1',
     'Z78515.1',
     'Z78514.1',
     'Z78513.1',
     'Z78512.1',
     'Z78511.1',
     'Z78510.1',
     'Z78509.1',
     'Z78508.1',
     'Z78507.1',
     'Z78506.1',
     'Z78505.1',
     'Z78504.1',
     'Z78503.1',
     'Z78502.1',
     'Z78501.1',
     'Z78500.1',
     'Z78499.1',
     'Z78498.1',
     'Z78497.1',
     'Z78496.1',
     'Z78495.1',
     'Z78494.1',
     'Z78493.1',
     'Z78492.1',
     'Z78491.1',
     'Z78490.1',
     'Z78489.1',
     'Z78488.1',
     'Z78487.1',
     'Z78486.1',
     'Z78485.1',
     'Z78484.1',
     'Z78483.1',
     'Z78482.1',
     'Z78481.1',
     'Z78480.1',
     'Z78479.1',
     'Z78478.1',
     'Z78477.1',
     'Z78476.1',
     'Z78475.1',
     'Z78474.1',
     'Z78473.1',
     'Z78472.1',
     'Z78471.1',
     'Z78470.1',
     'Z78469.1',
     'Z78468.1',
     'Z78467.1',
     'Z78466.1',
     'Z78465.1',
     'Z78464.1',
     'Z78463.1',
     'Z78462.1',
     'Z78461.1',
     'Z78460.1',
     'Z78459.1',
     'Z78458.1',
     'Z78457.1',
     'Z78456.1',
     'Z78455.1',
     'Z78454.1',
     'Z78453.1',
     'Z78452.1',
     'Z78451.1',
     'Z78450.1',
     'Z78449.1',
     'Z78448.1',
     'Z78447.1',
     'Z78446.1',
     'Z78445.1',
     'Z78444.1',
     'Z78443.1',
     'Z78442.1',
     'Z78441.1',
     'Z78440.1',
     'Z78439.1']




```python
record_iterator = SeqIO.parse("orchid.fasta.txt", "fasta")
```


```python
first_record = next(record_iterator)
```


```python
print(first_record.id)
```

    gi|2765658|emb|Z78533.1|CIZ78533



```python
#we first printed the ID and now we are printing the description
print(first_record.description)
```

    gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
second_record = next(record_iterator)
```


```python
print(second_record.description)
```

    gi|2765657|emb|Z78532.1|CCZ78532 C.californicum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
records = list(SeqIO.parse("orchid.gbk.txt", "genbank"))
```


```python
print("found %i records" % len(records))
```

    found 94 records



```python
print("the last record")
last_record = records[-1]
print(last_record.id)
print(repr(last_record.seq))
print(len(last_record))
```

    the last record
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
print("the first record")
first_record = records[0]
print(first_record.id)
print(repr(first_record.seq))
print(len(first_record))
```

    the first record
    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740



```python
from Bio import SeqIO
```


```python
record_iterator = SeqIO.parse("orchid.gbk.txt", "genbank")
```


```python
first_record = next(record_iterator)
```


```python
# this will print the rna gene where we can read it in an organized fashion
print(first_record)
```

    ID: Z78533.1
    Name: Z78533
    Description: C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA
    Number of features: 5
    /molecule_type=DNA
    /topology=linear
    /data_file_division=PLN
    /date=30-NOV-2006
    /accessions=['Z78533']
    /sequence_version=1
    /gi=2765658
    /keywords=['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2']
    /source=Cypripedium irapeanum
    /organism=Cypripedium irapeanum
    /taxonomy=['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium']
    /references=[Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')



```python
print(first_record.annotations)
```

    {'molecule_type': 'DNA', 'topology': 'linear', 'data_file_division': 'PLN', 'date': '30-NOV-2006', 'accessions': ['Z78533'], 'sequence_version': 1, 'gi': '2765658', 'keywords': ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'source': 'Cypripedium irapeanum', 'organism': 'Cypripedium irapeanum', 'taxonomy': ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], 'references': [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]}



```python
# this gives us the keys/ instead of all of the answers
print(first_record.annotations.keys())
```

    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references'])



```python
print(first_record.annotations.values())
```

    dict_values(['DNA', 'linear', 'PLN', '30-NOV-2006', ['Z78533'], 1, '2765658', ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'Cypripedium irapeanum', 'Cypripedium irapeanum', ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]])



```python
print(first_record.annotations["source"])
```

    Cypripedium irapeanum



```python
print(first_record.annotations["organism"])
```

    Cypripedium irapeanum



```python
all_species = []
```


```python
for seq_record in SeqIO.parse("orchid.gbk.txt", "genbank"):
    all_species.append(seq_record.annotations["organism"])
```


```python
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
all_species = [
    seq_record.annotations["organism"]
    for seq_record in SeqIO.parse("orchid.gbk.txt", "genbank")
]
```


```python
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
all_species = []
```


```python
for seq_record in SeqIO.parse("orchid.fasta.txt", "fasta"):
    all_species.append(seq_record.description.split()[1])
```


```python
# because the species are named differently in the fasta file, they are abbreviated here
print(all_species)
```

    ['C.irapeanum', 'C.californicum', 'C.fasciculatum', 'C.margaritaceum', 'C.lichiangense', 'C.yatabeanum', 'C.guttatum', 'C.acaule', 'C.formosanum', 'C.himalaicum', 'C.macranthum', 'C.calceolus', 'C.segawai', 'C.pubescens', 'C.reginae', 'C.flavum', 'C.passerinum', 'M.xerophyticum', 'P.schlimii', 'P.besseae', 'P.wallisii', 'P.exstaminodium', 'P.caricinum', 'P.pearcei', 'P.longifolium', 'P.lindenii', 'P.lindleyanum', 'P.sargentianum', 'P.kaiteurum', 'P.czerwiakowianum', 'P.boissierianum', 'P.caudatum', 'P.warszewiczianum', 'P.micranthum', 'P.malipoense', 'P.delenatii', 'P.armeniacum', 'P.emersonii', 'P.niveum', 'P.godefroyae', 'P.bellatulum', 'P.concolor', 'P.fairrieanum', 'P.druryi', 'P.tigrinum', 'P.hirsutissimum', 'P.barbigerum', 'P.henryanum', 'P.charlesworthii', 'P.villosum', 'P.exul', 'P.insigne', 'P.gratrixianum', 'P.primulinum', 'P.victoria', 'P.victoria', 'P.glaucophyllum', 'P.supardii', 'P.kolopakingii', 'P.sanderianum', 'P.lowii', 'P.dianthum', 'P.parishii', 'P.haynaldianum', 'P.adductum', 'P.stonei', 'P.philippinense', 'P.rothschildianum', 'P.glanduliferum', 'P.glanduliferum', 'P.sukhakulii', 'P.wardii', 'P.ciliolare', 'P.dayanum', 'P.hennisianum', 'P.callosum', 'P.tonsum', 'P.javanicum', 'P.fowliei', 'P.schoseri', 'P.bougainvilleanum', 'P.hookerae', 'P.papuanum', 'P.mastersianum', 'P.argus', 'P.venustum', 'P.acmodontum', 'P.urbanianum', 'P.appletonianum', 'P.lawrenceanum', 'P.bullenianum', 'P.superbiens', 'P.purpuratum', 'P.barbatum']



```python
record_iterator = SeqIO.parse("orchid.fasta.txt", "fasta")
```


```python
first_record = next(record_iterator)
```


```python
first_record.id
```




    'gi|2765658|emb|Z78533.1|CIZ78533'




```python
first_record.id = "new_id"
```


```python
first_record.id
```




    'new_id'




```python
first_record
```




    SeqRecord(seq=Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC'), id='new_id', name='gi|2765658|emb|Z78533.1|CIZ78533', description='gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA', dbxrefs=[])




```python
first_record.description = first_record.id + " " + "mutations induced randomly"
```


```python
print(first_record.format("fasta"[:200]))
```

    >new_id mutations induced randomly
    CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGGAATAAA
    CGATCGAGTGAATCCGGAGGACCGGTGTACTCAGCTCACCGGGGGCATTGCTCCCGTGGT
    GACCCTGATTTGTTGTTGGGCCGCCTCGGGAGCGTCCATGGCGGGTTTGAACCTCTAGCC
    CGGCGCAGTTTGGGCGCCAAGCCATATGAAAGCATCACCGGCGAATGGCATTGTCTTCCC
    CAAAACCCGGAGCGGCGGCGTGCTGTCGCGTGCCCAATGAATTTTGATGACTCTCGCAAA
    CGGGAATCTTGGCTCTTTGCATCGGATGGAAGGACGCAGCGAAATGCGATAAGTGGTGTG
    AATTGCAAGATCCCGTGAACCATCGAGTCTTTTGAACGCAAGTTGCGCCCGAGGCCATCA
    GGCTAAGGGCACGCCTGCTTGGGCGTCGCGCTTCGTCTCTCTCCTGCCAATGCTTGCCCG
    GCATACAGCCAGGCCGGCGTGGTGCGGATGTGAAAGATTGGCCCCTTGTGCCTAGGTGCG
    GCGGGTCCAAGAGCTGGTGTTTTGATGGCCCGGAACCCGGCAAGAGGTGGACGGATGCTG
    GCAGCAGCTGCCGTGCGAATCCCCCATGTTGTCGTGCTTGTCGGACAGGCAGGAGAACCC
    TTCCGAACCCCAATGGAGGGCGGTTGACCGCCATTCGGATGTGACCCCAGGTCAGGCGGG
    GGCACCCGCTGAGTTTACGC
    



```python
from Bio.Seq import Seq
```


```python
from Bio.SeqRecord import SeqRecord
```


```python
rec1 = SeqRecord(
Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD"
   "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
   "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
   "SSAC"
),
    id = "gi|14150838|gb|AAK54658.1|AF376133_1",
    description = "chalcone synthase [Cucumis sativus]"
```
## Seqeunce Alignment

```python
# sequence alignments take multiple versions of a sequence and compare them to determine how closely they are related.
```


```python
#https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/PF05371_seed.sth
```


```python
# Import AlignIO
from Bio import AlignIO
```


```python
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
```


```python
# Show alignment
print(alignment)
```

    Alignment with 7 rows and 52 columns
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
```


```python
# Show alignment length
print("Alignment length %i" % alignment.get_alignment_length())
```

    Alignment length 52



```python
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
for record in alignment:
    if record.dbxrefs:
        print("%s %s" % (record.id, record.dbxrefs))
```

    COATB_BPIKE/30-81 ['PDB; 1ifl ; 1-52;']
    COATB_BPM13/24-72 ['PDB; 2cpb ; 1-49;', 'PDB; 2cps ; 1-49;']
    Q9T0Q9_BPFD/1-49 ['PDB; 1nh4 A; 1-49;']
    COATB_BPIF1/22-73 ['PDB; 1ifk ; 1-50;']



```python
# show all information for all alignments
for record in alignment:
    print(record)
```

    ID: COATB_BPIKE/30-81
    Name: COATB_BPIKE
    Description: COATB_BPIKE/30-81
    Database cross-references: PDB; 1ifl ; 1-52;
    Number of features: 0
    /accession=P03620.1
    /start=30
    /end=81
    Per letter annotation for: secondary_structure
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA')
    ID: Q9T0Q8_BPIKE/1-52
    Name: Q9T0Q8_BPIKE
    Description: Q9T0Q8_BPIKE/1-52
    Number of features: 0
    /accession=Q9T0Q8.1
    /start=1
    /end=52
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA')
    ID: COATB_BPI22/32-83
    Name: COATB_BPI22
    Description: COATB_BPI22/32-83
    Number of features: 0
    /accession=P15416.1
    /start=32
    /end=83
    Seq('DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA')
    ID: COATB_BPM13/24-72
    Name: COATB_BPM13
    Description: COATB_BPM13/24-72
    Database cross-references: PDB; 2cpb ; 1-49;, PDB; 2cps ; 1-49;
    Number of features: 0
    /accession=P69541.1
    /start=24
    /end=72
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPZJ2/1-49
    Name: COATB_BPZJ2
    Description: COATB_BPZJ2/1-49
    Number of features: 0
    /accession=P03618.1
    /start=1
    /end=49
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA')
    ID: Q9T0Q9_BPFD/1-49
    Name: Q9T0Q9_BPFD
    Description: Q9T0Q9_BPFD/1-49
    Database cross-references: PDB; 1nh4 A; 1-49;
    Number of features: 0
    /accession=Q9T0Q9.1
    /start=1
    /end=49
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPIF1/22-73
    Name: COATB_BPIF1
    Description: COATB_BPIF1/22-73
    Database cross-references: PDB; 1ifk ; 1-50;
    Number of features: 0
    /accession=P03619.2
    /start=22
    /end=73
    Per letter annotation for: secondary_structure
    Seq('FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA')



```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
```


```python
# Align sequences of interest
align1 = MultipleSeqAlignment([
            SeqRecord(Seq("ACTGCTAGCTAG"), id="Alpha"),
            SeqRecord(Seq("ACT-CTAGCTAG"), id="Beta"),
            SeqRecord(Seq("ACTGCTAGDTAG"), id="Gamma"),
        ])
align2 = MultipleSeqAlignment([
            SeqRecord(Seq("GTCAGC-AG"), id="Delta"),
            SeqRecord(Seq("GACAGCTAG"), id="Epsilon"),
            SeqRecord(Seq("GTCAGCTAG"), id="Zeta"),
        ])
align3 = MultipleSeqAlignment([
            SeqRecord(Seq("ACTAGTACAGCTG"), id="Eta"),
            SeqRecord(Seq("ACTAGTACAGCT-"), id="Theta"),
            SeqRecord(Seq("-CTACTACAGGTG"), id="Iota"),
        ])
```


```python
my_alignments = [align1, align2, align3]
```


```python
print(my_alignments)
```

    [<<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 12) at 7faa11226750>, <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 9) at 7faa11226650>, <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 13) at 7faa11299b10>]



```python
AlignIO.write(my_alignments, "my_example.phy", "phylip")
```




    3




```python
# Add multiple alignments
alignments = AlignIO.parse("my_example.phy", "phylip")
```


```python
for alignment in alignments:
    print(alignment)
    print()
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma
    
    Alignment with 3 rows and 9 columns
    GTCAGC-AG Delta
    GACAGCTAG Epsilon
    GTCAGCTAG Zeta
    
    Alignment with 3 rows and 13 columns
    ACTAGTACAGCTG Eta
    ACTAGTACAGCT- Theta
    -CTACTACAGGTG Iota
    



```python
# print last multiple sequence alignments from phy file (3 records)
alignments = list(AlignIO.parse("my_example.phy", "phylip"))
```


```python
 last_align = alignments[-1]
```


```python
print(last_align)
```

    Alignment with 3 rows and 13 columns
    ACTAGTACAGCTG Eta
    ACTAGTACAGCT- Theta
    -CTACTACAGGTG Iota



```python
first_align = alignments[0]
```


```python
print(first_align)
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma



```python
count = AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.aln", "clustal")
```


```python
# Show converted alignments
print("Converted %i alignments" % count)
```

    Converted 1 alignments



```python
alignments = AlignIO.parse("PF05371_seed.sth", "stockholm")
```


```python
count = AlignIO.write(alignments, "PF05371_seed.aln", "clustal")
```


```python
 print("Converted %i alignments" % count)
```

    Converted 1 alignments



```python
AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.phy", "phylip")
```




    1




```python
AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.phy", "phylip-relaxed")
```




    1




```python
# Number the sequences
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
name_mapping = {}
for i, record in enumerate(alignment):
    name_mapping[i] = record.id
    record.id = "seq%i" % i
```


```python
print(name_mapping)
```

    {0: 'COATB_BPIKE/30-81', 1: 'Q9T0Q8_BPIKE/1-52', 2: 'COATB_BPI22/32-83', 3: 'COATB_BPM13/24-72', 4: 'COATB_BPZJ2/1-49', 5: 'Q9T0Q9_BPFD/1-49', 6: 'COATB_BPIF1/22-73'}



```python
AlignIO.write([alignment], "PF05371_seed.phy", "phylip")
```




    1




```python
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
```


```python
print("Number of rows: %i" % len(alignment))
```

    Number of rows: 7



```python
# exract numbered sequences
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
print(alignment[3:7])
```

    Alignment with 4 rows and 52 columns
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
print(alignment[2,6])
```

    T



```python
print(alignment[2].seq[6])
```

    T



```python
print(alignment[:, 6])
```

    TTT---T



```python
print(alignment[3:6, :6])
```

    Alignment with 3 rows and 6 columns
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49



```python
print(alignment[:, :6])
```

    Alignment with 7 rows and 6 columns
    AEPNAA COATB_BPIKE/30-81
    AEPNAA Q9T0Q8_BPIKE/1-52
    DGTSTA COATB_BPI22/32-83
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49
    FAADDA COATB_BPIF1/22-73



```python
print(alignment[:, 6:9])
```

    Alignment with 7 rows and 3 columns
    TNY COATB_BPIKE/30-81
    TNY Q9T0Q8_BPIKE/1-52
    TSY COATB_BPI22/32-83
    --- COATB_BPM13/24-72
    --- COATB_BPZJ2/1-49
    --- Q9T0Q9_BPFD/1-49
    TSQ COATB_BPIF1/22-73



```python
print(alignment[:, 9:])
```

    Alignment with 7 rows and 43 columns
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    ATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    AKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
edited = alignment[:, :6] + alignment[:, 9:]
```


```python
print(edited)
```

    Alignment with 7 rows and 49 columns
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
edited.sort()
```


```python
# Print based on ID
print(edited)
```

    Alignment with 7 rows and 49 columns
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49



```python
from Bio.Seq import Seq
```


```python
 from Bio.SeqRecord import SeqRecord
```


```python
from Bio.Align import MultipleSeqAlignment
```


```python
 alignment = MultipleSeqAlignment([
            SeqRecord(Seq("ACTCCTA"), id="seq1"),
            SeqRecord(Seq("AAT-CTA"), id="seq2"),
            SeqRecord(Seq("CCTACT-"), id="seq3"),
            SeqRecord(Seq("TCTCCTC"), id="seq4"),
])
```


```python
print(alignment)
```

    Alignment with 4 rows and 7 columns
    ACTCCTA seq1
    AAT-CTA seq2
    CCTACT- seq3
    TCTCCTC seq4



```python
substitutions = alignment.substitutions
```


```python
# Count all pairs and substitutions
print(substitutions)
```

        A    C    T
    A 2.0  4.5  1.0
    C 4.5 10.0  0.5
    T 1.0  0.5 12.0
    



```python
# Forcefully add G's
m = substitutions.select("ATCG")
```


```python
print(m)
```

        A    T    C   G
    A 2.0  1.0  4.5 0.0
    T 1.0 12.0  0.5 0.0
    C 4.5  0.5 10.0 0.0
    G 0.0  0.0  0.0 0.0
    



```python
m = substitutions.select("ACTG")
```


```python
print(m)
```

        A    C    T   G
    A 2.0  4.5  1.0 0.0
    C 4.5 10.0  0.5 0.0
    T 1.0  0.5 12.0 0.0
    G 0.0  0.0  0.0 0.0
    



```python
import Bio.Align.Applications
```


```python
# Show algorithms, some alignments can also be done in the terminal
dir(Bio.Align.Applications)
```




    ['ClustalOmegaCommandline',
     'ClustalwCommandline',
     'DialignCommandline',
     'MSAProbsCommandline',
     'MafftCommandline',
     'MuscleCommandline',
     'PrankCommandline',
     'ProbconsCommandline',
     'TCoffeeCommandline',
     '_ClustalOmega',
     '_Clustalw',
     '_Dialign',
     '_MSAProbs',
     '_Mafft',
     '_Muscle',
     '_Prank',
     '_Probcons',
     '_TCoffee',
     '__all__',
     '__builtins__',
     '__cached__',
     '__doc__',
     '__file__',
     '__loader__',
     '__name__',
     '__package__',
     '__path__',
     '__spec__']




```python
from Bio.Align.Applications import ClustalwCommandline
```


```python
#https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/opuntia.aln
```


```python
align = AlignIO.read("opuntia.aln", "clustal")
```


```python
print(align)
```

    Alignment with 7 rows and 906 columns
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191
    TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191



```python
from Bio import Phylo
```


```python
#https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/opuntia.dnd
```


```python
tree = Phylo.read("opuntia.dnd", "newick")
```


```python
# make a tree file showing evolutionary differences
Phylo.draw_ascii(tree)
```

                                 _______________ gi|6273291|gb|AF191665.1|AF191665
      __________________________|
     |                          |   ______ gi|6273290|gb|AF191664.1|AF191664
     |                          |__|
     |                             |_____ gi|6273289|gb|AF191663.1|AF191663
     |
    _|_________________ gi|6273287|gb|AF191661.1|AF191661
     |
     |__________ gi|6273286|gb|AF191660.1|AF191660
     |
     |    __ gi|6273285|gb|AF191659.1|AF191659
     |___|
         | gi|6273284|gb|AF191658.1|AF191658
         

```python

```
## Pairwise Alignments

