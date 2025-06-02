# from Bio import SeqIO # pip install biopython
# Try to create algorithms without using Biopython toolkit

# Prints total length of sequence and counts for each base
def base_count(input_seq:str)->int:
    caps_seq = input_seq.upper()
    count = {"A":0,
             "T":0,
             "C":0,
             "G":0}
    for base in caps_seq:
        if base in count:
            count[base] +=1
        else:
            pass

    print(count)
    print(f"Total length =", len(caps_seq))
    return len(caps_seq)

# Prints the GC concent as a percentage, also returns dictionary showing counts for each base
def gc_content(input_seq:str)->float:
    caps_seq = input_seq.upper()
    count = {"A":0,
             "T":0,
             "C":0,
             "G":0}
    for base in caps_seq:
        if base in count:
            count[base] +=1
        else:
            pass

    gc_content = ((count["C"] + count["G"]) / len(caps_seq)) *100
    #print(count)
    #print(f"GC content = {gc_content:.4f} %")
    return gc_content

# Takes DNA input and returns transcribed RNA sequence
def transcribe(input_seq:str)->str:
    caps_seq = input_seq.upper()
    mrna = [*caps_seq] # separates string into list of letters using unpack [*] method

    """print(len(mrna))
    print(mrna)
    used to test if mrna is the right length and a list"""

    for i in range(len(mrna)):
        if mrna[i] == "T":
            mrna[i] = "U"
        else:
            pass
    mrna = "".join(mrna) # joins each element of list into a string
    print("Transcribed mRNA sequence:", mrna)
    return mrna

# Reverses transcription (replacing all 'U's with 'T's)
def rev_transcribe(input_seq:str)->str:
    caps_seq = input_seq.upper()
    dna = [*caps_seq] # separates string into list of letters using unpack [*] method

    for i in range(len(dna)):
        if dna[i] == "U":
            dna[i] = "T"
        else:
            pass
    dna = "".join(dna) # joins each element of list into a string
    print("DNA sequence:", dna)
    return dna

#Reverse Complement:
def rev_comp(input_seq:str)->str:
    caps_seq = input_seq.upper()
    reverse = caps_seq[::-1] # reverses sequence
    reverse = [*reverse] # separates bases into list items

    for i in range(len(reverse)): # iterates through list and transcribes base
        if reverse[i] == "T":
            reverse[i] = "A"
        elif reverse[i] == "A": 
            reverse[i] = "T"
        elif reverse[i] == "C":
            reverse[i] = "G"
        elif reverse[i] == "G":
            reverse[i] = "C"
        else:
            pass
    reverse = "".join(reverse) # joins each element of list back into a single string
    return reverse

#Translate DNA or RNA into amino acid sequence (one letter code):
def translate(input_seq:str)->str:
    codon_table = {
    # Codon table in dictionary format
    # 'M' - START, '*' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "*", "TAG": "*", "TGA": "*"}
    caps_seq = input_seq.upper()
    codon_list = [(caps_seq[i:i+3]) for i in range(0, len(caps_seq), 3)]
    aa_list = codon_list
    for i in range(len(codon_list)):
        aa_list[i] = codon_table[codon_list[i]]

    aa_seq = "".join(aa_list)
    return aa_seq

    # Issues to fix and functionality to add: 
    # Issue 1: cannot handle DNA sequences that contain a number of bases that is not divisible by 3.
    # Issue 2: does not selectively begin translation at start codon (ATG)
    # Functionality 1: give predicted MW and pI of translated peptide
    # Functionality 2: identify all ORF and give alternative peptides that can be translated 

# Check DNA sequence:
bases = ["A", "C", "T", "G", "U", "N", "X"]
def is_dna_seq(input_seq:str)->bool:
    input_seq = input_seq.upper()
    for base in input_seq:
        if base in bases:
            continue
        else:
            print("Please input DNA sequence.")
            return False
    return True

"""
#Reads FASTA files and outputs python dictionary with each entry {>name1:sequence1,>name2:sequence2}
def read_fasta(input_file):
    #make sure you import SeqIO from Biopython with "from Bio import SeqIO"
    seq_dict = {}
    for seq_record in SeqIO.parse(input_file, "fasta"): # SeqIO.parse() is able to parse files
        print(seq_record.id)
        print(repr(seq_record.seq)) #print(repr(seq_record.seq)) 
        #repr() prints a smaller snippet so it doesn't clutter the terminal
        print(len(seq_record))
        seq_dict[seq_record.id] = seq_record.seq # adds to dictionary seq_record.id:seq_record.seq
    return seq_dict
"""

def parse_fasta(fasta_file):
    with open(fasta_file) as fasta_content:
        sequence = ""
        for line in fasta_content:
            line = line.rstrip()
            if line.startswith(">"):
                descriptor = line
            else:
                sequence = sequence + line 
        yield descriptor, line

#code to test the different functions:
if __name__ == "__main__":
    test = "AtGAACCCTTTGGGggggggggacacttggaATcGaTAAATTCCGGggctatat"
    base_count(test)
    gc_content(test)
    transcribe(test)
    rev_comp(test)
    translate(test)
    test1 = parse_fasta("rosalind_gc.txt")
    print(test1)

""" used to solve GC content problem in Rosalind Stronghold
my_dict = read_fasta("rosalind_gc.txt")
for seq in my_dict:
    print(f"{seq} = {gc_content(my_dict[seq])}")
"""
