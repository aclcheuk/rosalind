import mol_bio_tools as mbt
bases = ["A", "C", "T", "G", "U", "N", "X"]

def is_dna_seq(input_seq:str)->bool:
    input_seq = input_seq.upper()
    for base in input_seq:
        if base in bases:
            continue
        else:
            print("Please input ")
            return False
    return True


input_seq = "AAcctt"

print(is_dna_seq(input_seq))
#def clean_dna_seq(input_seq:str)->str:


