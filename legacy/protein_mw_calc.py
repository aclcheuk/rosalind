from mol_bio_tools import prot_mw_calc 

input_seq = str(input("What is your aa sequence? "))

mol_weight = prot_mw_calc(input_seq)
print(f"The molecular weight of this sequence is: {mol_weight} Da")