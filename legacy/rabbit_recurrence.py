# Initialise variables
adult_rabbits = 1
offspring_rabbits = 0 
n = 31
k = 3

# Functions:
def month_passing(adult_rabbits, offspring_rabbits, k):
    adult_rabbits_new = adult_rabbits + offspring_rabbits
    offspring_rabbits_new = adult_rabbits * k
    return adult_rabbits_new, offspring_rabbits_new

for i in range(n):
    print("Generation:", i+1, "Rabbit count (Adults, Offspring):", adult_rabbits, offspring_rabbits)
    adult_rabbits, offspring_rabbits = month_passing(adult_rabbits, offspring_rabbits, k)


print(adult_rabbits, offspring_rabbits)