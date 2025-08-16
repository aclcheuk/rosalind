# Initialise variables
adult_rabbits = 1
offspring_rabbits = 0 
n = 6
k = 3
m = 3 # mortality after m months
# rabbit longevity = m
# so 1/m = how much 1 month ages them as a proportion of total lifespan

# Functions:
def month_passing(adult_rabbits:int, offspring_rabbits:int, m:int)->int:
    adult_rabbits_start = adult_rabbits
    offsping_rabbits_start = offspring_rabbits
    for i in range(m):
        adult_rabbits = adult_rabbits + offspring_rabbits
        offspring_rabbits = offspring_rabbits + adult_rabbits
    adult_rabbits = adult_rabbits - adult_rabbits_start


test_offspring, test_adults = m_months_passing(1, 0, 8)

"""
for i in range(n):
    print("Generation:", i+1, "Rabbit count (Adults, Offspring):", adult_rabbits, offspring_rabbits)
    adult_rabbits, offspring_rabbits = month_passing(adult_rabbits, offspring_rabbits, k)
"""

print(adult_rabbits, offspring_rabbits)