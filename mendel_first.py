"""k, m, n individuals in a population"""

def calc_offspring_probs(homo_dominant:int, hetero:int, homo_recessive:int)->float:
    k = homo_dominant
    m = hetero
    n = homo_recessive
    total_pop = k + m + n

    prob_k1 = k/total_pop
    prob_m1 = m/total_pop
    prob_n1 = n/total_pop

    prob_k1_k2 = prob_k1 * ((k-1)/(total_pop-1))
    prob_k1_m2 = prob_k1 * ((m)/(total_pop-1))
    prob_k1_n2 = prob_k1 * ((n)/(total_pop-1))
    
    prob_m1_k2 = prob_m1 * ((k)/(total_pop-1))
    prob_m1_m2 = prob_m1 * ((m-1)/(total_pop-1))
    prob_m1_n2 = prob_m1 * ((n)/(total_pop-1))

    prob_n1_k2 = prob_n1 * ((k)/(total_pop-1))
    prob_n1_m2 = prob_n1 * ((m)/(total_pop-1))
    prob_n1_n2 = prob_n1 * ((n-1)/(total_pop-1))
    
    # check probabilities sum to 1
    parent_probs = [prob_k1_k2,prob_k1_m2,prob_k1_n2, prob_m1_k2, prob_m1_m2, prob_m1_n2, prob_n1_k2,prob_n1_m2,prob_n1_n2]
    sum_parent_probs = sum(parent_probs)
    print("Sum of parent probabilities: ", sum_parent_probs)
    

    # Offspring Genotype Probabilities
    prob_DD = (prob_k1_k2) + (0.5*prob_k1_m2) + (0.5*prob_m1_k2) + (0.25*prob_m1_m2)
    prob_Dd = (0.25*prob_k1_m2) + (prob_k1_n2) + (0.25*prob_m1_k2) + (0.5*prob_m1_m2) + (prob_n1_k2) + (0.5*prob_n1_m2)
    prob_dd = (0.25*prob_k1_m2) + (0.25*prob_m1_k2) + (0.25*prob_m1_m2) + (prob_m1_n2) +(0.5*prob_n1_m2) + (prob_n1_n2)
    offspring_probs = [prob_DD, prob_Dd, prob_dd]
    sum_offspring_probs = sum(offspring_probs)
    print("Sum of offspring probabilities: ", sum_offspring_probs)
    return (1-prob_dd)

k = 2
m = 2
n = 2

a = calc_offspring_probs(k,m,n)
print(a)

prob_dom = calc_offspring_probs(k,m,n)
print("Probability of Dominant offspring = ", prob_dom)


