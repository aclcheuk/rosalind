# Given an integer, n, return the total number of permutations
# and the list of all possible permutations of integers, 1 to n.
from itertools import permutations

def get_perms(input_list):
    perm_list = []
    for perm in permutations(input_list):
        perm_list.append(perm)
    return perm_list

ros_list = []
for i in range(6):
    ros_list.append(i+1)


perm_list = get_perms(ros_list)
print(len(perm_list))
for perm in perm_list:
    print(*perm)


"""
test_list =['a', 'b', 'c']
    perm_list = get_perms(test_list)
print(len(perm_list))
for perm in perm_list:
    print(perm)
"""
