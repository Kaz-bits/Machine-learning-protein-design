# Import modules
from goose import create
import numpy as np

# Create a variable with your sequence of reference
test = "PPPEDPPPAPPTPAELAALAAARAAAAARAAAAAPRLQRRKARRRKRERKKPGGKAKPARPRPPKRRRPKKRKGGKRARRKRKRKIKKKKKRKRRKKRKRAERRKAAREAARAAAAAALASA"

# Create empty list
temp_list = []

# Iterate over the same sequence with a list of parameters (kappa)
for i in np.arange(0.09, 0.7, 0.0003):
        temp = create.kappa_var(test, kappa = i)
        temp_list.append(temp)

# Save list of designed disorderd regions in txt file
with open("IDRBS_023_DN2_goose_kappa.txt", "w") as f:
    print(temp_list, file = f)