# Import modules 
from goose import create
import numpy as np

# Create an empty list
temp_list = []

# Generate different sequences with a fixed length
for i in np.arange(0.01, 0.9, 0.0001):
    temp = create.sequence(50, kappa = i)
    temp_list.append(temp)

# Save list of designed sequences in txt file
    with open("IDRBS_DN50_goose_kappa.txt", "w") as f:
        print(temp_list, file = f)