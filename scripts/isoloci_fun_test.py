from isoloci_fun import GiniSimpson, HindHe

print(GiniSimpson([50, 50]))
print(GiniSimpson([25, 25, 50]))
print(GiniSimpson([0,0]))
print(GiniSimpson([50,0]))

print(HindHe([[10, 5, 0, 0, 0],
              [3, 0, 10, 0, 0],
              [0, 3, 0, 20, 0]])) # got 0.3456409, same answer in R when done manually
