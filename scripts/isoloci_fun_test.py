from isoloci_fun import GiniSimpson, HindHe, SwapHap

print(GiniSimpson([50, 50]))
print(GiniSimpson([25, 25, 50]))
print(GiniSimpson([0,0]))
print(GiniSimpson([50,0]))

print(HindHe([[10, 5, 0, 0, 0],
              [3, 0, 10, 0, 0],
              [0, 3, 0, 20, 0]])) # got 0.3456409, same answer in R when done manually

# test out swapping algorithm
myNM = [[4, 5, 7, 1, 6, 0, 8, 3, 2, 3, 2, 1, 2, 2, 2, 2, 5, 2, 2, 2, 1, 2, 3, 2, 5, 5, 11, 2, 2, 3, 4, 5],
        [8, 1, 7, 5, 2, 4, 9, 4, 3, 4, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 4, 3, 2, 2, 8, 1, 3, 0, 1, 2]]
nhap = len(myNM[0])
myAssign = [[i for i in range(nhap) if myNM[0][i] < myNM[1][i]],
            [i for i in range(nhap) if myNM[0][i] >= myNM[1][i]]]
print(myAssign)

for i in range(10):
  myAssign = SwapHap(myNM, myAssign, 80)
  print(myAssign)

