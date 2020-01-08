#####################################
#          Auxiliary fns            #
#####################################
from itertools import count, islice
from functools import reduce
import random 
from math import gcd, log2
import time
import string
from math import sqrt
import numpy as np 
from scipy.linalg import logm, expm, sqrtm
import sys
import matplotlib.pyplot as plt
from sympy.physics.quantum import TensorProduct as tp
import sympy as sp
from sympy import Matrix, Rational, eye, latex, I
from sympy.physics.quantum.dagger import Dagger
from matplotlib import rc
import qutip

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

rc('text.latex', preamble=r'\usepackage{amsmath}')

import warnings
warnings.filterwarnings("ignore")

from sympy.assumptions.assume import global_assumptions
from sympy import AppliedPredicate, Q

verbose=False


strProb = lambda st, probs: sum([probs[char] for char in st]) #This function returns, for a given string, the sum of the probabilities
#associated to each string symbol. Useful for calculating the Huffman code, as we'll later see. As a lambda function

def getMinProbability(probs, current, k): #As in the huffman algorithm for the k-ary code, we'll need to find 
	# the k elements with the lowest probability in order to assign them a new digit in the code, this function looks for the k elements with
	# the min probability value (in the sense as the function above) and returns them
	sortCurr=sorted(current,key=lambda x: strProb(x,probs)) #To that extent, it sorts the "current" array passed as argument
	return sortCurr[:k] #And returns the first k elements after sorting

def swapAccordingTo2(a,b,probs): #Auxiliary fn for the huffman code method.
	#Sort the two provided elements in order to assign 0 to the element with the symbol that appears first in
	# the alphabet A. (Conventions)
	for char in probs.keys():
		if char in a:
			return a,b
		elif char in b:
			return b,a
	assert False # Hopping not to reach

def swapAccordingTo3(a,b,c,probs): 
	#The same as the previous one but w/ three elements, bc now we want to assign "0" to the 1st element, "1" to the 2nd and "2" to the 3rd.
	ret=[]
	for char in probs.keys():
		if char in a and a not in ret:
			ret.append(a)
		elif char in b and b not in ret:
			ret.append(b)
		elif char in c and c not in ret:
			ret.append(c)
		if len(ret)==3:
			break 

	return ret

from math import log2

#We'll know, using the MR algorithm, that the number is composite BUT, it may be a power of a prime, and then we'd have only one prime factor in m.
# To be sure that it has at least two of them, let's check, in low time complexity if that is the case.
def isPowerOfPrime(n): 
	for i in range(2,int(log2(n))+1): #If p^i=n, then i<=log2(n)
		possiblePrime=int(pow(n,1/i)) #Get the integer closest to the power n^(1/i)
		if pow(possiblePrime+1, i)==n: #We assume floating point errors, so we check for x= p-1,p and p+1. If it turns that x^i=n
		#then we've found that n is a power. If x is prime (we know it to incredible precission using MR algorithm), then it's a power of prime
			if MillerRabinTest(possiblePrime+1): return True
		elif pow(possiblePrime, i)==n:
			if MillerRabinTest(possiblePrime): return True
		elif pow(possiblePrime-1, i)==n:
			if MillerRabinTest(possiblePrime-1): return True
	return False #Else it's not

binaryArrayToInt = lambda x: int("".join([str(i) for i in x]),2)

multiKron = lambda *args : reduce(np.kron, args) #Nested kronecker products are ugly. this function performs the kronecker product between more than two matrices. 
multiProd = lambda *args : reduce(np.matrix.matmul, args) #Same with ordinary matrix product

def decomposeProductTensor(mat, dims): #Decompose a nmxnm matrix into a nxnxmxm matrix. Matrix, n and m as the function input.
	a,b=dims
	return mat.reshape([a, b, a, b]) 

def myPartialTrace4x3in2ndSpace(mat): #Calculate the partial trace over the 3rd (C) 3x3 system.
	reshaped_matrix = decomposeProductTensor(mat, [4,3]) #Decompose it
	tracedMatrix= np.einsum('jiki->jk', reshaped_matrix) #Sums over repeated indices, calculates the trace on the 2nd subsystem. (Returns a 4by4 density matrix on the AB system)
	return tracedMatrix

def myPartialTraspose4x3in2ndSpace(mat): #Calculate the partial traspose over the 3rd (C) 3x3 system.
	reshaped_matrix = decomposeProductTensor(mat, [4,3]) #Decompose it
	trMatrix= np.einsum('ijkl->ilkj', reshaped_matrix) 
	return trMatrix.reshape([12,12])

def myPartialTraspose2x2in2ndSpace(mat): #Calculate the partial traspose over the 2nd (B) 2x2 system.
	reshaped_matrix = decomposeProductTensor(mat, [2,2]) #Decompose it
	trMatrix= np.einsum('ijkl->ilkj', reshaped_matrix) 
	return trMatrix.reshape([4,4])











