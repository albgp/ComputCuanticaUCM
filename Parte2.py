
from auxiliaryFns import *


probs={} #We'll store the probabilities in a dictionary with the symbols as keys.
probs["-"]=0.1927;probs["A"]=0.0575;probs["B"]=0.0128;probs["C"]=0.0263;probs["D"]=0.0285;probs["E"]=0.0912
probs["F"]=0.0173;probs["G"]=0.0133;probs["H"]=0.0313;probs["I"]=0.0599;probs["J"]=0.0006;probs["K"]=0.0084
probs["L"]=0.0335;probs["M"]=0.0235;probs["N"]=0.0596;probs["O"]=0.0689;probs["P"]=0.0192;probs["Q"]=0.0008
probs["R"]=0.0508;probs["S"]=0.0567;probs["T"]=0.0706;probs["U"]=0.0334;probs["V"]=0.0069;probs["W"]=0.0119
probs["X"]=0.0073;probs["Y"]=0.0164;probs["Z"]=0.0007

assert len(probs)==27 #Checking input is ok.


#Function to compute entropy given a probability dictionary (only probabilities are needed)
calculateEntropy = lambda probs: sum([-np.log2(p)*p for _,p in probs.items()]) 

print("La entropía calculada para las probabilidades siguientes\n{}\nes {}".format(probs,calculateEntropy(probs))) #Just priting









#####################################
#          Base 2 code              #
#####################################
 

def myHuffman2(probs): #This function computes the huffman binary code for a given symbol list alongside their probability, and returns it in a dict 
	# indexed by symbol.
	current=list(probs.keys()) #We start with a exhaustive list of symbols
	code={i:"" for i in current} #And the dictionary that will later store the code (now empty).

	while len(current)>=2: #While there are at least two elements left
		a,b=getMinProbability(probs,current,2) #We get the two elements with minimum probability
		a,b=swapAccordingTo2(a,b,probs) #And order them to match the convention stated in class and in the provided slides.
		for char in a: code[char]+="0"  #To the symbols in the first element, append a "0" char to the code of that symbol
		for char in b: code[char]+="1" #and for the ones in the 2nd, a "1" char
		current.remove(a); current.remove(b) #Once they're processed, we delete them from the list
		current.append(a+b) #and append the concatenation of them both to the list
		#Note: we do not calculate the sum of probabilities for the new element, as we'll recalculate it every time
		# it's inefficient af but we want to keep it simple.

	for k,v in code.items(): code[k]=v[::-1] #Last, for every code, we must reverse it, as we've created it backwards
	return code #And that's it.



code2=myHuffman2(probs) #Calculating the code

print ("Símbolo & Probab & HuffmanCod\\\\") #Printing it
for k,v in code2.items():
    print ("{}&{}&{}\\\\".format(k,probs[k],v))

### Longitud Media. Input: dict indexed by char containing arrays w/ its Huffman code and the probability dict w/ ocurrence probability

length = lambda probs, codes: sum(len(v)*probs[k] for k,v in codes.items()) #Function to compute the "mean length" of a huffman code

print("La longitud media del código Huffman es: {}".format(length(probs,code2))) #Just printing, again

### Checking the Kraft-McMillan identity for base n encoding

tol=1E-12
KMid = lambda base,codes : sum([base**(-len(v)) for _,v in codes.items()])<=1+tol #Returns true if the K-M identity holds for the code "codes" in base "base"
#(up to floating point precission with tolerance).

print(sum([2**(-len(v)) for _,v in code2.items()]))

if KMid(2,code2):
	print("La igualdad de Kraft-McMillan se satisface para el caso binario.")
else:
	print("La igualdad de Kraft-McMillan NO se satisface para el caso binario.")












####################################
#            Base 3 code           #
####################################


def myHuffman3(probs): #Ternary huffman code generator function
	current=list(probs.keys())
	code={i:"" for i in current}


	while len(current)>1:
		if len(current)==2: #Special case in which there are only two groups of symbols left to process.
			a,b=getMinProbability(probs,current,2) #Which is just the same as the binary treatment of the array
			a,b=swapAccordingTo2(a,b,probs)
			for char in a: code[char]+="0"
			for char in b: code[char]+="1"
			current.remove(a); current.remove(b)
			current.append(a+b)

		else:
			a,b,c=getMinProbability(probs,current,3)
			a,b,c=swapAccordingTo3(a,b,c,probs)
			for char in a: code[char]+="0"
			for char in b: code[char]+="1"
			for char in c: code[char]+="2"
			current.remove(a); current.remove(b); current.remove(c)
			current.append(a+b+c)

	for k,v in code.items():
		code[k]=v[::-1] #And finally, reverse it

	return code 

code3=myHuffman3(probs)

print ("Símbolo&Probab&HuffmanCod\\\\")
for k,v in code3.items():
    print ("{}&{}&{}\\\\".format(k,probs[k],v))

print("La longitud media del código Huffman es: {}".format(length(probs,code3))) #Just printing, again

print(sum([3**(-len(v)) for _,v in code3.items()]))

if KMid(3,code3):
	print("La igualdad de Kraft-McMillan se satisface para el caso ternario.")
else:
	print("La igualdad de Kraft-McMillan NO se satisface para el caso ternario.")











####################################
#       Encoding and so on         #
####################################

encodeMsg=lambda code, msg: "".join([code[char] for char in msg]) #Given a code and a msg, it returns the encoded msg, C(msg)

msgToEncode="A-SINGLE-LOOK-AT-THEM-IS-ENOUGH-TO-SHOW-THEY-COULD-ONLY-BE-WRITTEN-DOWN-BY-A-MATHEMATICIAN-OF-THE-HIGHEST-CLASS"
encodedMsg=encodeMsg(code3,msgToEncode) 
print("El mensaje codificado es: {}".format( encodedMsg )) #Easy 

n=int(encodedMsg,3) #Previous huffman sequence read as a base 3 integer
assert len(str(n))==141 #Checking everything is going ok

print("El número anterior corresponde al número decimal n={}".format(n))














###################################
#           Decoding              #
###################################

def decodeMsg(code, msg): #Method to decode a msg according to a huffman code for its alphabet
  #Inefficient, but ok, giving up the cool efficient tree structure decomposition.

  #Assuming the input msg corresponds to a code-encoded original, otherwise the algorithm may never end.
  decodedmsg=""
  while msg!="": #While not fully decoded
    for k,v in code.items(): #We look for the matching huffman code at the begining of the remaining sequence
      if msg.startswith(v): #Found it!
        decodedmsg+=k #Append the char to the decoded msg
        msg=msg[len(v):] # and delete the already decoded part from the original seq
        break # Go find the next one
  return decodedmsg 





msgToDecode="110011101001100010111111100110111011000101000\
1000101111111001101110011011110110011000111001000010011101\
0001100101111001011000111011100000100111010011010111000101\
1011011100111001101111011000111100111001111110000111001110\
11100111010110110111001000010010111111110010101111110110111"

print("El mensaje decodificado es: {}".format(decodeMsg(code2,msgToDecode)))














###################################
#           Factoring             #
###################################




def MillerRabinTest(n, verbose=False):# Miller-Rabin primality test. To avoid complications, works well for n>10
	#Returns: False if composite, True if maybe a prime.
    assert n==int(n)

    #We look for the decomposition of n-1=d*2^s with larger s
    s = 0
    d = n-1
    while d%2==0:
        d>>=1 #d//=2 but faster
        s+=1
    assert(2**s * d == n-1)
 
 	#The following conditions hold for primes bc if n is prime, then by Fermat's little theorem, a^(n-1)=1 (mod n). By the definition of prime in a field,
 	#and the decomposiion x^2=(x-1)(x+1), if n=2^rd, as calculated, then a^d=1 (mod n) or a^(d2^s)=-1. If it is not satisfied by any of the s and d that define
 	# the decomposition, then it cannot be a prime.
    def isComposite(a): #Checking if any prime-conditioned constraint holds, then we cannot conclude that a is prime
    	# The iterative powering must be done in a modular way (otherwise it'd be exponential in time). Pow allows it.
        if pow(a, d, n) == 1: return False #Equality a^d (mod n) =1 => n may be prime 
        for i in range(s): #As s was the larger, we only need to try i-values up to s-1.
            if pow(a, 2**i * d, n) == n-1: return False #a^(d2^i) mod(n)=n-1=-1 (mod n) => may be prime
        if verbose:
        	print("El test de Miller-Rabin ha hallado un compuesto con semilla a={}".format(a))
        return True  #If it does not satisfy any of the above contraints,
        #then, by contrapositive, it can't be prime.
 
    for i in range(15):# If n is probably prime, number of chances to find the composite answer. Success probability: 1-1/4^(15). Fckin high.
        a = random.randrange(2, n) #Choosing a seed
        if isComposite(a): #If we know for sure that it's not a prime, we return composite
            return False
 
    return True  #Else, returns probably prime. (Deterministic for negative answers, probabilistic for positive ones)



primesFound=[2,5,653,348527,551569] #Prime divisors obtained using the Brent algorithm
# https://maths-people.anu.edu.au/~brent/pd/rpb051i.pdf
# alongside the sympy implementation of the Pollard-rho1 algorithm
# the brute force method would also work for these ones. 

#Not the coolest algorithm for prime checking, as it works in O(2^(n/2)) in the worst case scenario for a n-bit int, but it's okay.
isPrime = lambda n: not any(n%k==0 for k in range(2,int(sqrt(n)-1)+2)) 

for i in range(2,10000000):
	if (n%i==0) and isPrime(i):
		print(i)

m=n
for p in primesFound: #Getting rid of the anoying prime divisors that we already know
	while (m%p==0):
		m//=p 

print("m={}".format(m))



for p in primesFound: #Checking the primes are correct, and that they are, indeed, primes
	assert isPrime(p)
	assert n%p==0

print("El número anterior es divisible por los siguientes primos: ", end="")
for prime in primesFound:
	print(prime, end=", ")
print()


if not MillerRabinTest(m, verbose=True) and not isPowerOfPrime(m):
	print("El número que obtenemos dividiendo iterativamente por los primos que hemos hallado hasta el momento, m,\
 no es primo, lo que obtenemos aplicando el algoritmo de Miller-Rabin. Además, como los primos mencionados no pueden dividir a ningún factor\
  de m, tenemos que deben existir al menos dos otros primos (puesto que m no es potencia de primo)\
  que dividan al número m y, por tanto, al número original n. n tiene, por tanto, al menos 7 divisores primos distintos.")
else:
	assert False
	#Hoping not to reach this point. In that case, we'd need to use the AKS algorithm, or the deterministic MR, for primality test
	#which can be implemented in polynomial time.










###################################
#           Hamming codes         #
###################################


# The code is H_4 which is a linear code of type [n=15=2^4-1,k=n-4=11] over the field F2

H=np.array([
	[0,0,0,0,0,0,0,1,1,1,1,1,1,1,1],
	[0,0,0,1,1,1,1,0,0,0,0,1,1,1,1],
	[0,1,1,0,0,1,1,0,0,1,1,0,0,1,1],
	[1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]]) 

G=np.array([
	[1,1,0,1,0,0,0,1,0,0,0,0,0,0,1],
	[0,1,0,1,0,0,0,1,0,0,0,0,0,1,0],
	[1,0,0,1,0,0,0,1,0,0,0,0,1,0,0],
	[0,0,0,1,0,0,0,1,0,0,0,1,0,0,0],
	[1,1,0,0,0,0,0,1,0,0,1,0,0,0,0],
	[0,1,0,0,0,0,0,1,0,1,0,0,0,0,0],
	[1,0,0,0,0,0,0,1,1,0,0,0,0,0,0],
	[1,1,0,1,0,0,1,0,0,0,0,0,0,0,0],	
	[0,1,0,1,0,1,0,0,0,0,0,0,0,0,0],
	[1,0,0,1,1,0,0,0,0,0,0,0,0,0,0],
	[1,1,1,0,0,0,0,0,0,0,0,0,0,0,0]])


D="01000101011010001011001111101011110011100111100111001100\
11110011011111110100010100000001101101100001011011001110\
11001111011110001110000101101010111111100110111111011100\
01011011110010110100101110101001011010111010110100101011\
1101101010101011"

nBlocks=len(D)//15 #How many 15-bits codes are in the sequence

D=np.array([int(i) for i in D ])

D=np.reshape(D, (nBlocks,15)) #We split D into nBlocks vecs with length 15


#Syndromes

syndromes=[H.dot(v)%2 for v in D] #For each vector v in D, we calculate de syndrome Hv, where the product is performed in F_2~Z_2.
#these will be, in binary, the positions where the transmition errors occurred.

Dlimpias=D.copy() 

for i, syn in enumerate(syndromes): #For each calculated syndrome
	position=binaryArrayToInt(syn) #Let's take the integer that it represents in Z^4_2
	Dlimpias[i][position-1]+=1 #and correct it in Dlimpias. The -1 stands for the possition correction, as python indexes from 0.

Dlimpias%=2 # The operation was in Z_2 so we need to assure that we didn't leave Z_2. This corrects it if necessary

print("Dlimpias=",Dlimpias.reshape(nBlocks*15)) #Print the new D vector (Dlimpias) in a human readable way.
print("Dlimpias="+"".join([str(i) for i in Dlimpias.reshape(nBlocks*15)]))

#Apartado B 

okCols=np.array([2,4,5,6,8,9,10,11,12,13,14]) #Indexed from 0! Slice for slicing, taking out the 2^k-th value for k=0,1,2,3.
#Just bc all power-of-two positions are parity (error correction) bits in Hamming coding.

Dlimpias=np.array([i[okCols] for i in Dlimpias]) #Take the interesting columns from the vectors in Dlimpias

Dlimpias=np.array([i[::-1] for i in Dlimpias]) #We reverse them because the msgs are encoded in inverse order.

n=binaryArrayToInt(Dlimpias.reshape(nBlocks*11)) #Read it in binary
print("El número codificado en la secuencia original, recuperado con el control de errores basado en códigos Hamming y en base 10 es: {}".format(n)) 









########################################
#            Ejercicio 2               #
########################################

# rho0=3/4*multiKron(Pa,(1-Pb),Pc)+1/4*multiKron(1-Pa,Pb,Pc)

# TOL=1E-7
# if np.abs(np.matrix.trace(np.matmul(rho0,rho0))-1)<TOL:
# 	print("El estado rho(0) es puro.")
# else:
# 	print("El estado rho(0) es un estado mixto.")


# vonNeumanEntropy = lambda M : -np.matrix.trace(np.matmul(M, logm(M))) #Natural log, isnt it?


# traceDistance = lambda M,N : 1/2*np.matrix.trace( sqrtm( matmul(M-N,M-N) ) )

# fidelity = lambda M,N: np.matrix.trace(sqrtm( multiProd(sqrtm(M), N, sqrtm(M)) ) )






