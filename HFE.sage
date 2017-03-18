from random import randint
import hashlib
import time
import sys

global hashLen


def HFEkeyGen(n, q):
	
	global R, K, P, L, Q, variables, g


	K = GF(q)
	variables = ['x'+str(i) for i in range(1, n+1)]
	R = PolynomialRing(K, variables, n) 
	J = ideal([x^q-x for x in R.gens()])	
	H = R.quotient_ring(J,variables) 
	

	
	# S = (Ms, Vs)
	Ms = random_matrix(K, n, n) 
	while Ms.is_singular():
		Ms = random_matrix(K, n, n)

	Vs = random_vector(K, n)


	# T = (Mt, Vt)	
	Mt = random_matrix(K, n, n) 
	while Mt.is_singular():
		Mt = random_matrix(K, n, n)

	Vt = random_vector(K, n)

	privateKey = ((Ms,Vs),(Mt,Vt))


	# To construct the public keys
	Mesg = vector(R, n, R.gens()) ;
	
	# Apply S
	Mesg = Ms*Mesg + Vs

	
	mesgList = Mesg.list()
	mesgList.reverse()
	
	P.<x> = PolynomialRing(K)
	g = P.irreducible_element(n)
	L.<t> = PolynomialRing(R)
	poly = L(mesgList)
	#print "before: ", poly

	g = L(g)
	I = L.ideal([g])
	Q = L.quotient_ring(I)
	
	poly = Q(poly^(q+q)+poly^q+1)
	
	#print "after: ", poly
	polyList = poly.list()
	polyList.reverse()
	publicKey = vector(R, n, polyList)

	publicKey = Mt * publicKey + Vt

	return (publicKey, privateKey)	



def HFEencryption(plaintext, publicKey, n):

	if len(plaintext) != n:
		sys.exit("invalid usage of encryption")

	return publicKey(plaintext).list()	




def coefSpliter(poly):

	global n

	coefList = [0]*n
	termList = poly.split(' + ')

	for term in termList:

		term.strip()

		if term.startswith('t'):
			
			if '^' in term:
				coefList[n-1-int(term[2:])] = 1
			else:
				coefList[n-2] = 1


		elif 't' in term:
			
			coef = int(term[0:term.index('t')-1])

			if '^' in term:
				coefList[n-1-int(term[term.index('t')+2:])] = coef	
			else:
				coefList[n-2] = coef


		else:
			coefList[n-1] = int(term)	

	return coefList		



def HFEdecryption(ciphertext, privateKey, q):
	
	((Ms,Vs),(Mt,Vt)) = privateKey

	# apply T^-(1)
	cipher = Mt^(-1)*(vector(R, n, ciphertext) - Vt) 


	# solve P(x) = cipher	
	extensionField.<t> = GF(q^n, modulus=g)
	ring.<x> = extensionField[]
	
	cipherElement = 0
	for i in range(n):
			cipherElement += extensionField(cipher[i])*t^(n-1-i) 		

	p = x^(q+q) + x^q + 1 - cipherElement

	roots = p.roots()

	if len(roots) == 0:
		sys.exit("empty baby")

	rootsList = []	
	for root in roots:
		rootsList.append( coefSpliter(str(root[0])) ) 



	# apply S^(-1)
	plaintextSet = []
	for root in rootsList:
		plaintextSet.append( Ms^(-1)*(vector(R, n, root) - Vs) )


	print "possible plaintext: "
	for plaintext in plaintextSet:
		print plaintext	



		
def Authenticator(solutionSet): 		
	
	global hashLen
	
	for dict in solutionSet:
		List = sorted(dict.items(), reverse=True)
		hashCode = [ pair[1] for pair in List[:hashLen] ]
		plaintext = [ pair[1] for pair in List[hashLen:] ]
		hex = hashlib.sha256( str(plaintext) ).hexdigest()
		
		if map(int, bin(int(hex[:3],16))[2:]) == hashCode:
			return plaintext


def DegreeReduction( sysEquations ):
	
	possibleExp = ["^2", "^3", "^4"]
	
	for removedItem in possibleExp:
		sysEquations = sysEquations.replace(removedItem, "") 
			
	return sysEquations


			
		
def GrobnerBasisBreaker(pubKey,ciphertext,q):

	equations = ""
	variables = ""
	n = len(pubKey)

	for i in range(n):
		variables += 'x' + str(i+1) + ','
		equations += str(pubKey[i]) + "-" + str(ciphertext[i]) + ","

	equations = equations[:-1]
	variables = variables[:-1]
		
	Pring = PolynomialRing(GF(q),variables,order='lex')
	Pring.inject_variables(verbose=False)
	
	# Optimization
	
	
	equations = DegreeReduction( equations ) 
	
	I = Ideal( sage_eval(equations, locals=Pring.gens_dict()) )
	J = I + sage.rings.ideal.FieldIdeal( Pring )

	return Authenticator( J.variety() )

 


 
 
# Main
if len(sys.argv) != 2:
	sys.exit("Usage: ./sage HFE.sage [plaintext length]")


# generating the plaintext
	
n = int(sys.argv[1])
plaintext = []
for i in range(n):
	plaintext.append( randint(0,1) )

#print "Random generated plaintext:"
#print plaintext

# appending the hash code (using SHA256)
hex = hashlib.sha256( str(plaintext) ).hexdigest()
hashCode = map(int, bin(int(hex[:3],16))[2:])
sendText = hashCode + plaintext

hashLen = len( hashCode )

#print "Send text len:", len(sendText)

(pubKey, prvKey) = HFEkeyGen(len(sendText), 2)
ciphertext = HFEencryption(sendText,pubKey,len(sendText))


start_time = time.time()

brokenText = GrobnerBasisBreaker(pubKey, ciphertext, 2)

elapsed_time = time.time() - start_time


if brokenText == plaintext:
	print "Success, time =", elapsed_time
else:
	print "Fails"	

#HFEdecryption(ciphertext, prvKey, p^k)

