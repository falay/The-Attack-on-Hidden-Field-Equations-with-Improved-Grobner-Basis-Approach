from scipy.linalg import lu
import itertools
import sys


def selectPairs(polySet):

	B = set()
	for pair in itertools.combinations(polySet,2):
		B.add( pair )	

	return B		



def degreeParser(term):

	if '^' not in term:
		return 1
	else:
		degree = term[term.index('^')+1:]
	
	return int(degree)	



def LCM(monomial1, monomial2):

	var1 = monomial1.split('*')
	var2 = monomial2.split('*')

	Dict = {}

	for term in var1:
		Dict[term[0]] = degreeParser(term)

	for term in var2:
		if term[0] not in Dict:
			Dict[term[0]] = degreeParser(term)	
		else:
			Dict[term[0]] = max(degreeParser(term),Dict[term[0]])	

	lcm = ""
	for var in sorted(Dict):
		
		lcm += var
		if Dict[var] > 1:
			lcm += "^" + str(Dict[var])

		lcm += "*"	


	return lcm[:-1]	









def Lconstructor(Btmp):

	global q

	L = set()

	for (f1,f2) in Btmp:
		
		largestCommon = lcm(f1.lt(), f2.lt())
		L.add( CAST(largestCommon/f1.lt()*f1) )
		L.add( CAST(largestCommon/f2.lt()*f2) )

	return L	



def select(B):

	dMin = 100
	for (f1,f2) in B:
		dMin = min( dMin, lcm(f1.lt(),f2.lt()).degree() )


	Bselect = set()	
	for (f1, f2) in B:
		
		if lcm( f1.lt(), f2.lt() ).degree() == dMin:
			Bselect.add((f1,f2))

	return Bselect ;		




def vectorConstructor(termList, monSet):

	global Pring

	n = len(monSet)
	coefList = [0]*n

	pos1 = pos2 = 0

	realTermList = []
	for term in termList:
		realTermList.append(sage_eval(term, locals=Pring.gens_dict()))
	realTermList = sorted(realTermList, reverse=True)	

	while pos1 < len(realTermList) and pos2 < n:

		(coef, mon) = coefSpliter(str(realTermList[pos1]))
		curMon = sage_eval(mon, locals=Pring.gens_dict())

		if monSet[pos2] == curMon:
			coefList[pos2] = coef
			pos1 += 1
			pos2 += 1	
		else:
			pos2 += 1

	return coefList			

	

def coefSpliter(term):

	if '*' not in term:
		if term.isdigit():
			return (int(term),'1')
		elif '-' in term	:
			return (-1, term)
		else:
			return (1, term)	
	
	pos = term.index('*')
	coef = term[:pos]

	if coef.isdigit():
		return (int(coef), term[pos+1:])
	elif '-' in term:
		return (-1, term)
	else:
		return (1, term)		




def matrixConstructor(F):

	global Pring, monSet

	termSet = []
	monSet = []

	for poly in F:

		termList = str(poly).split('+')
		termSet.append( termList )

		for term in termList:

			(coef, monomial) = coefSpliter(term)
			mon = sage_eval(monomial, locals=Pring.gens_dict())
			if mon not in monSet:
				monSet.append( mon )

	monSet = sorted(monSet, reverse=True)	

	matrixSet = []
	for termList in termSet:
		matrixSet.append( vectorConstructor(termList, monSet) )

	return matrix(matrixSet)	



def polymonialer(coefList):

	global monSet

	poly = 0
	for i in range(len(coefList)):
		poly += int(coefList[i])*monSet[i]

	return poly
	


# Original version
def FplusConstructor(U, F):

	leadingTermF = set()
	for poly in F:
		leadingTermF.add(poly.lt())

	Fplus = set()	
	for row in U:
		poly = polymonialer(row)
		if poly != 0 and poly.lt() not in leadingTermF:
			Fplus.add(poly)

	return Fplus
	


def Term(polySet):
	
	global Pring 
	
	termSet = set()
	for poly in polySet:
		termList = str(poly).split('+')
		for term in termList:
			termSet.add( sage_eval(term, locals=Pring.gens_dict()) )
	
	return termSet
	



	
def SymbolicPreprocessing(F, Gtmp):

	D = set()
	for poly in F:
		D.add(poly.lt())
	P = set()
	
	while True:
		
		T = Term(F.union(P))
		
		if T == D:
			break
		
		for m in T.difference(D):
			D.add(m)
			
			for g in Gtmp:
				try:
					if m % g.lt() == int(0): 
						P.add(CAST(g*m/g.lt()))
				except TypeError:
					continue
				
	return F.union(P)	


def ReductionF4(L, Gtmp):

	global Pring

	F = SymbolicPreprocessing(L,Gtmp)

	M = matrixConstructor(F)

	pl, U = lu(M, permute_l=True)

	return FplusConstructor(U, F)




def realSet(polyList):
	
	global Pring
	
	polySet = set()
	for poly in polyList:	
		polySet.add(sage_eval(poly, locals=Pring.gens_dict()))

	return polySet
		

def F4algorithm(polyList):

	Gtmp = realSet( polyList )
	Fplus = Gtmp
	B = selectPairs( Gtmp )

	while len(B) > 0:

		Btmp = select(B)
		
		if len(Btmp) == 0:
			continue
		
		B = B.difference(Btmp)
		L = Lconstructor(Btmp)
		Fplus = ReductionF4(L, Gtmp)	

		for poly in Fplus:
			for g in Gtmp:
				B.add( (poly,g) )
			Gtmp.add( poly )


	return Gtmp	
	
