from sage.all import *

def discreteLog(P, Q, E, p, a, b):
    # Let's create fields, compute lifts, find psi, and finish to get an answer
    Field = GF(p)
    Qandp = Qp(p, prec = 3, type = 'capped-rel', print_mode = 'series')
    x1 = P[0]
    y1 = P[1]
    integerX = Q(Integer(x1))
    integerY = Integer(y1)
    h0 = integerY
    numH1 = ((y1 ** 2) - (x1 ** 3) - (a * x1) - b) / p
    h1 = int( - Field(numH1) / Field(2 * y1) )
    numH2 = (( (y1 + h1*p) ** 2) - (x1 ** 3) - a * x1 - b) / (p ** 2)
    h2 = int( - Field(numH2) / Field(2 * (y1 + h1*p)) )
    newY = Q(h0 + h1*p + h2*(p**2))
    liftedP = [integerX, newY]
    x2 = Q[0]
    y2 = Q[1]
    integerX = Q(Integer(x2))
    integerY = Integer(y2)
    h0 = integerY
    numH1 = ((y2 ** 2) - (x2 ** 3) - (a * x2) - b) / p
    h1 = int( - Field(numH1) / Field(2 * y2) )
    numH2 = (( (y2 + h1*p) ** 2) - (x2 ** 3) - a * x - b) / (p ** 2)
    h2 = int( - Field(numH2) / Field(2 * (y2 + h1*p)) )
    newY = Q(h0 + h1*p + h2*(p**2)) 
    liftedQ = [integerX, newY]
    newP = successiveSquare(p, liftedP, E, a)
    newQ = successiveSquare(p, liftedQ, E, a)
    psiP = - (newP[0] / newP[1])
    psiQ = - (newQ[0] / newQ[1])
    pval = Field( psiP.add_bigoh(2) / p )
    qval = Field( psiQ.add_bigoh(2) / p )
    log = qval / pval
   
def successiveSquare(m, P, E, a):
    m = Integer(m)
    if m==0: return "INF"
    elif m==1: return P
    elif (m % 2 == 0):
        return successiveSquare(m//2, double(P, E, a), E, a)
    else:
        return add(P, successiveSquare(m//2, double(P, E, a), E, a), E)

def double(P, E, a):
    if P == "INF": return P
    elif P == [P[0], -P[1]]: return "INF"
    x1 = P[0]
    y1 = P[1]
  
    N = ((3 * (x1 ** 2)) + a) / (2*y1)
    x2 = (N ** 2) - 2 * x1
    newY = N * (x1 - x2) - y1
    return [newx, newY]

def add(P, Q, E):
    if P == "INF": return Q
    elif Q == "INF": return P
    elif P == Q: return double(P, E)
    elif P == [Q[0], -Q[1]]: return "INF"
    xp = P[0]
    yp = P[1]
    xq = Q[0]
    yq = Q[1]
    M = (yq - yp) / (xq - xp)
    xr = (M ** 2) - xp - xq
    yr = M*(xp-xr) - yp
    return [xr, yr]
   

p = 730750818665451459112596905638433048232067471723
a = 425706413842211054102700238164133538302169176474
b = 203362936548826936673264444982866339953265530166
P = (24, 310224165475973298147806088269428225647703826034)
Q = (50, 66275076336461442314332875274548217082100991151)

E = EllipticCurve(GF(p), [a, b])
discrete_log(P, Q, E, p, a, b)