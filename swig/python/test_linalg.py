import sys
import numpy as np
import scipy
import scipy.sparse as ssp
import spams
import time
from test_utils import *

if not ('rand' in ssp.__dict__):
    import myscipy_rand
    ssprand = myscipy_rand.rand
else:
    ssprand = ssp.rand

def test_Sort():
    n = 2000000
    X = np.random.normal(size = (n,))
    return Xtest('np.sort(X)','spams.Sort(X,True)',locals())


def test_CalcAAt():
    """
    test A * A'
    """
    m=200; n = 200000; d= 0.05
    A = ssprand(m,n,density=d,format='csc',dtype=np.float64)
    return Xtest('np.dot(A,A.T)','spams.CalcAAt(A)',locals())

def test_CalcXAt():
    m=200; n = 200000; d= 0.05
    A = ssprand(m,n,density=d,format='csc',dtype=np.float64)
    X = np.asfortranarray(np.random.normal(size = (64,n)))

    # dot is very very slow betewwen a full and a sparse matrix
    return Xtest('np.dot(X,A.T.todense())','spams.CalcXAt(X,A)',locals())

def test_CalcXY():
    X = np.asfortranarray(np.random.normal(size = (64,200)))
    Y = np.asfortranarray(np.random.normal(size = (200,20000)))
    return Xtest('np.dot(X,Y)','spams.CalcXY(X,Y)',locals())

def test_CalcXYt():
    X = np.asfortranarray(np.random.normal(size = (64,200)))
    Y = np.asfortranarray(np.random.normal(size = (20000,200)))
    return Xtest('np.dot(X,Y.T)','spams.CalcXYt(X,Y)',locals())

def test_CalcXtY():
    X = np.asfortranarray(np.random.normal(size = (200,64)))
    Y = np.asfortranarray(np.random.normal(size = (200,20000)))
    return Xtest('np.dot(X.T,Y)','spams.CalcXtY(X,Y)',locals())

def test_Bayer():
    n = 2000000
    X = np.random.normal(size = (n,))

    Z = Xtest1('spams','spams.Bayer(X,0)',locals())
    return None

def test_ConjGrad():
    A = np.asfortranarray(np.random.normal(size = (5000,500)))
#    np.random.seed(0)
#    A = np.asfortranarray(np.random.normal(size = (10,5)))
    A = np.asfortranarray(np.dot(A.T,A))
    b = np.ones((A.shape[1],),dtype=np.float64,order="FORTRAN")
    x0 = b
    tol = 1e-4
    itermax = int(0.5 * len(b))

    tic = time.time()
    for i in xrange(0,20):
        y1 = np.linalg.solve(A,b)
    tac = time.time()
    print "  Time (numpy): ", tac - tic
    x1 = np.abs(b - np.dot(A,y1))
    print "Mean error on b : %f" %(x1.sum() / b.shape[0])

    tic = time.time()
    for i in xrange(0,20):
        y2 = spams.ConjGrad(A,b,x0,tol,itermax)
#        y2 = spams.ConjGrad(A,b)
    tac = time.time()
    print "  Time (spams): ", tac - tic
    x1 = np.dot(A,y2)
    x2 = np.abs(b - x1)
    print "Mean error on b : %f" %(x2.sum() / b.shape[0])

    err = abs(y1 - y2)
    return err.max()

def test_InvSym():
    A = np.asfortranarray(np.random.random(size = (1000,1000)))
    A =np.asfortranarray( np.dot(A.T,A))
    return Xtest('np.linalg.inv(A)','spams.InvSym(A)',locals())

def test_Normalize():
    A = np.asfortranarray(np.random.random(size = (100,1000)))
    res2 = Xtest1('spams','spams.Normalize(A)',locals())
    return None

tests = {
    'Sort' : test_Sort,
    'CalcAAt' : test_CalcAAt,
    'CalcXAt' : test_CalcXAt,
    'CalcXY' : test_CalcXY,
    'CalcXYt' : test_CalcXYt,
    'CalcXtY' : test_CalcXtY,
    'Bayer' : test_Bayer,
    'ConjGrad' : test_ConjGrad,
    'InvSym' : test_InvSym,
    'Normalize' : test_Normalize,
    }
