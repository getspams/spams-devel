import sys
import numpy as np
import scipy
import scipy.sparse

import spams
import time
from test_utils import *

if not ('rand' in scipy.sparse.__dict__):
    import myscipy_rand as ssp
else:
    import scipy.sparse as ssp

def test_sparseProject():
    np.random.seed(0)
    X = np.asfortranarray(np.random.normal(size = (20000,100)))
    # matlab : X=X./repmat(sqrt(sum(X.^2)),[size(X,1) 1]);
    X = np.asfortranarray(X / np.tile(np.sqrt((X*X).sum(axis=0)),(X.shape[0],1)))
    param = {'numThreads' : -1, # number of processors/cores to use (-1 => all cores)
             'pos' : False,
             'mode': 1, # projection on the l1 ball
             'thrs' : 2}
    print "\n  Projection on the l1 ball"
    tic = time.time()
    X1 = spams.sparseProject(X,**param)
    tac = time.time()
    t = tac - tic
    print "  Time : ", t
    if (t != 0):
        print "%f signals of size %d projected per second" %((X.shape[1] / t),X.shape[0])
    s = np.abs(X1).sum(axis=0)
    print "Checking constraint: %f, %f" %(min(s),max(s))

    print "\n  Projection on the Elastic-Net"
    param['mode'] = 2  # projection on the Elastic-Net
    param['lambda1'] = 0.15
    tic = time.time()
    X1 = spams.sparseProject(X,**param)
    tac = time.time()
    t = tac - tic
    print "  Time : ", t
    if (t != 0):
        print "%f signals of size %d projected per second" %((X.shape[1] / t),X.shape[0])
    constraints = (X1*X1).sum(axis=0) + param['lambda1'] * np.abs(X1).sum(axis=0)
    print 'Checking constraint: %f, %f (Projection is approximate : stops at a kink)' %(min(constraints),max(constraints))
    
    print "\n  Projection on the FLSA"
    param['mode'] = 6       # projection on the FLSA
    param['lambda1'] = 0.7
    param['lambda2'] = 0.7
    param['lambda3'] = 1.0
    X = np.asfortranarray(np.random.random(size = (2000,100)))
    # matlab : X=X./repmat(sqrt(sum(X.^2)),[size(X,1) 1]);
    X = np.asfortranarray(X / np.tile(np.sqrt((X*X).sum(axis=0)),(X.shape[0],1)))
    tic = time.time()
    X1 = spams.sparseProject(X,**param)
    tac = time.time()
    t = tac - tic
    print "  Time : ", t
    if (t != 0):
        print "%f signals of size %d projected per second" %((X.shape[1] / t),X.shape[0])
    constraints = 0.5 * param['lambda3'] * (X1*X1).sum(axis=0) + param['lambda1'] * np.abs(X1).sum(axis=0) + \
    param['lambda2'] * np.abs(X1[2:,] - X1[1:-1,]).sum(axis=0)
    print 'Checking constraint: %f, %f (Projection is approximate : stops at a kink)' %(min(constraints),max(constraints))
    return None

def test_cd():
    return None

def test_L1L2BCD():
    return None

def test_lasso():
    np.random.seed(0)
    print "test Lasso"
##############################################
# Decomposition of a large number of signals
##############################################
# data generation
    X = np.asfortranarray(np.random.normal(size=(100,100000)))
    # X=X./repmat(sqrt(sum(X.^2)),[size(X,1) 1]);
    X = np.asfortranarray(X / np.tile(np.sqrt((X*X).sum(axis=0)),(X.shape[0],1)))
    D = np.asfortranarray(np.random.normal(size=(100,200)))
    D = np.asfortranarray(D / np.tile(np.sqrt((D*D).sum(axis=0)),(D.shape[0],1)))
    # parameter of the optimization procedure are chosen
#param.L=20; # not more than 20 non-zeros coefficients (default: min(size(D,1),size(D,2)))
    param = {
        'lambda1' : 0.15, # not more than 20 non-zeros coefficients
        'numThreads' : -1, # number of processors/cores to use; the default choice is -1
        # and uses all the cores of the machine
        'mode' : 2}        # penalized formulation

    tic = time.time()
    alpha = spams.lasso(X,D = D,return_reg_path = False,**param)
    tac = time.time()
    t = tac - tic
    print "%f signals processed per second\n" %(float(X.shape[1]) / t)
########################################
# Regularization path of a single signal 
########################################
    X = np.asfortranarray(np.random.normal(size=(64,1)))
    D = np.asfortranarray(np.random.normal(size=(64,10)))
    D = np.asfortranarray(D / np.tile(np.sqrt((D*D).sum(axis=0)),(D.shape[0],1)))
    (alpha,path) = spams.lasso(X,D = D,return_reg_path = True,**param)
    return None

def test_lassoMask():
    return None

def test_lassoWeighted():
    return None

def test_omp():
    return None

def test_ompMask():
    return None

def test_somp():
    return None


tests = {
    'sparseProject' : test_sparseProject,
    'cd' : test_cd,
    'L1L2BCD' : test_L1L2BCD,
    'lasso' : test_lasso,
    'lassoMask' : test_lassoMask,
    'lassoWeighted' : test_lassoWeighted,
    'omp' : test_omp,
    'ompMask' : test_ompMask,
    'somp' : test_somp,
    }
