import sys
import numpy as np
import scipy
import scipy.sparse as ssp

import spams
import time
from test_utils import *

def test_fistaFlat():
    param = {'numThreads' : 1,'verbose' : False,
             'lambda1' : 0.05, 'it0' : 10, 'max_it' : 200,
             'L0' : 0.1, 'tol' : 1e-3, 'intercept' : False,
             'pos' : False}
    np.random.seed(0)
    m = 100;n = 200
    X = np.asfortranarray(np.random.normal(size = (m,n)))
    X = np.asfortranarray(X - np.tile(np.mean(X,0),(X.shape[0],1)))
    X = spams.normalize(X)
    Y = np.asfortranarray(np.random.normal(size = (m,1)))
    Y = np.asfortranarray(Y - np.tile(np.mean(Y,0),(Y.shape[0],1)))
    Y = spams.normalize(Y)
    W0 = np.zeros((X.shape[1],Y.shape[1]),dtype=np.float64,order="FORTRAN")
    # Regression experiments 
    # 100 regression problems with the same design matrix X.
    print '\nVarious regression experiments'
    param['compute_gram'] = True
    param['verbose'] = True
    print '\nFISTA + Regression l1'
    param['loss'] = 'square'
    param['regul'] = 'l1'
    # param.regul='group-lasso-l2';
    # param.size_group=10;
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,X,W0,True,**param)',locals())
    print 'mean loss: %f, mean relative duality_gap: %f, number of iterations: %f' %(np.mean(optim_info[0,:],0),np.mean(optim_info[2,:],0),np.mean(optim_info[3,:],0))
###
    print '\nISTA + Regression l1'
    param['ista'] = True
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,X,W0,True,**param)',locals())
    print 'mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n' %(np.mean(optim_info[0,:]),np.mean(optim_info[2,:]),np.mean(optim_info[3,:]))
##
    print '\nSubgradient Descent + Regression l1'
    param['ista'] = False
    param['subgrad'] = True
    param['a'] = 0.1
    param['b'] = 1000 # arbitrary parameters
    max_it = param['max_it']
    it0 = param['it0']
    param['max_it'] = 500
    param['it0'] = 50
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,X,W0,True,**param)',locals())
    print 'mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n' %(np.mean(optim_info[0,:]),np.mean(optim_info[2,:]),np.mean(optim_info[3,:]))
    param['subgrad'] = False
    param['max_it'] = max_it
    param['it0'] = it0

###
    print '\nFISTA + Regression l2'
    param['regul'] = 'l2'
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,X,W0,True,**param)',locals())
    print 'mean loss: %f, mean relative duality_gap: %f, number of iterations: %f\n' %(np.mean(optim_info[0,:]),np.mean(optim_info[2,:]),np.mean(optim_info[3,:]))
###
    print '\nFISTA + Regression l2 + sparse feature matrix'
    param['regul'] = 'l2';
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,ssp.csc_matrix(X),W0,True,**param)',locals())
    print 'mean loss: %f, mean relative duality_gap: %f, number of iterations: %f' %(np.mean(optim_info[0,:]),np.mean(optim_info[2,:]),np.mean(optim_info[3,:]))
###########

    print '\nFISTA + Regression Elastic-Net'
    param['regul'] = 'elastic-net'
    param['lambda2'] = 0.1
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,X,W0,True,**param)',locals())
    print 'mean loss: %f, number of iterations: %f' %(np.mean(optim_info[0,:]),np.mean(optim_info[3,:]))

    print '\nFISTA + Group Lasso L2'
    param['regul'] = 'group-lasso-l2'
    param['size_group'] = 2
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,X,W0,True,**param)',locals())
    print 'mean loss: %f, number of iterations: %f' %(np.mean(optim_info[0,:]),np.mean(optim_info[3,:]))
    
    print '\nFISTA + Trace Norm'
    param['regul'] = 'trace-norm-vec'
    param['size_group'] = 5
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,X,W0,True,**param)',locals())
    print 'mean loss: %f, number of iterations: %f' %(np.mean(optim_info[0,:]),np.mean(optim_info[3,:]))
    
####    
   
    print '\nFISTA + Regression Fused-Lasso'
    param['regul'] = 'fused-lasso'
    param['lambda2'] = 0.1
    param['lambda3'] = 0.1; #
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,X,W0,True,**param)',locals())
    print 'mean loss: %f, number of iterations: %f' %(np.mean(optim_info[0,:]),np.mean(optim_info[3,:]))
    
    print '\nFISTA + Regression no regularization'
    param['regul'] = 'none'
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,X,W0,True,**param)',locals())
    print 'mean loss: %f, number of iterations: %f' %(np.mean(optim_info[0,:]),np.mean(optim_info[3,:]))
    
    
    print '\nFISTA + Regression l1 with intercept '
    param['intercept'] = True
    param['regul'] = 'l1'
    x1 = np.asfortranarray(np.concatenate((X,np.ones((X.shape[0],1))),1))
    W01 = np.asfortranarray(np.concatenate((W0,np.zeros((1,W0.shape[1]))),0))
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,x1,W01,True,**param)',locals()) # adds a column of ones to X for the intercept,True)',locals())
    print 'mean loss: %f, mean relative duality_gap: %f, number of iterations: %f' %(np.mean(optim_info[0,:]),np.mean(optim_info[2,:]),np.mean(optim_info[3,:]))
    
    print '\nFISTA + Regression l1 with intercept+ non-negative '
    param['pos'] = True
    param['regul'] = 'l1'
    x1 = np.asfortranarray(np.concatenate((X,np.ones((X.shape[0],1))),1))
    W01 = np.asfortranarray(np.concatenate((W0,np.zeros((1,W0.shape[1]))),0))
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,x1,W01,True,**param)',locals())
    print 'mean loss: %f, number of iterations: %f' %(np.mean(optim_info[0,:]),np.mean(optim_info[3,:]))
    param['pos'] = False
    param['intercept'] = False

    print '\nISTA + Regression l0'
    param['regul'] = 'l0'
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,X,W0,True,**param)',locals())
    print 'mean loss: %f, number of iterations: %f' %(np.mean(optim_info[0,:]),np.mean(optim_info[3,:]))
    
# Classification
    
    print '\nOne classification experiment'
#    Y = 2 * double(randn(100,1) > 0)-1
    Y = np.asfortranarray(2 * np.asarray(np.random.normal(size = (100,1)) > 1,dtype='float64') - 1)
    print '\nFISTA + Logistic l1'
    param['regul'] = 'l1'
    param['loss'] = 'logistic'
    param['lambda1'] = 0.01
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,X,W0,True,**param)',locals())
    print 'mean loss: %f, mean relative duality_gap: %f, number of iterations: %f' %(np.mean(optim_info[0,:]),np.mean(optim_info[2,:]),np.mean(optim_info[3,:]))
# can be used of course with other regularization functions, intercept,...
    param['regul'] = 'l1'
    param['loss'] = 'weighted-logistic'
    param['lambda1'] = 0.01
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,X,W0,True,**param)',locals())
    print 'mean loss: %f, mean relative duality_gap: %f, number of iterations: %f' %(np.mean(optim_info[0,:]),np.mean(optim_info[2,:]),np.mean(optim_info[3,:]))
# can be used of course with other regularization functions, intercept,...
#!    pause
    
    

    print '\nFISTA + Logistic l1 + sparse matrix'
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,ssp.csc_matrix(X),W0,True,**param)',locals())
    print 'mean loss: %f, mean relative duality_gap: %f, number of iterations: %f' %(np.mean(optim_info[0,:]),np.mean(optim_info[2,:]),np.mean(optim_info[3,:]))
# can be used of course with other regularization functions, intercept,...
    

# Multi-Class classification
#    Y = double(ceil(5*rand(100,1000))-1)
    Y = np.asfortranarray(np.ceil(5 * np.random.random(size = (100,1000))) - 1)
    param['loss'] = 'multi-logistic'
    print '\nFISTA + Multi-Class Logistic l1'
    nclasses = np.max(Y[:])+1
    W0 = np.zeros((X.shape[1],nclasses * Y.shape[1]),dtype=np.float64,order="FORTRAN")
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,X,W0,True,**param)',locals())

    print 'mean loss: %f, mean relative duality_gap: %f, number of iterations: %f' %(np.mean(optim_info[0,:]),np.mean(optim_info[2,:]),np.mean(optim_info[3,:]))
# can be used of course with other regularization functions, intercept,...
    
    
# Multi-Task regression
    Y = np.asfortranarray(np.random.normal(size = (100,100)))
    Y = np.asfortranarray(Y - np.tile(np.mean(Y,0),(Y.shape[0],1)))
    Y = spams.normalize(Y)
    param['compute_gram'] = False
    param['verbose'] = True;   # verbosity, False by default
    W0 = np.zeros((X.shape[1],Y.shape[1]),dtype=np.float64,order="FORTRAN")
    param['loss'] = 'square'
    print '\nFISTA + Regression l1l2 '
    param['regul'] = 'l1l2'
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,X,W0,True,**param)',locals())
    print 'mean loss: %f, mean relative duality_gap: %f, number of iterations: %f' %(np.mean(optim_info[0,:]),np.mean(optim_info[2,:]),np.mean(optim_info[3,:]))
    
    print '\nFISTA + Regression l1linf '
    param['regul'] = 'l1linf'
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,X,W0,True,**param)',locals())
    print 'mean loss: %f, mean relative duality_gap: %f, number of iterations: %f' %(np.mean(optim_info[0,:]),np.mean(optim_info[2,:]),np.mean(optim_info[3,:]))
    
    
    print '\nFISTA + Regression l1l2 + l1 '
    param['regul'] = 'l1l2+l1'
    param['lambda2'] = 0.1
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,X,W0,True,**param)',locals())
    print 'mean loss: %f, number of iterations: %f' %(np.mean(optim_info[0,:]),np.mean(optim_info[3,:]))
    
    
    print '\nFISTA + Regression l1linf + l1 '
    param['regul'] = 'l1linf+l1'
    param['lambda2'] = 0.1
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,X,W0,True,**param)',locals())
    print 'mean loss: %f, number of iterations: %f' %(np.mean(optim_info[0,:]),np.mean(optim_info[3,:]))
    
    
    print '\nFISTA + Regression l1linf + row + columns '
    param['regul'] = 'l1linf-row-column'
    param['lambda2'] = 0.1
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,X,W0,True,**param)',locals())
    print 'mean loss: %f, mean relative duality_gap: %f, number of iterations: %f' %(np.mean(optim_info[0,:]),np.mean(optim_info[2,:]),np.mean(optim_info[3,:]))
    
# Multi-Task Classification
    
    print '\nFISTA + Logistic + l1l2 '
    param['regul'] = 'l1l2'
    param['loss'] = 'logistic'
#    Y = 2*double(randn(100,100) > 0)-1
    Y = np.asfortranarray(2 * np.asarray(np.random.normal(size = (100,100)) > 1,dtype='float64') - 1)
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,X,W0,True,**param)',locals())
    print 'mean loss: %f, mean relative duality_gap: %f, number of iterations: %f' %(np.mean(optim_info[0,:]),np.mean(optim_info[2,:]),np.mean(optim_info[3,:]))
# Multi-Class + Multi-Task Regularization
    
    
    param['verbose'] = False
    print '\nFISTA + Multi-Class Logistic l1l2 '
#    Y = double(ceil(5*rand(100,1000))-1)
    Y = np.asfortranarray(np.ceil(5 * np.random.random(size = (100,1000))) - 1)
    Y = spams.normalize(Y)
    param['loss'] = 'multi-logistic'
    param['regul'] = 'l1l2'
    nclasses = np.max(Y[:])+1
    W0 = np.zeros((X.shape[1],nclasses * Y.shape[1]),dtype=np.float64,order="FORTRAN")
    (W, optim_info) = Xtest1('spams','spams.fistaFlat(Y,X,W0,True,**param)',locals())
    print 'mean loss: %f, mean relative duality_gap: %f, number of iterations: %f' %(np.mean(optim_info[0,:]),np.mean(optim_info[2,:]),np.mean(optim_info[3,:]))
# can be used of course with other regularization functions, intercept,...
    
    
#############
def test_fistaGraph():
    return None

def test_fistaTree():
    return None

def test_proximalFlat():
    param = {'numThreads' : -1,'verbose' : True,
             'lambda1' : 0.1 }
    m = 100;n = 1000
    U = np.asfortranarray(np.random.normal(size = (m,n)))
    
    # test L0
    print "\nprox l0"
    param['regul'] = 'l0'
    param['pos'] = False       # false by default
    param['intercept'] = False # false by default
    alpha = Xtest1('spams','spams.proximalFlat(U,False,**param)',locals())

    # test L1
    print "\nprox l1, intercept, positivity constraint"
    param['regul'] = 'l1'
    param['pos'] = True       # can be used with all the other regularizations
    param['intercept'] = True # can be used with all the other regularizations
    alpha = Xtest1('spams','spams.proximalFlat(U,False,**param)',locals())

    # test L2
    print "\nprox squared-l2"
    param['regul'] = 'l2'
    param['pos'] = False
    param['intercept'] = False
    alpha = Xtest1('spams','spams.proximalFlat(U,False,**param)',locals())

# test elastic-net
    print "\nprox elastic-net"
    param['regul'] = 'elastic-net'
    param['lambda2'] = 0.1
    alpha = Xtest1('spams','spams.proximalFlat(U,**param)',locals())

# test fused-lasso
    print "\nprox fused lasso"
    param['regul'] = 'fused-lasso'
    param['lambda2'] = 0.1
    param['lambda3'] = 0.1
    alpha = Xtest1('spams','spams.proximalFlat(U,**param)',locals())

# test l1l2
    print "\nprox mixed norm l1/l2"
    param['regul'] = 'l1l2'
    alpha = Xtest1('spams','spams.proximalFlat(U,**param)',locals())

# test l1linf
    print "\nprox mixed norm l1/linf"
    param['regul'] = 'l1linf'
    alpha = Xtest1('spams','spams.proximalFlat(U,**param)',locals())

# test l1l2+l1
    print "\nprox mixed norm l1/l2 + l1"
    param['regul'] = 'l1l2+l1'
    param['lambda2'] = 0.1
    alpha = Xtest1('spams','spams.proximalFlat(U,**param)',locals())

# test l1linf+l1
    print "\nprox mixed norm l1/linf + l1"
    param['regul'] = 'l1linf+l1'
    param['lambda2'] = 0.1
    alpha = Xtest1('spams','spams.proximalFlat(U,**param)',locals())

# test l1linf-row-column
    print "\nprox mixed norm l1/linf on rows and columns"
    param['regul'] = 'l1linf-row-column'
    param['lambda2'] = 0.1
    alpha = Xtest1('spams','spams.proximalFlat(U,**param)',locals())

# test none
    print "\nprox no regularization"
    param['regul'] = 'none'
    alpha = Xtest1('spams','spams.proximalFlat(U,**param)',locals())


    return None

def test_proximalGraph():
    return None

def test_proximalTree():
    return None



tests = {
    'fistaFlat' : test_fistaFlat,
    'fistaGraph' : test_fistaGraph,
    'fistaTree' : test_fistaTree,
    'proximalFlat' : test_proximalFlat,
    'proximalGraph' : test_proximalGraph,
    'proximalTree' : test_proximalTree,
    }
