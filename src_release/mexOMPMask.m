% 
% Usage:   A=mexOMPMask(X,D,B,param);
% or       [A path]=mexOMPMask(X,D,B,param);
%
% Name: mexOMPMask
%
% Description: mexOMPMask is a variant of mexOMP that allow using
%     a binary mask B
%     
%     for all columns x of X, and columns beta of B, it computes a column 
%         alpha of A by addressing
%         min_{alpha} ||alpha||_0  s.t  ||diag(beta)*(x-Dalpha)||_2^2 
%                                                               <= eps*||beta||_0/m
%         or
%         min_{alpha} ||diag(beta)*(x-Dalpha)||_2^2  s.t. ||alpha||_0 <= L
%         
%
% Inputs: X:  double m x n matrix   (input signals)
%            m is the signal size
%            n is the number of signals to decompose
%         D:  double m x p matrix   (dictionary)
%            p is the number of elements in the dictionary
%            All the columns of D should have unit-norm !
%         B:  boolean m x n matrix   (mask)
%               p is the number of elements in the dictionary
%         param: struct
%            param.L (maximum number of elements in each decomposition)
%            param.eps (threshold on the squared l2-norm of the residual
%            param.numThreads (optional, number of threads for exploiting
%            multi-core / multi-cpus. By default, it takes the value -1,
%            which automatically selects all the available CPUs/cores).
%
% Output: A: double sparse p x n matrix (output coefficients)
%         path (optional): double dense p x L matrix 
%                                     (regularization path of the first signal)
%
% Note: this function admits a few experimental usages, which have not
%     been extensively tested:
%      - single precision setting (even though the output alpha is double 
%        precision)
%      - Passing an int32 vector of length n to param.L allows to provide
%        a different parameter L for each input signal x_i
%      - Passing a double vector of length n to param.eps allows to provide
%        a different parameter eps for each input signal x_i
%
% Author: Julien Mairal, 2010


