clear all;
n=40000;
p=100;
density=0.01;
density=1;

% generate random data
format compact;
randn('seed',0);
rand('seed',0);
X=sprandn(p,n,density);
mean_nrm=mean(sqrt(sum(X.^2)));
X=X/mean_nrm;
X=full(X);

% generate some true model
z=double(sign(full(sprandn(p,1,0.05))));  
y=X'*z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPERIMENT 1: Lasso
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('EXPERIMENT FOR LASSO\n');
nrm=sqrt(sum(y.^2));
y=y+0.01*nrm*randn(n,1);    % add noise to the model
nrm=sqrt(sum(y.^2));
y=y*(sqrt(n)/nrm);

% set optimization parameters
clear param;
param.regul='log-dc';        % many other regularization functions are available
param.eps=0.001;
param.loss='square';     % only square and log are available
param.numThreads=-1;    % uses all possible cores
param.normalized=false;  % if the columns of X have unit norm, set to true.
param.strategy=3;        % MISO with all heuristics
                         % 0: no heuristics, slow  (only for comparison purposes)
                         % 1: adjust the constant L on 5% of the data 
                         % 2: adjust the constant L on 5% of the data + unstable heuristics (this strategy does not work)
                         % 3: adjust the constant L on 5% of the data + stable heuristic (this is by far the best choice)
param.verbose=true;
param.minibatches=min(n,ceil(1/density));  % size of the minibatches, requires to store twice the size of X 

% set grid of lambda
max_lambda=max(abs(X*y))/n;
tablambda=max_lambda*(2^(1/8)).^(0:-1:-50);  % order from large to small
tablambda=tablambda/50;
tabepochs=[1 2 3 5 10];  % in this script, we compare the results obtained when changing the number of passes on the data.

param.lambda=tablambda;
%%%% EXP with ISTA
nlambdas=length(param.lambda);
obj=zeros(length(tabepochs),nlambdas);
spar=[];
for ii=1:length(tabepochs)
   param.max_iter=tabepochs(ii);   % one pass over the data
   fprintf('EXP WITH %d PASS OVER THE DATA\n',tabepochs(ii));
   Beta0=zeros(p,nlambdas);
   paramfista.regul=param.regul;
   paramfista.loss=param.loss;
   paramfista.max_it=tabepochs(ii);
   paramfista.tol=1e-10;
   paramfista.verbose=false;
   paramfista.a=param.eps;
   paramfista.ista=true;
   paramfista.it0=1;
   for jj=1:nlambdas
      paramfista.lambda=n*param.lambda(jj);
      [Beta(:,jj) tmp]=mexFistaFlat(y,X',Beta0(:,jj),paramfista);
      obj(ii,jj)=tmp(1)/n;
   end
   fprintf('Objective functions: \n');
   obj
   fprintf('Sparsity: \n');
   spar=[spar; sum(Beta ~= 0)];
   spar
end
objista=obj;
sparista=spar;
BetaIsta=Beta;


%%% exp 1: batch algorithm 
%% The problem which will be solved is
%fprintf('EXPERIMENT: ALL LAMBDAS WITHOUT WARM RESTART\n');
%param.minibatches=n;
%param.strategy=0;
%param.warm_restart=false;
%param.verbose=false;
%obj=[];
%spar=[];
%for ii=1:length(tabepochs)
%   param.epochs=tabepochs(ii);   % one pass over the data
%   fprintf('EXP WITH %d PASS OVER THE DATA\n',tabepochs(ii));
%   nlambdas=length(param.lambda);
%   Beta0=zeros(p,nlambdas);
%   tic
%   [Beta tmp]=mexIncrementalProx(y,X,Beta0,param);
%   toc
%   fprintf('Objective functions: \n');
%   obj=[obj; tmp(1,:)];
%   obj
%   fprintf('Sparsity: \n');
%   spar=[spar; sum(Beta ~= 0)];
%   spar
%end
%obj3=obj;


%% The problem which will be solved is
param.minibatches=min(n,ceil(1/density));  % size of the minibatches, requires to store twice the size of X 
param.strategy=3;
fprintf('EXPERIMENT: ALL LAMBDAS WITHOUT WARM RESTART\n');
param.warm_restart=false;
param.verbose=false;
obj=[];
spar=[];
for ii=1:length(tabepochs)
   param.epochs=tabepochs(ii);   % one pass over the data
   fprintf('EXP WITH %d PASS OVER THE DATA\n',tabepochs(ii));
   nlambdas=length(param.lambda);
   Beta0=zeros(p,nlambdas);
   tic
   [Beta tmp]=mexIncrementalProx(y,X,Beta0,param);
   toc
   fprintf('Objective functions: \n');
   obj=[obj; tmp(1,:)];
   obj
   fprintf('Sparsity: \n');
   spar=[spar; sum(Beta ~= 0)];
   spar
end
obj2=obj;
spar2=spar;
Beta2=Beta;

% the problems are here solved sequentially with warm restart
% this seems to be the prefered choice.
fprintf('EXPERIMENT: SEQUENTIAL LAMBDAS WITH WARM RESTART\n');
fprintf('A SINGLE CORE IS USED\n');
param.warm_restart=true;
param.num_threads=1;
obj=[];
spar=[];
for ii=1:length(tabepochs)
   param.epochs=tabepochs(ii);   % one pass over the data
   fprintf('EXP WITH %d PASS\n',tabepochs(ii));
   nlambdas=length(param.lambda);
   Beta0=zeros(p,nlambdas);
   tic
   [Beta tmp]=mexIncrementalProx(y,X,Beta0,param);
   toc
   fprintf('Objective functions: \n');
   obj=[obj; tmp(1,:)];
   obj
   fprintf('Sparsity: \n');
   spar=[spar; sum(Beta ~= 0)];
   spar
end


