clear all;
n=4000;
p=40000;
density=0.01;

% generate random data
format compact;
randn('seed',0);
rand('seed',0);
X=sprandn(p,n,density);
mean_nrm=mean(sqrt(sum(X.^2)));
X=X/mean_nrm;

% generate some true model
z=double(sign(full(sprandn(p,1,0.05))));  
y=X'*z;
nrm=sqrt(sum(y.^2));
y=y+0.01*nrm*randn(n,1);    % add noise to the model
nrm=sqrt(sum(y.^2));
y=y*(sqrt(n)/nrm);

% set optimization parameters
clear param;
param.regul='l1';        % many other regularization functions are available
param.loss='square';     % only square and log are available
param.num_threads=1;    % uses all possible cores
param.normalized=false;  % if the columns of X have unit norm, set to true.
param.strategy=3;
param.verbose=false;
param.minibatches=ceil(1/density);  

% set grid of lambda
max_lambda=max(abs(X*y))/n;
param.lambda=max_lambda*(2^(1/8)).^(0:-1:-50);    % best to order from large to small

tabepochs=[1 2 3 5];
%% The problem which will be solved is
%%   min_beta  1/(2n) ||y-X' beta||_2^2 + lambda ||beta||_1
% the problems for different lambdas are solve INDEPENDENTLY in parallel
fprintf('EXPERIMENT: ALL LAMBDAS IN PARALLEL\n');
param.warm_restart=false;
for ii=1:length(tabepochs)
   param.epochs=tabepochs(ii);   % one pass over the data
   fprintf('EXP WITH %d PASS\n',tabepochs(ii));
   nlambdas=length(param.lambda);
   Beta0=zeros(p,nlambdas);
   tic
   [Beta tmp]=mexIncrementalProx(y,X,Beta0,param);
   toc
   yR=repmat(y,[1 nlambdas]);
   fprintf('Objective functions: \n');
   0.5*sum((yR-X'*Beta).^2)/n+param.lambda.*sum(abs(Beta))
   fprintf('Sparsity: \n');
   sum(Beta ~= 0)
   pause
end

% the problems are here solved sequentially with warm restart
% this seems to be the prefered choice.
fprintf('EXPERIMENT: SEQUENTIAL LAMBDAS WITH WARM RESTART\n');
fprintf('A SINGLE CORE IS USED\n');
param.warm_restart=true;
param.num_threads=1;
for ii=1:length(tabepochs)
   param.epochs=tabepochs(ii);   % one pass over the data
   fprintf('EXP WITH %d PASS\n',tabepochs(ii));
   nlambdas=length(param.lambda);
   Beta0=zeros(p,nlambdas);
   tic
   [Beta tmp]=mexIncrementalProx(y,X,Beta0,param);
   toc
   yR=repmat(y,[1 nlambdas]);
   fprintf('Objective functions: \n');
   0.5*sum((yR-X'*Beta).^2)/n+param.lambda.*sum(abs(Beta))
   fprintf('Sparsity: \n');
   sum(Beta ~= 0)
   pause
end


