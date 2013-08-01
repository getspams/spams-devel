n=200000;
p=20000;
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
param.num_threads=-1;    % uses all possible cores
param.normalized=false;  % if the columns of X have unit norm, set to true.
param.averaging_mode=0;
param.weighting_mode=0;  
param.verbose=false;

% set grid of lambda
max_lambda=max(abs(X*y))/n;
param.lambda=max_lambda*(2^(1/8)).^(0:-1:-25);    % best to order from large to small

tabepochs=[1 5 20 50];
%% The problem which will be solved is
%%   min_beta  1/(2n) ||y-X' beta||_2^2 + lambda ||beta||_1
fprintf('EXPERIMENT: ALL LAMBDAS IN PARALLEL\n');
% we try different experiments when varying the number of epochs.
% the problems for different lambdas are solve INDEPENDENTLY in parallel
for ii=1:length(tabepochs)
   param.iters=tabepochs(ii)*n;   % one pass over the data
   fprintf('EXP WITH %d PASS\n',tabepochs(ii));
   nlambdas=length(param.lambda);
   Beta0=zeros(p,nlambdas);
   tic
   [Beta tmp]=mexStochasticProx(y,X,Beta0,param);
   toc
   yR=repmat(y,[1 nlambdas]);
   fprintf('Objective functions: \n');
   0.5*sum((yR-X'*Beta).^2)/n+param.lambda.*sum(abs(Beta))
   fprintf('Sparsity: \n');
   sum(Beta ~= 0)
end


