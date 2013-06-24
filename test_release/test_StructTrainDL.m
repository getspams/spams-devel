clear all; 

I=double(imread('data/lena.png'))/255;
% extract 8 x 8 patches
X=im2col(I,[8 8],'sliding');
X=X-repmat(mean(X),[size(X,1) 1]);
X=X ./ repmat(sqrt(sum(X.^2)),[size(X,1) 1]);

param.K=64;  % learns a dictionary with 100 elements
param.lambda=0.05;
param.numThreads=4; % number of threads
param.batchsize=400;
param.tol = 1e-3

param.iter=20;  %

param.regul = 'l1';
fprintf('with Fista Regression %s\n',param.regul);
tic
D = mexStructTrainDL(X,param);
t=toc;
fprintf('time of computation for Dictionary Learning: %f\n',t);
param.approx=0;
fprintf('Evaluating cost function...\n');
alpha=mexLasso(X,D,param);
R=mean(0.5*sum((X-D*alpha).^2)+param.lambda*sum(abs(alpha)));
ImD=displayPatches(D);
subplot(1,3,1);
imagesc(ImD); colormap('gray');
fprintf('objective function: %f\n',R);
drawnow;
%
param.regul = 'l2';
fprintf('with Fista Regression %s\n',param.regul);
tic
D = mexStructTrainDL(X,param);
t=toc;
fprintf('time of computation for Dictionary Learning: %f\n',t);

param.approx=0;
fprintf('Evaluating cost function...\n');
alpha=mexLasso(X,D,param);
R=mean(0.5*sum((X-D*alpha).^2)+param.lambda*sum(abs(alpha)));
ImD=displayPatches(D);
subplot(1,3,2);
imagesc(ImD); colormap('gray');
fprintf('objective function: %f\n',R);
drawnow;

%
param.regul = 'elastic-net';
fprintf('with Fista %s\n',param.regul);
param.lambda2=0.1;
tic
D = mexStructTrainDL(X,param);
t=toc;
fprintf('time of computation for Dictionary Learning: %f\n',t);

param.approx=0;
fprintf('Evaluating cost function...\n');
alpha=mexLasso(X,D,param);
R=mean(0.5*sum((X-D*alpha).^2)+param.lambda*sum(abs(alpha)));
ImD=displayPatches(D);
subplot(1,3,3);
imagesc(ImD); colormap('gray');
fprintf('objective function: %f\n',R);
drawnow;
%%%%%%
%?pause
%%% GRAPH
param.lambda=0.1; % regularization parameter
param.tol=1e-5;
param.K = 10
graph.eta_g=[1 1 1 1 1];
graph.groups=sparse([0 0 0 1 0;
                     0 0 0 0 0;
                     0 0 0 0 0;
                     0 0 0 0 0;
                     0 0 1 0 0]);   % g5 is included in g3, and g2 is included in g4
graph.groups_var=sparse([1 0 0 0 0; 
                         1 0 0 0 0; 
                         1 0 0 0 0 ; 
                         1 1 0 0 0; 
                         0 1 0 1 0;
                         0 1 0 1 0;
                         0 1 0 0 1;
                         0 0 0 0 1;
                         0 0 0 0 1;
                         0 0 1 0 0]); % represents direct inclusion relations 

param.graph = graph

param.regul = 'graph';
fprintf('with Fista %s\n',param.regul);
tic
D = mexStructTrainDL(X,param);
t=toc;
fprintf('time of computation for Dictionary Learning: %f\n',t);

param.approx=0;
fprintf('Evaluating cost function...\n');
alpha=mexLasso(X,D,param);
R=mean(0.5*sum((X-D*alpha).^2)+param.lambda*sum(abs(alpha)));
ImD=displayPatches(D);
subplot(1,3,1);
imagesc(ImD); colormap('gray');
fprintf('objective function: %f\n',R);
drawnow;
%
param.regul = 'graph-ridge';
fprintf('with Fista %s\n',param.regul);
tic
D = mexStructTrainDL(X,param);
t=toc;
fprintf('time of computation for Dictionary Learning: %f\n',t);

param.approx=0;
fprintf('Evaluating cost function...\n');
alpha=mexLasso(X,D,param);
R=mean(0.5*sum((X-D*alpha).^2)+param.lambda*sum(abs(alpha)));
ImD=displayPatches(D);
subplot(1,3,2);
imagesc(ImD); colormap('gray');
fprintf('objective function: %f\n',R);
drawnow;
%

%%
%%% TREE
%?pause;
param = rmfield(param,'graph');

param.lambda=0.1; % regularization parameter
param.tol=1e-5;
param.K = 10

tree.own_variables=  int32([0 0 3 5 6 6 8 9]);   % pointer to the first variable of each group
tree.N_own_variables=int32([0 3 2 1 0 2 1 1]); % number of "root" variables in each group
tree.eta_g=[1 1 1 2 2 2 2.5 2.5];       
tree.groups=sparse([0 0 0 0 0 0 0 0; ...
                    1 0 0 0 0 0 0 0; ...
                    0 1 0 0 0 0 0 0; ...
                    0 1 0 0 0 0 0 0; ...
                    1 0 0 0 0 0 0 0; ...
                    0 0 0 0 1 0 0 0; ...
                    0 0 0 0 1 0 0 0; ...
                    0 0 0 0 0 0 1 0]);  % first group should always be the root of the tree

param.tree = tree

param.regul = 'tree-l0';
fprintf('with Fista %s\n',param.regul);

tic
D = mexStructTrainDL(X,param);
t=toc;
fprintf('time of computation for Dictionary Learning: %f\n',t);

param.approx=0;
fprintf('Evaluating cost function...\n');
alpha=mexLasso(X,D,param);
R=mean(0.5*sum((X-D*alpha).^2)+param.lambda*sum(abs(alpha)));
ImD=displayPatches(D);
subplot(1,3,1);
imagesc(ImD); colormap('gray');
fprintf('objective function: %f\n',R);
drawnow;
%
param.regul = 'tree-l2';
fprintf('with Fista %s\n',param.regul);

tic
D = mexStructTrainDL(X,param);
t=toc;
fprintf('time of computation for Dictionary Learning: %f\n',t);

param.approx=0;
fprintf('Evaluating cost function...\n');
alpha=mexLasso(X,D,param);
R=mean(0.5*sum((X-D*alpha).^2)+param.lambda*sum(abs(alpha)));
ImD=displayPatches(D);
subplot(1,3,3);
imagesc(ImD); colormap('gray');
fprintf('objective function: %f\n',R);
drawnow;

%
param.regul = 'tree-linf';
fprintf('with Fista %s\n',param.regul);

tic
D = mexStructTrainDL(X,param);
t=toc;
fprintf('time of computation for Dictionary Learning: %f\n',t);

param.approx=0;
fprintf('Evaluating cost function...\n');
alpha=mexLasso(X,D,param);
R=mean(0.5*sum((X-D*alpha).^2)+param.lambda*sum(abs(alpha)));
ImD=displayPatches(D);
subplot(1,3,3);
imagesc(ImD); colormap('gray');
fprintf('objective function: %f\n',R);
drawnow;

