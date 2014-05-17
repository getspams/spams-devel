clear all; 

I=double(imread('data/lena.png'))/255;
% extract 8 x 8 patches
X=im2col(I,[8 8],'sliding');
X=X-repmat(mean(X),[size(X,1) 1]);
X=X ./ repmat(sqrt(sum(X.^2)),[size(X,1) 1]);
X=X(:,1:1000); 

param.p=64;  % learns a dictionary with 64 elements
param.robust=true;
param.epsilon=1e-3;  % width for Huber loss
param.computeXtX=false;
param.stepsFISTA=0;
param.stepsAS=10;

%%%%%%%%%% FIRST EXPERIMENT %%%%%%%%%%%
tic
Z = mexArchetypalAnalysis(X,param);
t=toc;
fprintf('time of computation for Dictionary Learning: %f\n',t);

fprintf('Evaluating cost function...\n');
alpha=mexDecompSimplex(X,Z);
if param.robust
  R=sum(sqrt(sum((X-Z*alpha).^2)));
else
  R=sqrt(sum(sum((X-Z*alpha).^2)));
end
fprintf('objective function: %f\n',R);

%%%%%%%%%% FIRST EXPERIMENT %%%%%%%%%%%
tic
Z2 = mexArchetypalAnalysisContinue(X,Z,param);
t=toc;
fprintf('time of computation for Dictionary Learning: %f\n',t);

fprintf('Evaluating cost function...\n');
alpha=mexDecompSimplex(X,Z2);
if param.robust
  R=sum(sqrt(sum((X-Z2*alpha).^2)));
else
  R=sqrt(sum(sum((X-Z2*alpha).^2)));
end
fprintf('objective function: %f\n',R);
