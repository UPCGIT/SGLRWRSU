
function [U,SRE_A5,RMSE_5] = RSU_g(Y,M,opt)
%%
% Parameter settings

lambda_21 = opt.lambda_21;% Regularization parameter
lambda_g = opt.lambda_g
Lg = opt.Lg;
XT=opt.XT;
parw=opt.parw;
group= opt.group; 
mu=opt.mu ; % ADMM parameter

max_iters = opt.max_iters; % Maximum number of iterations
tolerance = 1e-6; % Convergence tolerance
%Y=V;

% initialization A, Z, U


IF=inv(M' * M + 3*eye(size(M, 2)));
U = IF*M'*Y;

%U = max(U,0);
V1 = M*U-Y;
D1= zeros(size(Y));
V2 = U;
D2 = zeros(size(U));
V3 = U;
D3 = zeros(size(U));
V4 = U;
D4 = zeros(size(U));
V5 = U;
D5 = zeros(size(U));
[L,p]=size(M);
[L,N] = size(Y);
[~, n] = size(M);
tol=1e-4;
iter=1;
 tol1 = sqrt(N)*tol;
 res=inf;

while (iter<max_iters)&&(res>tol1)
    % Update A
   
    U = IF*(M'*(V1+D1+Y)+V2+D2+V3+D3+V4+D4);
    V1=prox_group_lasso1(M*U-Y-D1,group,1/mu);
    V2  =prox_group_lasso1(U-D2,group,parw*lambda_21/mu);
    nu_aux_3=U-D3;
    for i=1:size(group,1)
        V3(:,group{i,1})=mu*nu_aux_3(:,group{i,1})*(2*lambda_g*Lg{1,i}+mu*eye(size(Lg{1,i})))^-1;
    end
       V4 = max(U-D4,0);
    D1 = D1 - (M*U-Y-V1);
    D2 = D2 - (U-V2);
    %D3 = D3 - (U-V3);
    D3 = D3 - (U-V3);
    D4 = D4 - (U-V4);
  %  D5 = D5 - (U-V5);
X_cal = U.*(U>=0);
X_cal(X_cal>1)=1;                               
r5=RMSE1(XT,X_cal);
SRE_A5(iter+1)=cal_SRE(XT,X_cal);
RMSE_5(iter+1)=RMSE1(XT,X_cal);
fprintf('Iteration %d: SRE = %f\n', iter, SRE_A5(iter+1));
iter=iter+1;
end
end

