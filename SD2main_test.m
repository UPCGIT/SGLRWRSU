 
close all
clear all
clc

mkdir('examples/DC4')
% mkdir('mat_data')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of end members
p = 9;  % fixed for this demo

%SNR in dB
SNR =  20;
% noise bandwidth in pixels of the noise  low pass filter (Gaussian)
bandwidth = 10000; % 10000 == iid noise
%bandwidth = 5*pi/224; % colored noise

% define random states
rand('state',10);
randn('state',10);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load fractional abundances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load spatial2.mat

%  Size of the images
nl = size(Xim,1);
nc = size(Xim,2);
np = nl*nc;     % number of pixels



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% buid the dictionary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load USGS_1995_Library.mat
%  order bands by increasing wavelength
[dummy index] = sort(datalib(:,1));
A =  datalib(index,4:end);
names = names(4:end,:);

% prune the library
% min angle (in degres) between any two signatures
% the larger min_angle the easier is the sparse regression problem
min_angle = 4.44;
[A, index] = prune_library2(A,min_angle); % 240  signature
names = names(index',:);

% order  the columns of A by decreasing angles
[A, index, angles] = sort_library_by_angle(A);
names = names(index',:);
namesStr = char(names);

% Names of the first 10 ordered materials, with 4.44 deg. prunning:
% 1 - Jarosite GDS99 K,Sy 200C
% 2 - Jarosite GDS101 Na,Sy 200
% 3 - Anorthite HS349.3B
% 4 - Calcite WS272
% 5 - Alunite GDS83 Na63
% 6 - Howlite GDS155
% 7 - Corrensite CorWa-1
% 8 - Fassaite HS118.3B
% 9 - Adularia GDS57 Orthoclase
% 10 - Andradite NMNH113829


%% select p endmembers  from A

% angles (a_1,a_j) \sisizemeq min_angle)
% supp = 1:p;
supp = [2 3 4 5 6 7 8 9 10]; % dont take 2 Jarosites

% % Sample endmembers at random
% supp = randsample(size(A,2), p);

M = A(:,supp);
[L,p] = size(M);  % L = number of bands; p = number of material


%%
%---------------------------------
% generate  the observed  data X
%---------------------------------

% set noise standard deviation
sigma = sqrt(sum(sum((M*X).^2))/np/L/10^(SNR/10));
% generate Gaussian iid noise
noise = sigma*randn(L,np);


% make noise correlated by low pass filtering
% low pass filter (Gaussian)
filter_coef = exp(-(0:L-1).^2/2/bandwidth.^2)';
scale = sqrt(L/sum(filter_coef.^2));
filter_coef = scale*filter_coef;
noise = idct(dct(noise).*repmat(filter_coef,1,np));

%  observed spectral vector
Y = M*X + noise;


% create  true X wrt  the library A
n = size(A,2);
N = nl*nc;
XT = zeros(n,N);
XT(supp,:) = X;

%%
%Parameter settings

if SNR==20  
 opt_RSU_g.lambda_21=1;
 opt_RSU_g.lambda_g=0.1 ;
 opt_RSU_g.mu=1; % ADMM parameters
   slic_size  = 6;
slic_reg   = 0.004;
lambda = 0.01;
u=0.5;
lambda = 0.005;
end
%%
%Superpixel segmentation and construction of Laplacian matrix
X_3D=reshape(Y',100,100,224);
[group]=seg(X_3D,slic_size,slic_reg);
    for count1=1:size(group,1)        
          [row, col]=ind2sub([100,100],group{count1});
         Y1=Y(:,group{count1});
         Y2=[row;col];
        % Y2=Y2/max(max(Y2));
         n=size(Y1,2);
         sigma =0.8;  % Standard deviation in Gaussian kernel
         W=zeros(n,n);
        S1=zeros(n,n);
         for i=1:size(Y1,2)
               for j=1:size(Y1,2)
                   S1(i,j)=two_SAD(Y1(:,i),Y1(:,j));
               end
         end
        % S1=S1/max(max(S1));
         S2=pdist2(Y2', Y2', 'euclidean');
        % S2=S2/max(max(S2));
         dd1=S1;
         dd2=S2.*S2;
         W1=u*exp( -(dd1)/(sigma) ); % The distance between the i-th sample and the j-th sample
         W2 =(1-u)*exp( -(dd2)/(2*5.^2) );
          W=W1+W2;
         n = length(W);
         D = diag(sum(W, 2));
        Lg1{1,count1} = D-W;
   
    end 
%%
%%Prior weight calculation
Y_MID=zeros(size(Y));
for i=1:size(group,1)
     Y_MID(:,group{i})=repmat(sum(Y(:,group{i}),2)/size(group{i},2),1,size(group{i},2));
end
lambda = 0.01;
[U] =  sunsal(A,Y_MID,'lambda',lambda,'ADDONE','NO','POSITIVITY','yes', ...
                    'TOL',1e-4, 'AL_iters',2000,'verbose','yes');
           
       parw=zeros(size(U,1),1);
       for j=1:size(U,1)
               parw(j)= 1/(norm(U(j,:),2)+eps);
       end
RMSE_2=RMSE1(XT,U)
SRE_2=cal_SRE(XT,U)
%%
opt_RSU_g.XT = XT;
opt_RSU_g.max_iters=200;
opt_RSU_g.Lg=Lg1;
opt_RSU_g.parw=parw;
opt_RSU_g.group=group; 
[U5,DC2_20dB_SRE,DC2_20dB_RMSE] = RSU_g(Y,A,opt_RSU_g);
U5=U5.*(U5>0);
U5(U5>1)=1; 
RMSE_5=RMSE1(XT,U5)
SRE_A5=cal_SRE(XT,U5)
