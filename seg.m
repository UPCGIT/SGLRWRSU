function [Group]=seg(Y, slic_size, slic_reg)

% addpath 'F:\vlfeat-0.9.21\toolbox'
% run('F:\vlfeat-0.9.21\toolbox\vl_setup')

% % Extract features : SVD
% % apply SVD on img_noisy, then use first p components to compute similarity between two patches.
% r = im_size(1);
% c = im_size(2);
% 
% p   = 1;  
% [U_ss,D]=svd(Y,'econ');
% U_ss(:,p+1:end) = [];                         % remain first p component
% Y_feature =U_ss'*Y;                           % truncated SVD
% img_feature = reshape(Y_feature', r, c, p);
% 
% 
% [L, num] = superpixels(img_feature,num,'Compactness',10);

% reorder and rescale data into 2-D array
[numRows,numCols,numSpectra] = size(Y);
scfact = mean(reshape(sqrt(sum(Y.^2,3)), numRows*numCols, 1));
Y = Y./scfact;
imgVec = reshape(Y, [numRows*numCols numSpectra]);

% compute superpixels
disp('Computing SLIC Superpixels...');
L = vl_slic(single(Y), slic_size, slic_reg);
num = double(max(L(:)))+1; 

% BW = boundarymask(L);
% 
% % d=mat2gray(img_feature);
% % imshow(d);
% % imshow(imoverlay(d,BW,'cyan'),'InitialMagnification',67);

im_idx=reshape(1:numRows*numCols,numRows,numCols)';
im_idx=im_idx';
Group=cell(num, 1);
for i = 1:numRows
    for j = 1:numCols
        Group{L(i,j)+1,1}(end+1)=im_idx(i,j);
    end
end

end