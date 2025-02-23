function Y = vector_soft_row(X,tau)
%
%  computes the vector soft columnwise


NU = sqrt(sum(X.^2,2));
A = max(0, NU-tau);
% Y=(A./(A+tau)).* X;
Y = repmat((A./(A+tau)),1,size(X,2)).* X;
% figure;
% for count1=1:size(GROUP,2)
%     ref1=Y(:,GROUP{count1});
%     subplot(6,7,count1)
%     plot(ref1);
% end