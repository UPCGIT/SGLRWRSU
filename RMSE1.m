function rmse=RMSE(Y,Y1)

[band,N]=size(Y);
rmse=0;
for i=1:N
   rmse=rmse+norm(Y(:,i)-Y1(:,i),2)^2;
    
end
rmse=sqrt(rmse/(band*N));