function TT=two_SAD(spect_temp,Class_spectral)

for s=1:size(spect_temp,2)
    T=spect_temp(:,s);
    for t=1:size(Class_spectral,2)
        temp_spectralbundle=Class_spectral(:,t);
        if isequal(T, temp_spectralbundle)
            TT(s,t)=0;
        else
            for i=1:size(T,2)
                for j=1:size(temp_spectralbundle,2)
                    temp_s_sad(i,j)=acos(dot(T(:,i),temp_spectralbundle(:,j))/norm(T(:,i),2)/norm(temp_spectralbundle(:,j),2));
                end
            end
            TT(s,t)=min(temp_s_sad(:));
            clear temp_s_sad
        end
    end
end