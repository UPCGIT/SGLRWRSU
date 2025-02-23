function SRE=cal_SRE(XT,X)
SRE=20*log10(norm(XT,'fro')/norm(X-XT,'fro'));
