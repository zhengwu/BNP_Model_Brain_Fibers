function [c,cn,thrd] = BtoCluster_zw(B,N,realCN)
% function for converting partition matrix B to cluster vector c and cn
c = zeros(1,N);
idx = 1:N;
cn = [];

max_elp = B(1,1);
Bstar = zeros(N,N);


thrd = max_elp;
tt =1;
while (thrd > 0)
    thrd = thrd-1;
    M = 1;
    c = zeros(1,N);
    idx = 1:N;
    cn = [];
    Btmp = B;

for i=1:N
    if(idx(i)>0)
       tmp_v = Btmp(idx(i),:);
       c(find(tmp_v>thrd)) = M;
       Btmp(:,find(tmp_v>thrd)) = 0;
       Btmp(find(tmp_v>thrd),:) = 0;
       idx(find(tmp_v>thrd)) = 0; % delete current clustering;
       cn(M) = length(find(tmp_v>thrd));
       M = M + 1;
    end;
end;
if(length(cn)==realCN)
    return 
end;

end;

disp('Can not find a thrd @_@ ... ');
cn = NaN;
c = NaN;
thrd = NaN;