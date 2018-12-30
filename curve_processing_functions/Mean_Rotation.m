function Omean = Mean_Rotation(O)

% compute the Karcher mean of a set of rotational matrices

epsilon = 1;

N = size(O,1);
M = size(O,3);

Omean = eye(N);

for iter = 1:10
    
    S = zeros(N); 
    for i=1:M
        A(:,:,i) = real(logm(Omean'*O(:,:,i)));
        S = S + A(:,:,i);
    end
    A_bar  = S/M;
    
    Omean = Omean*expm(epsilon*A_bar);   
    if norm(A_bar,'fro') < 1e-10
        break;
    end
end
return; 