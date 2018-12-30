% pairwise distnace subroutine

function Dis = pwd_sub(q1, q2)

N = 3;
[K, M] = size(q1);
t = (0:M-1)/(M-1);

% use DP to find the optimal warping from q2 to q1
q2new = q2;
for i = 1:N
    [q2n,Rt] = Find_Best_Rotation(q1,q2new);
    
    norm1 = norm(q1);
    norm2 = norm(q2n);
    gam0 = DynamicProgrammingQ(q2n/norm2,q1/norm1,0,0);
    gam = (gam0-gam0(1))/(gam0(end)-gam0(1));    % normalization
    gam_dev = gradient(gam, 1/(M-1));
    q2new = interp1(t', q2n', (t(end)-t(1)).*gam' + t(1), 'linear', 'extrap')' ...
        .*(ones(K,1)*sqrt(gam_dev));
    clear dd
    for k = 1:K
        dd(k) = trapz(t, (q1(k,:)-q2new(k,:)).^2);
    end
    dis(i) = sqrt(sum(dd));
    
    if (i > 1) & (abs(dis(i) - dis(i-1))/dis(i-1) < 1e-2)
        Dis = dis(i);
        %disp(['early stop, i = ' num2str(i)]);
        break;
    end
end
Dis = dis(i);
% 
% figure;
% set(gcf, 'position', [300 300 1200 600]);
% for k = 1:K
%     subplot(K,1,k);
%     plot(t, q1(k,:), 'b', 'linewidth', 2);
%     hold on;
%     plot(t, q2(k,:), 'r', 'linewidth', 2);
%     plot(t, q2new(k,:), 'g', 'linewidth', 2);
% end

