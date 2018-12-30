function B = calculate_B(c,N)

Ncluster = unique(c);
B = zeros(N,N);

for i=1:N
   B(i,i) = 1; 
end

for i=1:length(Ncluster)
    curr_idx = find(c == Ncluster(i));
%     for j=1:length(curr_idx)
%     B(curr_idx(j),curr_idx) = 1;
%     end
    B(curr_idx,curr_idx) = 1;
end

%  for row =1:N
%         for col = 1:N
%             if(c(row)==c(col))
%                 A(row,col) = 1;
%             end
%         end
%  end
%  
%  keyboard;