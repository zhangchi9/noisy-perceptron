function count = motif2struct_bin(WE)
WW = WE.*WE';
count = zeros(2,1);
n_uc = sum(sum(WE))-sum(WW(:)); %number unidirectional connections
n_bc = (sum(WW(:)) - sum(diag(WW)))/2;%number of bidirectional connections
count(1,1) = n_uc;
count(2,1) = n_bc;
end
