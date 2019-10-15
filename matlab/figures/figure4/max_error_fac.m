function max_error_fac()
load prob_retrieval_map_load85_same_noise.mat

rin_range = [-12:0.5:-2];
% a_range = [0.1,5:5:100];
%rin_range = [-9.5];
a_range = [0.1,5:5:100];
th = 0.5;
max_err_map = nan(length(rin_range),length(a_range));
max_err_fac = nan(length(rin_range),length(a_range));

for i = 1:length(rin_range)
	tmp = nan(1,length(a_range));
    tmp2 = nan(1,length(a_range));
    for j = 1:length(a_range)
        [rin_range(i),a_range(j)]
        retrieval_prob = nan(20,2);
        retrieval_prob(:,1) = [5:5:100];
        retrieval_prob(:,2) = [noise5(i,j),noise10(i,j),noise15(i,j),...
            noise20(i,j),noise25(i,j),noise30(i,j),noise35(i,j),noise40(i,j),...
            noise45(i,j),noise50(i,j),noise55(i,j),noise60(i,j),noise65(i,j),...
            noise70(i,j),noise75(i,j),noise80(i,j),noise85(i,j),noise90(i,j),noise95(i,j),noise100(i,j)];
        tmp(j) = find_fac(retrieval_prob,th);
        tmp2(j) = find_fac(retrieval_prob,th)/a_range(j);
    end
         max_err_map(i,:) = tmp;
         max_err_fac(i,:) = tmp2;
end
save max_error.mat
end

function capacity = find_fac(retrieval_prob,th)

retrieval_prob = retrieval_prob(all(~isnan(retrieval_prob),2),:);
capacity = nan;
if max(retrieval_prob(:,2))>th && min(retrieval_prob(:,2))<th
    ind = retrieval_prob(:,2)-th > 0;
    close_point = [nnz(ind),nnz(ind)+1];
    capacity = interp1([retrieval_prob(close_point(1),2),retrieval_prob(close_point(2),2)],[retrieval_prob(close_point(1),1),retrieval_prob(close_point(2),1)],th);
else
    1
end


end


