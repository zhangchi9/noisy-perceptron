function capacity = numerical_capacity()
rin_range = [0.01:0.01:0.2];
a_range = [0.1,5:5:100];
% rin_range = [0.12];
% a_range = [35];
loadfac = [1,0.95,0.9,0.8,0.7];
n_capacity = nan(length(rin_range),length(a_range));
for i = 1:length(rin_range)
    %tmp = nan(1,length(a_range));
    for j = 1:length(a_range)
        [rin_range(i),a_range(j)]
        prob_learning = nan(5,2);
        for k = 1:5
            filename = ['N_400_rj_',num2str(rin_range(i)),'_a_',num2str(a_range(j)),'_loadfac_',num2str(loadfac(k)),'.mat'];
            if isfile(filename)
                file_variable = load(filename);
                prob_learning(k,1) = file_variable.m;
                prob_learning(k,2) = file_variable.prob_learning;
            end
        end
        n_capacity(i,j) = get_numerical_capacity(prob_learning);
    end
     %tmp;
end
save numerical_capacity_map.mat
%figure,imagesc(a_range,rin_range,n_capacity),axis xy
end

function capacity = get_numerical_capacity(prob_learning)
prob_learning = prob_learning(all(~isnan(prob_learning),2),:);
capacity = nan;
if max(prob_learning(:,2))>0.5 && min(prob_learning(:,2))<0.5
    ind = prob_learning(:,2)-0.5 > 0;
    [~,close_point(1)] = min(prob_learning(ind,2));
    [~,close_point(2)] = max(prob_learning(~ind,2));
    close_point(1) = close_point(1) + nnz(~ind);
    capacity = interp1([prob_learning(close_point(1),2),prob_learning(close_point(2),2)],[prob_learning(close_point(1),1),prob_learning(close_point(2),1)],0.5)/400;
end
end