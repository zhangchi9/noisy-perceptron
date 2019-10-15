function capacity = get_numerical_capacity(prob_learning)
prob_learning = prob_learning(all(~isnan(prob_learning),2),:);
capacity = nan;
if max(prob_learning(:,2))>=0.5 && min(prob_learning(:,2))<=0.5
    ind = prob_learning(:,2)-0.5 >= 0;
    [~,close_point(1)] = max(prob_learning(ind,1));
    [~,close_point(2)] = min(prob_learning(~ind,1));
    close_point(2) = close_point(2) + nnz(ind);
    capacity = interp1([prob_learning(close_point(1),2),prob_learning(close_point(2),2)],[prob_learning(close_point(1),1),prob_learning(close_point(2),1)],0.5);
end
end