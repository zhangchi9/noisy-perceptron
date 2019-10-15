function plot_current_distribution(W,p1)
close all
N = length(W);
f = 0.2;
p2 = (1-f)*p1/f;
X = rand(N,1)<f;
n = 3000;
x = nan(N,n);
%W = rand(N,1);
for j = 1 : n
    for i = 1 : N
        if X(i) == 0
            ranum = rand(1);
            if ranum < p1
                x(i,j) = 1;
            else
                x(i,j) = 0;
            end
        else
            ranum = rand(1);
            if ranum < p2
                x(i,j) = 0;
            else
                x(i,j) = 1;
            end
        end
    end
end
unper = W'*X/N-1;
per = W'*x/N-1;
normper = (per - mean(per)) / std(per);
[ks,p] =  kstest(normper)

h = histogram(per,'Normalization','pdf','EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5]);
bincenter = (h.BinEdges(1:(end-1)) + h.BinEdges(2:end))/2;
f = fit(bincenter',h.Values','gauss1');
hold on
plot(f,bincenter',h.Values')
%plot([unper,unper],[0,0.2],'r')
xlabel('Current')
ylabel('Probability')
axis square


