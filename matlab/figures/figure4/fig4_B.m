clear
clc
load('rin_-6a_40TrialNum_1loadfac_0.85.mat')

if_retrieval = nan(1,100);
figure
for ind = 1:10
    hd = zeros(1,size(X,2)-1);
    for i = 1 : (size(X,2)-1)
        if i == 1
            x = X(:,i);
            %                     p1 = rj/2/(1-f);
            %                     p2 = rj/2/f;
            %                     err1 = rand(N,1)<p1;    % 0 -> 1
            %                     err2 = rand(N,1)<p2;    % 1 -> 0
            %                     x(x==0) = x(x==0) + err1(x==0);
            %                     x(x==1) = x(x==1) - err2(x==1);
            %                     raster(:,i) = x;
        else
            x = y;
        end
        y = W'*x/N + randn(N,1)*38/sqrt(N) > 1;
        hd(i) = nnz(X(:,i+1)-y)/N;
    end
    plot(hd),hold on
end
xlim([0,64])
axis square 
plot([0 64],[0.15,0.15],'--')