% This function trains a perceptron using the association that one gives.The training algorithm is hebbian learning rule
% X is the input matrix with dimension N*m
% y is the out matrix with dimension a*m
% runmax is the max running steps
% If this association can be learned, then suc = 1, otherwise suc = 0
% runsteps is the total run steps for this training. 
% w is the conection weights with dimension N*a
% the ith coloum is the projection from N neuron to neuron i

function [suc,runsteps,w]=hebbian(X,Y,runmax)

[N,m] = size(X);
[a,b] = size(Y);
J = zeros(N,a);

if b ~= m
    error('The column dimension of X and Y does not equal');
end

for i = 1:a
    
    y = Y(i,:);
    j = J(:,a);
    suc = 0 ;
    runs=0;
    Learned=((j'*X).*y)>0;
    
    while runs<runmax && nnz(Learned)<m
        mu=find(Learned==0,1,'first');
        dj=X(:,mu).*y(mu)./N^0.5;
        j = j + dj;
        runs=runs+1;
        Learned=((j'*X).*y)>0;
    end
    
    J(:,i) = j;
    if nnz(Learned)==m
        suc = 1;
    else
        break;
    end
    
end

if nargout > 1
    runsteps = runs;
    if nargout > 2
        w = J;
    end
end