function [S,SS,L1,L2,final] = Chaosidx(W,I_ext,N)
% [S,SS,h,hh,L1,L2] = Chaosidx(W,N)
W = W';
Nruns = 500000;
S0 = rand(N,1)>rand;
tmp = randi(N);
S1 = S0;
S1(tmp) = -S0(tmp)+1; %flip
L1 = nan;
L2 = nan;
S(:,1) = S0;
SS(:,1) = S1;
for i = 1:Nruns-1
    idx1 = ((W*S(:,i)./N-1 + I_ext)>0);
    idx2 = ((W*SS(:,i)./N-1 + I_ext)>0);
    S(idx1,i+1) = 1;
    SS(idx2,i+1) = 1;
    
    if(mod(i,500)==0)
        idx = [];
        idx = find(sum(bsxfun(@minus,S(:,1:i),S(:,i+1))~=0)==0);
        if(~isempty(idx))
            L2 = i- idx(end)+1;
            L1 = idx(1);
            for j=1:499
                idxy = find(sum(bsxfun(@minus,S(:,1:i-j),S(:,i-j+1))~=0)==0);
                if(length(idxy)==length(idx))
                    L1=L1-1;
                else
                    break;
                end
            end
            if(L2==1)
                final = mean(S(:,i));
            else
                final = mean(mean(S(:,idx(end):i),2));
            end
            return;
        else
            L1=Nruns;
            L2=Nruns;
            final = nan;
        end
    end
    
end



end