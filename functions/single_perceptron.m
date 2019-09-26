
clear
clc
N = 50;
fin = 0.1;
fout = 0.1;
total_training = 200;
suc = zeros(1,total_training);
mm = 50:10:100;
rate = zeros(1,length(mm));
for i = 1:length(mm)
    m = mm(i);
    
    for j = 1: total_training
        [m,j]
        
        X = 2*(rand(N,m)<fin)-1;
        Y = 2*(rand(1,m)<fout)-1;
        %Y = -Y;
        [suc(j)]=convex(X,Y);
        %[suc(j)]=hebbian(X,Y,500000);
    end
    suc(find(suc==-1))=[];
    rate(i) = sum(suc)/(length(suc));
    
end
plot(mm/N,rate,'*-')



