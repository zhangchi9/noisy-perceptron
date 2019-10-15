clear
clc
rout = 0.1;
rou_range = [0.1,1:1:60];
fac_silent = zeros(1,length(rou_range));
for i=1:length(rou_range)
    rou = rou_range(i);
    [capacity,Pconexc] = Brunel_figure(rou,rout);
    fac_silent(i) = 1-Pconexc;
end
plot(rou_range*20/sqrt(150000),fac_silent)
ylim([0.5,1])
hold on 