clear
clc
for loadfac = 0.5:0.1:1.5
    tmp = [];
    flagtmp = [];
    for num = 1:1:5
        filename = ['N_800_rj_-6_beta_',num2str(40/14),'_loadfac_',num2str(loadfac),'tri_num',num2str(num),'.mat'];
        if isfile(filename)
            load([filename]);
            tmp = [tmp,epsilonsum];
            flagtmp = [flagtmp,exitflag];
        end
    end
    epsilonsum = tmp;
    exitflag = flagtmp;
    save(['N_800_rj_-6_a_',num2str(40),'_loadfac_',num2str(loadfac),'.mat']);
end