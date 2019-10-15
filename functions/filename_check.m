function filename = filename_check(savepath,rj_ind,a,TrialNum,N)
filename = [savepath,'/rin_',num2str(rj_ind),'a_',num2str(a),'TrialNum_' num2str(TrialNum),'N_',num2str(N),'.mat'];
while isfile(filename)
    TrialNum = TrialNum + 1;
    if TrialNum > 100
        error('100 network has generated')
    end
    filename = ['/home/zhang.chi9/research/logscale/network_different_size/rin_',num2str(rj_ind),'a_',num2str(a),'TrialNum_' num2str(TrialNum),'N_',num2str(N),'.mat'];
end
end