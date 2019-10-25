function filename = filename_check(savepath,model,rj_ind,a,TrialNum,N)
    filename = [savepath,'/',model,'_rin_',num2str(rj_ind),'a_',num2str(a),'TrialNum_' num2str(TrialNum),'N_',num2str(N),'.mat'];
while isfile(filename)
    TrialNum = TrialNum + 1;
    if TrialNum > 100
        error('100 network has generated')
    end
    filename = [savepath,'/',model,'_rin_',num2str(rj_ind),'a_',num2str(a),'TrialNum_' num2str(TrialNum),'N_',num2str(N),'.mat'];
end
end