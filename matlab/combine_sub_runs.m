clear
clc
cd /scratch/zhang.chi9/perceptron/data/tmp_N_800/

model = {'perceptron','fmincon'};

network_number_range = 6:10;

W_all = nan(800,800);

for mod = model
    for i = network_number_range
        
        filename_all = sprintf('/scratch/zhang.chi9/perceptron/data/different_network_size_load_at_capacity/%s_rin_-6a_40TrialNum_%dN_800.mat',mod{1},i);
        for j = 0:7
            filename = [mod{1},'_',num2str(i),'_',num2str(j),'.mat'];
            s = load(filename);
            W_all(:,100*j + [1:100]) = s.W(:,100*j + [1:100]);
        end
        load(filename);
        W = W_all;
        save(filename_all)
    end
end