clear
clc
load motif_data_2neuron.mat
close all
for i = 1:5
    for j = 1:10
        subplot(5,10,10*(i-1)+j),errorbar(zscore_mean_map{i,j},zscore_std_map{i,j})
        axis square
        grid on 
        set(gca,'xtick',[],'ytick',[])
        
    end
end
drawnow
