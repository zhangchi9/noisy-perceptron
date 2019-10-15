clear
clc
close all
if ~ispc
    load '/home/chi/Dropbox/Project_Perceptron/codes/figures/figure2/fig2data.mat'
else
    load 'C:\Users\chizhang\Dropbox\Project_Perceptron\codes\figures\figure2\fig2data.mat'
end
figure,imagesc(rou_range,rin_range,capacity')
hold on 
[C,h] = contour(rou_range,rin_range,capacity',[0.3 0.3], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
[C,h] = contour(rou_range,rin_range,capacity',[0.2 0.2], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
colorbar
axis xy
xlabel('rou')
ylabel('r_{in}')
title('Theoretical Capacity')
axis square
caxis([0.1 0.4])

load numerical_capacity_map.mat
figure,imagesc(a_range,rin_range,n_capacity)
hold on 
[C,h] = contour(a_range,rin_range,n_capacity,[0.3 0.3]*0.9, 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
[C,h] = contour(a_range,rin_range,n_capacity,[0.2 0.2]*0.9, 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
colorbar 
axis xy
xlabel('rou')
ylabel('r_{in}')
title('Numerical Capacity')
axis square
caxis([0.1 0.4])