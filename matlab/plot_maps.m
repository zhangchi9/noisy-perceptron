clear
clc
%close all
beta = 0;
rou_range = [0.1,5:5:100];

f = 0.2;

load fmincon_structure_map_load_at_numerical_capacity_th120and_140.mat
rin_range = rin_ind_range;
figure,imagesc(a_range,rin_range,Pcon_inh)
hold on 
[C,h] = contour(a_range,rin_range,Pcon_inh,[0.25,0.56], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
colorbar 
axis xy
xlabel('rou')
ylabel('r_{in}')
title('Pconinh')
axis square
caxis([0 0.7])

figure,imagesc(a_range,rin_range,Pcon_exc)
hold on 
[C,h] = contour(a_range,rin_range,Pcon_exc,[0.1,0.19], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
colorbar 
axis xy
xlabel('rou')
ylabel('r_{in}')
title('Pconexc')
axis square
caxis([0 0.3])

figure,imagesc(a_range,rin_range,CV_inh)
hold on 
[C,h] = contour(a_range,rin_range,CV_inh,[0.78,0.78], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
colorbar 
axis xy
xlabel('rou')
ylabel('r_{in}')
title('CVinh')
axis square
caxis([0.7 0.95])

figure,imagesc(a_range,rin_range,CV_exc)
hold on 
[C,h] = contour(a_range,rin_range,CV_exc,[0.85,1.1], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
colorbar 
axis xy
xlabel('rou')
ylabel('r_{in}')
title(['CVexc'])
axis square
caxis([0.75 0.95])
