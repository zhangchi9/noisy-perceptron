clear
clc
close all
max_error_fac()
% load prob_retrieval_map_load85_different_fac.mat
% figure,imagesc(a_range,rin_range,fac0),axis xy,axis square, colorbar, title('no noise'),caxis([0,1])
% figure,imagesc(a_range,rin_range,fac6),axis xy,axis square, colorbar, title('60%'),caxis([0,1])
% figure,imagesc(a_range,rin_range,fac8),axis xy,axis square, colorbar, title('80%'),caxis([0,1])
% figure,imagesc(a_range,rin_range,fac10),axis xy,axis square, colorbar, title('100%'),caxis([0,1])
% figure,imagesc(a_range,rin_range,fac12),axis xy,axis square, colorbar, title('120%'),caxis([0,1])
% figure,imagesc(a_range,rin_range,fac14),axis xy,axis square, colorbar, title('140%'),caxis([0,1])
% figure,imagesc(a_range,rin_range,fac16),axis xy,axis square, colorbar, title('160%'),caxis([0,1])
% figure,imagesc(a_range,rin_range,fac18),axis xy,axis square, colorbar, title('180%'),caxis([0,1])
% figure,imagesc(a_range,rin_range,fac20),axis xy,axis square, colorbar, title('200%'),caxis([0,1])

load max_error.mat
figure,imagesc(a_range,rin_range,max_err_map),axis xy,axis square, colorbar, title('max amount noise')
figure,imagesc(a_range,rin_range,max_err_fac),axis xy,axis square, colorbar, title('max amount noise fac'),caxis([0.9,1.2])
hold on 
[C,h] = contour(a_range,rin_range,max_err_fac,[1 1.5], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)


load('retrieval_length_map_load85.mat')
load fig2data.mat
figure,imagesc(a_range,rin_range,retrieval_length),axis xy,axis square, colorbar, title('retrieval length')
figure,imagesc(a_range,rin_range,retrieval_length./round(capacity'*400)),axis xy,axis square, colorbar, title('retrieval fraction')
info = retrieval_length.*I'./capacity';
[~,b] = max(info);
figure,imagesc(a_range,rin_range,info),axis xy,axis square, colorbar, title('Information per neuron per retrieved sequence')
hold on 
[C,h] = contour(a_range,rin_range,info,[5 10 15 20 30], 'LineColor', 'k','LineStyle','--','LineWidth',2); 
clabel(C,h)
plot(a_range,rin_range(b),'r')

% load prob_retrieval_map_load85_same_noise.mat
% figure,imagesc(a_range,rin_range,noise5),axis xy,axis square, colorbar, title('noise5'),caxis([0,1])
% figure,imagesc(a_range,rin_range,noise10),axis xy,axis square, colorbar, title('noise10'),caxis([0,1])
% figure,imagesc(a_range,rin_range,noise15),axis xy,axis square, colorbar, title('noise15'),caxis([0,1])
% figure,imagesc(a_range,rin_range,noise20),axis xy,axis square, colorbar, title('noise20'),caxis([0,1])
% figure,imagesc(a_range,rin_range,noise25),axis xy,axis square, colorbar, title('noise25'),caxis([0,1])
% figure,imagesc(a_range,rin_range,noise30),axis xy,axis square, colorbar, title('noise30'),caxis([0,1])
% figure,imagesc(a_range,rin_range,noise35),axis xy,axis square, colorbar, title('noise35'),caxis([0,1])
% figure,imagesc(a_range,rin_range,noise40),axis xy,axis square, colorbar, title('noise40'),caxis([0,1])
% figure,imagesc(a_range,rin_range,noise45),axis xy,axis square, colorbar, title('noise45'),caxis([0,1])
% figure,imagesc(a_range,rin_range,noise50),axis xy,axis square, colorbar, title('noise50'),caxis([0,1])
% figure,imagesc(a_range,rin_range,noise55),axis xy,axis square, colorbar, title('noise55'),caxis([0,1])
% figure,imagesc(a_range,rin_range,noise60),axis xy,axis square, colorbar, title('noise60'),caxis([0,1])
% figure,imagesc(a_range,rin_range,noise65),axis xy,axis square, colorbar, title('noise65'),caxis([0,1])
% figure,imagesc(a_range,rin_range,noise70),axis xy,axis square, colorbar, title('noise70'),caxis([0,1])
% figure,imagesc(a_range,rin_range,noise75),axis xy,axis square, colorbar, title('noise75'),caxis([0,1])
% figure,imagesc(a_range,rin_range,noise80),axis xy,axis square, colorbar, title('noise80'),caxis([0,1])
% figure,imagesc(a_range,rin_range,noise85),axis xy,axis square, colorbar, title('noise85'),caxis([0,1])
% figure,imagesc(a_range,rin_range,noise90),axis xy,axis square, colorbar, title('noise90'),caxis([0,1])
% figure,imagesc(a_range,rin_range,noise95),axis xy,axis square, colorbar, title('noise95'),caxis([0,1])
% figure,imagesc(a_range,rin_range,noise100),axis xy,axis square, colorbar, title('noise100'),caxis([0,1])