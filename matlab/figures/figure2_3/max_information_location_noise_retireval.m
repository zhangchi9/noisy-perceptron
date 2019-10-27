clear
clc

load('/home/chi/Dropbox/Research/Project_Perceptron/results/retrieval_prob_length_map_load_at_numerical_capacity_noise_all.mat')
%load('retrieval_prob_length_map_load85_no_noise.mat')
load('/home/chi/Dropbox/Research/Project_Perceptron/results/m_capacity.mat')

capacity = m_capacity/400;
fout = 0.2;

rin_range = 2.^rin_range_ind;

Is = fout*(rin_range/2/fout.*log(rin_range/2/fout) + (1-rin_range/2/fout).*log(1-rin_range/2/fout) -log(fout)) + ...
    (1-fout)*(rin_range/2/(1-fout).*log(rin_range/2/(1-fout)) + (1-rin_range/2/(1-fout)).*log(1-rin_range/2/(1-fout)) -log((1-fout)));
I = capacity.* (Is'*ones(1,length(a_range)));

rin_ind_range = (log2(rin_range));

% for i = 1:length(noise_added)
% 
% retrieval_prob_tmp = retrieval_prob(:,:,i);
% retrieval_length_tmp = retrieval_length(:,:,i); 
% 
% W = retrieval_prob_tmp.*I;
% Wmean = mean(W,1);
% grid = 0:0.1:100;
% Wmeanq = interp1(a_range,Wmean,grid,'spline');
% 
% % figure,plot(a_range,Wmean)
% % hold on 
% % plot(0:100,Wmeanq)
% 
% [a1,b1] = max(Wmeanq);
% compelete_infor_max(i) = grid(b1);
% 
% 
% V = retrieval_length_tmp./m_capacity.*I;
% Vmean = mean(V,1);
% Vmeanq = interp1(a_range,Vmean,grid,'spline');
% 
% [a2,b2] = max(Vmeanq);
% retrievable_infor_max(i) = grid(b2);
% end

for i = 1:length(noise_added)

retrieval_prob_tmp = retrieval_prob(:,:,i);
retrieval_length_tmp = retrieval_length(:,:,i); 

W = retrieval_prob_tmp.*I;
[Xq,Yq] = meshgrid(0:1:100,-2:-0.1:-12) ;
Wq = interp2(a_range,rin_range_ind,W,Xq,Yq);

% figure,imagesc(0:1:100,-2:-0.1:-12,Wq),axis xy
% figure,plot(a_range,mean(W,1),'o-')
% hold on 
% plot(0:1:100,mean(Wq,1),'*-')

[a1,b1] = max(mean(Wq,1));
compelete_infor_max(i) = b1;

V = retrieval_length_tmp./m_capacity.*I;
Vq = interp2(a_range,rin_range_ind,V,Xq,Yq);
[a2,b2] = max(mean(Vq,1));
retrievable_infor_max(i) = b2;
end

figure,plot(noise_added,compelete_infor_max)
hold on, plot(noise_added,retrievable_infor_max)
plot(noise_added,noise_added,'k')
axis square
xlabel('noise in retrieval')
ylabel('noise in learning')
xlim([0 70])
ylim([0 70])


figure
retrieval_prob_tmp = retrieval_prob(:,:,1);
plot(capacity(:),retrieval_prob_tmp(:),'o','color','r')
xlim([0 0.6])
axis square

figure
retrieval_prob_tmp = retrieval_prob(:,:,7);
plot(capacity(:),retrieval_prob_tmp(:),'o','color','b')
xlim([0 0.6])
axis square

figure
retrieval_prob_tmp = retrieval_prob(:,:,13);
plot(capacity(:),retrieval_prob_tmp(:),'o','color','k')
xlim([0 0.6])
axis square
%legend({'retrieval noise 0','retrieval noise 30','retrieval noise 60'})

