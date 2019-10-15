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

for i = 1:length(noise_added)

retrieval_prob_tmp = retrieval_prob(:,:,i);
retrieval_length_tmp = retrieval_length(:,:,i); 

W = retrieval_prob_tmp.*I;
[a1,b1] = max(mean(W,1));
compelete_infor_max(i) = noise_added(b1);
V = retrieval_length_tmp./m_capacity.*I;
[a2,b2] = max(mean(V,1));
retrievable_infor_max(i) = noise_added(b2);
end

figure,plot(noise_added,compelete_infor_max)
hold on, plot(noise_added,retrievable_infor_max)
plot(noise_added,noise_added,'k')
axis square
xlabel('noise in retrieval')
ylabel('noise in learning')
