% this function will combine the strutural data and dynamical data

clear
clc

perceptron200 = nan(12,1);
perceptron400 = nan(12,1);
perceptron800 = nan(12,1);
fmincon200 = nan(12,1);
fmincon400 = nan(12,1);
fmincon800 = nan(12,1);
replica = nan(12,1);

perceptron200(1) = 0.1243;
perceptron400(1) = 0.1354;
perceptron800(1) = 0.1466;
fmincon200(1) = 0.1492;
fmincon400(1) = 0.1639;
fmincon800(1) = 0.1764;
replica(1) = 0.1872;

s = load('structural data.mat');
replica(6:9) = [s.Pconinh_theory; s.Pconexc_theory; s.CVinh_theory; s.CVexc_theory];


s = load('/home/chi/Dropbox/Research/Project_Perceptron/results/perceptron_200_structures.mat');
d = load('/home/chi/Dropbox/Research/Project_Perceptron/results/perceptron_200_dynamics.mat');
r = load('/home/chi/Dropbox/Research/Project_Perceptron/results/perceptron_200_retrieval.mat');
comp_retrieval_info = perceptron200(1)*0.6109*r.prob;
retrieval_info = perceptron200(1)*0.6109*r.length_retrieval_fraction;
perceptron200(2:12) = [r.prob;r.length_retrieval_fraction;comp_retrieval_info;retrieval_info;s.Pcon_inh;s.Pcon_exc;s.CV_inh;s.CV_exc;d.CV_ISI_mean2;d.SPKS_COR_mean;d.COR_I_mean];

s = load('/home/chi/Dropbox/Research/Project_Perceptron/results/perceptron_400_structures.mat');
d = load('/home/chi/Dropbox/Research/Project_Perceptron/results/perceptron_400_dynamics.mat');
r = load('/home/chi/Dropbox/Research/Project_Perceptron/results/perceptron_400_retrieval.mat');
comp_retrieval_info = perceptron400(1)*0.6109*r.prob;
retrieval_info = perceptron400(1)*0.6109*r.length_retrieval_fraction;
perceptron400(2:12) = [r.prob;r.length_retrieval_fraction;comp_retrieval_info;retrieval_info;s.Pcon_inh;s.Pcon_exc;s.CV_inh;s.CV_exc;d.CV_ISI_mean2;d.SPKS_COR_mean;d.COR_I_mean];

s = load('/home/chi/Dropbox/Research/Project_Perceptron/results/perceptron_800_structures.mat');
d = load('/home/chi/Dropbox/Research/Project_Perceptron/results/perceptron_800_dynamics.mat');
r = load('/home/chi/Dropbox/Research/Project_Perceptron/results/perceptron_800_retrieval.mat');
comp_retrieval_info = perceptron800(1)*0.6109*r.prob;
retrieval_info = perceptron800(1)*0.6109*r.length_retrieval_fraction;
perceptron800(2:12) = [r.prob;r.length_retrieval_fraction;comp_retrieval_info;retrieval_info;s.Pcon_inh;s.Pcon_exc;s.CV_inh;s.CV_exc;d.CV_ISI_mean2;d.SPKS_COR_mean;d.COR_I_mean];

s = load('/home/chi/Dropbox/Research/Project_Perceptron/results/fmincon_200_structures.mat');
d = load('/home/chi/Dropbox/Research/Project_Perceptron/results/fmincon_200_dynamics.mat');
r = load('/home/chi/Dropbox/Research/Project_Perceptron/results/fmincon_200_retrieval.mat');
comp_retrieval_info = fmincon200(1)*0.6109*r.prob;
retrieval_info = fmincon200(1)*0.6109*r.length_retrieval_fraction;
fmincon200(2:12) = [r.prob;r.length_retrieval_fraction;comp_retrieval_info;retrieval_info;s.Pcon_inh;s.Pcon_exc;s.CV_inh;s.CV_exc;d.CV_ISI_mean2;d.SPKS_COR_mean;d.COR_I_mean];

s = load('/home/chi/Dropbox/Research/Project_Perceptron/results/fmincon_400_structures.mat');
d = load('/home/chi/Dropbox/Research/Project_Perceptron/results/fmincon_400_dynamics.mat');
r = load('/home/chi/Dropbox/Research/Project_Perceptron/results/fmincon_400_retrieval.mat');
comp_retrieval_info = fmincon400(1)*0.6109*r.prob;
retrieval_info = fmincon400(1)*0.6109*r.length_retrieval_fraction;
fmincon400(2:12) = [r.prob;r.length_retrieval_fraction;comp_retrieval_info;retrieval_info;s.Pcon_inh;s.Pcon_exc;s.CV_inh;s.CV_exc;d.CV_ISI_mean2;d.SPKS_COR_mean;d.COR_I_mean];

s = load('/home/chi/Dropbox/Research/Project_Perceptron/results/fmincon_800_structures.mat');
d = load('/home/chi/Dropbox/Research/Project_Perceptron/results/fmincon_800_dynamics.mat');
r = load('/home/chi/Dropbox/Research/Project_Perceptron/results/fmincon_800_retrieval.mat');
comp_retrieval_info = fmincon800(1)*0.6109*r.prob;
retrieval_info = fmincon800(1)*0.6109*r.length_retrieval_fraction;
fmincon800(2:12) = [r.prob;r.length_retrieval_fraction;comp_retrieval_info;retrieval_info;s.Pcon_inh;s.Pcon_exc;s.CV_inh;s.CV_exc;d.CV_ISI_mean2;d.SPKS_COR_mean;d.COR_I_mean];

data = [perceptron200, perceptron400, perceptron800,fmincon200,fmincon400,fmincon800, replica];

ax = axes;
pbaspect([3 1 1])
h = bar(data,'BarWidth',1);
h(1).FaceColor = [0 1 1];
h(2).FaceColor = [0,0.5,0.5];
h(3).FaceColor = [0,0.4,0.6];
h(4).FaceColor = [1 0 1];
h(5).FaceColor = [0.5,0,0.5];
h(6).FaceColor = [0.4,0,0.6];
h(7).FaceColor = [0.5 0.5 0.5];
ax.YGrid = 'on';
ax.GridLineStyle = '-';
xticks(ax,[1:12]);
% Naming each of the bar groups
xticklabels(ax,{'Capacity','Retrieval prob','retrieval length fraction','comp information','retrieved information','P_{con}^{inh}', 'P_{con}^{exc}', 'CV^{inh}', 'CV^{exc}', 'CV ISI', 'SPKS COR', 'COR I'});
lg = legend('perceptron N=200','perceptron N = 400','perceptron N = 800','fmincon N=200','fmincon N=400','fmincon N=800','Theory','AutoUpdate','off');
%lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;
% Finding the number of groups and the number of bars in each group
[ngroups,nbars] = size(data);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
xlim([0,13])
ylim([0 1.15])
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
 
% for i = 1:nbars
%     % Calculate center of each bar
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%     errorbar(x, data(:,i), errbar(:,i), 'k', 'linestyle', 'none');
% end

