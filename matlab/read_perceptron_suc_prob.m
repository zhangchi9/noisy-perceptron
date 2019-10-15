clear
clc

file_obj = dir('/scratch/zhang.chi9/perceptron/data/perceptron_suc_prob/N_200*');
filename = {file_obj.name};
folder = {file_obj.folder};
sucprob200 = zeros(1,length([0.4:0.1:1.5]));
for i = 1:length(filename)
    load([folder{i},'/',filename{i}])
    sucprob200(i) = mean(routmin(1:200)<=rout);
end
prob_learning200 = [[0.4:0.1:1.5]*0.1872;sucprob200]';

file_obj = dir('/scratch/zhang.chi9/perceptron/data/perceptron_suc_prob/N_400*');
filename = {file_obj.name};
folder = {file_obj.folder};
sucprob400 = zeros(1,length(0.5:0.1:1.5));
for i = 1:length(filename)
    load([folder{i},'/',filename{i}])
    sucprob400(i) = mean(routmin(1:200)<=rout);
end
prob_learning400 = [[0.5:0.1:1.5]*0.1872;sucprob400]';

file_obj = dir('/scratch/zhang.chi9/perceptron/data/perceptron_suc_prob/N_800*');
filename = {file_obj.name};
folder = {file_obj.folder};
sucprob800 = zeros(1,length([0.5:0.1:1.0,1.2:0.1:1.4]));
for i = 1:length(filename)
    load([folder{i},'/',filename{i}])
    sucprob800(i) = mean(routmin(1:200)<=rout);
end
prob_learning800 = [[0.5:0.1:1.0,1.2:0.1:1.4]*0.1872;sucprob800]';

get_numerical_capacity(prob_learning200)
get_numerical_capacity(prob_learning400)
get_numerical_capacity(prob_learning800)


figure
plot([0.4:0.1:1.5]*0.1872,sucprob200,'*-') 
hold on 
plot((0.5:0.1:1.5)*0.1872,sucprob400,'*-') 
plot(([0.5:0.1:1.0,1.2:0.1:1.4])*0.1872,sucprob800,'*-') 