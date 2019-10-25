clear
clc
data1 = load('feasible_not_solved_by_perceptron.mat');
sample_interval = data1.sample_interval;
rout = data1.rout;
W1 = data1.W;
w1 = data1.w;
n_epochs1 = data1.n_epochs;
rout_perceptron1 = data1.rout_perceptron;

data2 = load('feasible_solved_by_perceptron.mat');
W2 = data2.W;
w2 = data2.w;
n_epochs2 = data2.n_epochs;
rout_perceptron2 = data2.rout_perceptron;

data3 = load('nonfeasible_not_solved_by_perceptron.mat');
W3 = data3.W;
w3 = data3.w;
n_epochs3 = data3.n_epochs;
rout_perceptron3 = data3.rout_perceptron;

figure(3),hold on
plot(sample_interval:sample_interval:n_epochs2,rout_perceptron1(1:length(rout_perceptron2)),'b')
plot(sample_interval:sample_interval:n_epochs2,rout_perceptron2,'g')
plot(sample_interval:sample_interval:n_epochs2,rout_perceptron3(1:length(rout_perceptron2)),'r')
%plot(sample_interval:sample_interval:n_epochs,rout_W*ones(1,length(rout_perceptron)),'r--')
plot(sample_interval:sample_interval:n_epochs2,rout*ones(1,length(rout_perceptron2)),'k')
axis square
legend('feasible not solved by perceptron','feasible solved by perceptron','nonfeasible not solved by perceptron')
xlim([0 5*10^5])
ylim([0 0.2])
box on 

% Use matlab fit function 
% fo = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[0,-Inf],...
%                'Upper',[Inf,Inf],...
%                'StartPoint',[1 1]);
% ft = fittype('a*x+b','options',fo);
% [fitobject,gof,output] = fit(W1,w1,ft)
corrcoef(W1,w1)
corrcoef(W2,w2)
corrcoef(W3,w3)

% Use ols to fit 
X1 = [ones(length(W1),1),W1];
X2 = [ones(length(W2),1),W2];
X3 = [ones(length(W3),1),W3];
coeff1 = (X1'*X1)^-1*X1'*w1;
coeff2 = (X2'*X2)^-1*X2'*w2;
coeff3 = (X3'*X3)^-1*X3'*w3;

x = -1100:1100;
figure(4),hold on
plot(W1,w1,'bo'),hold on
plot(W2,w2,'go')
plot(W3,w3,'ro')
plot(x,coeff1(1)+coeff1(2)*x,'b')
plot(x,coeff2(1)+coeff2(2)*x,'g')
plot(x,coeff3(1)+coeff3(2)*x,'r')
%plot(1.1*[min(W1),max(W1)],1.1*[min(W1),max(W1)])
xlim([-1100,1100])
ylim([-1100,1100])
axis square
legend('feasible not solved by perceptron','feasible solved by perceptron','nonfeasible not solved by perceptron')
box on 