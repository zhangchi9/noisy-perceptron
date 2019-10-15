clear
clc

N = 200;
f = 0.2;
filepath = ['E:\Dropbox\Project_Perceptron\codes\N',num2str(N),'_homo_rin_rout_load\'];
rjRange = 10.^-(1:0.5:4);
routRange = 10.^-(1:0.5:4);
loadind = 1:10;
trialRange = 1:100;



index = 1;
for i = 1 : length(rjRange)
    for j = 1:length(routRange) 
        [i,j]
        WcrossTrail = [];
        rin = rjRange(i);
        rout = routRange(j);
        LoadRange = nan(1,10);
        sucProb = nan(1,10);
        
        
        JInhMean = nan(1,10);
        JInhMeanStd = nan(1,10);
        JExcMean = nan(1,10);
        JExcMeanStd = nan(1,10);
        PconInhMean = nan(1,10);
        PconInhStd = nan(1,10);
        PconExcMean= nan(1,10);
        PconExcStd= nan(1,10);
        CVInhMean= nan(1,10);
        CVInhStd= nan(1,10);
        CVExcMean= nan(1,10);
        CVExcStd = nan(1,10);
        KappaMean = nan(1,10);
        KappaStd = nan(1,10);
     
        for k = 1:length(loadind)
            Load = loadind(k);
            WCrossTrail = [];
            suc_single_trail=nan(1,100);
            mean_inh= nan(1,100);
            mean_exc= nan(1,100);
            Pcon_inh= nan(1,100);
            Pcon_exc= nan(1,100);
            CV_inh= nan(1,100);
            CV_exc= nan(1,100);
            std_exc= nan(1,100);
            std_inh= nan(1,100);
            kappa = nan(1,100);
            memory_load = nan;
            for ind = 1 : length(trialRange)
               
                trial = trialRange(ind);
                filename = [filepath,'N_',num2str(N),'_rj_',num2str(rin),'_rout_',num2str(rout),'_load_',num2str(Load),'_trial_',num2str(trial),'.mat'];
                if exist(filename,'file')
                    load(filename)
                    if exitflag > 0
                        p1 = rj/2/(1-f)*ones(N,1);
                        p2 = (1-f)*p1/f;
                        %WcrossTrail = [WcrossTrail,W];
                        suc_single_trail(ind)= epsilonsum<0.001;
                        [mean_inh(ind),mean_exc(ind),Pcon_inh(ind),Pcon_exc(ind),CV_inh(ind),CV_exc(ind),std_exc(ind),std_inh(ind)] = connection_prob_CV(W,20);
                        kappa(ind) = sum(((1-f)*p1.*(1-p1)+f*p2.*(1-p2)).*W.*W)/N;
                    end
                else
                    rerun(index,:) = [i,j,k,ind];
                    index = index + 1;
                end
            end
%             figure
%             hist(WcrossTrail(:),1000)
            JInhMean(k) = nanmean(mean_inh);
            JInhMeanStd(k)= nanstd(mean_inh)/sqrt(nnz(~isnan(mean_inh)));
            JExcMean(k) = nanmean(mean_exc);
            JExcMeanStd(k)= nanstd(mean_exc)/sqrt(nnz(~isnan(mean_exc)));
            PconInhMean(k)= nanmean(Pcon_inh);
            PconInhStd(k)= nanstd(Pcon_inh)/sqrt(nnz(~isnan(Pcon_inh)));
            PconExcMean(k)= nanmean(Pcon_exc);
            PconExcStd(k)= nanstd(Pcon_exc)/sqrt(nnz(~isnan(Pcon_exc)));
            CVInhMean(k)= nanmean(CV_inh);
            CVInhStd(k)= nanstd(CV_inh)/sqrt(nnz(~isnan(CV_inh)));
            CVExcMean(k) = nanmean(CV_exc);
            CVExcStd(k) = nanstd(CV_exc)/sqrt(nnz(~isnan(CV_exc)));
            KappaMean(k) = nanmean(kappa);
            KappaStd(k) = nanstd(kappa)/sqrt(nnz(~isnan(kappa)));
            LoadRange(k) = memory_load;
            sucProb(k) = nanmean(suc_single_trail);
        end
        SucRate{i,j} = [LoadRange;sucProb;PconInhMean;PconInhStd;PconExcMean;PconExcStd;CVInhMean;CVInhStd;CVExcMean;CVExcStd;KappaMean;KappaStd;JInhMean;JInhMeanStd;JExcMean;JExcMeanStd];
        %         figure,plot(LoadRange,sucProb)
        %         hold on
        %         capacity =  Robust_theoretical_solution_rin_rout(rin,rout);
        %         plot([capacity,capacity],[0,1])
    end
end

%unique(rerun(:,1:3),'rows');

%%% plot results %%%%%
for i = 2;
    for j = 1:7;
        rin = rjRange(i);
        rout = routRange(j);
        figure
        plot(SucRate{i,j}(1,:),SucRate{i,j}(2,:),'k')
        hold on
        plot(SucRate{i,j}(1,:),SucRate{i,j}(3,:),'--','color','r')
        plot(SucRate{i,j}(1,:),SucRate{i,j}(5,:),'color','b')
        plot(SucRate{i,j}(1,:),SucRate{i,j}(7,:),'--','color','r')
        plot(SucRate{i,j}(1,:),SucRate{i,j}(9,:),'color','b')
        legend('Suc Prob','Pcon inh','Pcon exc','CV inh','CV exc')
        [capacityTheory(i,j),exitflag(i,j),PconinhTheory(i,j),PconexcTheory(i,j),CVinhTheory(i,j),CVexcTheory(i,j),JmeanExcTheory(i,j),JmeanInhTheory(i,j),kappaTheory(i,j)] = Robust_homo_theoretical_solution_rin_rout(800,rin,rout);
        plot([capacityTheory(i,j),capacityTheory(i,j)],[0,1],'--')
        
        
        PconInh(i,j)= SucRate{i,j}(3,5);
        PconInh_errbar(i,j)= SucRate{i,j}(4,5);
        PconExc(i,j)= SucRate{i,j}(5,5);
        PconExc_errbar(i,j)= SucRate{i,j}(6,5);
        CVInh(i,j)= SucRate{i,j}(7,5);
        CVInh_errbar(i,j)= SucRate{i,j}(8,5);
        CVExc(i,j)= SucRate{i,j}(9,5);
        CVExc_errbar(i,j)= SucRate{i,j}(10,5);
        Kappa(i,j)= SucRate{i,j}(11,5);
        Kappa_errbar(i,j)= SucRate{i,j}(12,5);
        JInh(i,j)=SucRate{i,j}(13,5);
        JInh_errbar(i,j)=SucRate{i,j}(14,5);
        JExc(i,j)=SucRate{i,j}(15,5);
        JExc_errbar(i,j)=SucRate{i,j}(16,5);

    end
end

figure,
subplot(4,2,1),semilogx(log10(routRange),capacityTheory),xlabel('rout'),ylabel('capacity'),axis square
subplot(4,2,2),errorbar(log10(routRange),Kappa,Kappa_errbar),xlabel('rout'),ylabel('Kappa')
hold on, plot(log10(routRange),kappaTheory),axis square
subplot(4,2,3),errorbar(log10(routRange),PconInh,PconInh_errbar),xlabel('rout'),ylabel('Pcon inh')
hold on, plot(log10(routRange),PconinhTheory),axis square
subplot(4,2,4),errorbar(log10(routRange),PconExc,PconExc_errbar),xlabel('rout'),ylabel('Pcon exc')
hold on, plot(log10(routRange),PconexcTheory),axis square
subplot(4,2,5),errorbar(log10(routRange),CVInh,CVInh_errbar),xlabel('rout'),ylabel('CV inh')
hold on, plot(log10(routRange),CVinhTheory),axis square
subplot(4,2,6),errorbar(log10(routRange),CVExc,CVExc_errbar),xlabel('rout'),ylabel('CV exc')
hold on, plot(log10(routRange),CVexcTheory),axis square
subplot(4,2,7),errorbar(log10(routRange),JInh,JInh_errbar),xlabel('rout'),ylabel('Jinh Mean')
hold on, plot(log10(routRange),JmeanInhTheory),axis square
subplot(4,2,8),errorbar(log10(routRange),JExc,JExc_errbar),xlabel('rout'),ylabel('Jexc Mean')
hold on, plot(log10(routRange),JmeanExcTheory),axis square



% h = errorbar(routRange,PconInh,PconInh_errbar);
% set(gca,'YScale','log');


