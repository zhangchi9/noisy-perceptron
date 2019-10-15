clear
clc
close all
load dynamics_data.mat
load fig2data.mat


[C1, ~] = contour(rou_range,rin_range,Pconinh',[0.25 0.25]);
C1 = [C1(:,2:end),[100;0.2]];
[C2, ~] = contour(rou_range,rin_range,Pconinh',[0.56 0.56]);
C2 = [C1,fliplr(C2(:,2:end))];
close all
figure, fill(C2(1,:),C2(2,:),[1 0.9 0.1],'linestyle','none')
xlim([0,100])
ylim([0.01,0.2])
axis square
alpha(0.5)
set(gca, 'Color', 'none');
%set(gca,'xtick',[],'ytick',[])
export_fig 1.png -transparent 

[C1, ~] = contour(rou_range,rin_range,Pconexc',[0.1 0.1]);
C1 = [C1(:,2:end),[100;0.2]];
[C2, ~] = contour(rou_range,rin_range,Pconexc',[0.19 0.19]);
C2 = [C1,fliplr(C2(:,2:end))];
figure, fill(C2(1,:),C2(2,:),[1 0.7 0.3],'linestyle','none')
xlim([0,100])
ylim([0.01,0.2])
axis square
alpha(0.5)
set(gca, 'Color', 'none');
%set(gca,'xtick',[],'ytick',[])
export_fig 2.png -transparent 

[C1, ~] = contour(rou_range,rin_range,CVinh',[0.78,0.78]);
C1 = [C1(:,2:end),[100,100;0.01,0.2]];
figure, fill(C1(1,:),C1(2,:),[0.5 0.7 0.3],'linestyle','none')
xlim([0,100])
ylim([0.01,0.2])
axis square
alpha(0.5)
set(gca, 'Color', 'none');
%set(gca,'xtick',[],'ytick',[])
export_fig 3.png -transparent 

[C1, ~] = contour(rou_range,rin_range,CVexc',[0.85,0.85]);
C1 = [C1(:,2:end),[100,100;0.01,0.2]];
figure, fill(C1(1,:),C1(2,:),[0.5 1 0.3],'linestyle','none')
xlim([0,100])
ylim([0.01,0.2])
axis square
alpha(0.5)
set(gca, 'Color', 'none');
%set(gca,'xtick',[],'ytick',[])
export_fig 4.png -transparent 

% CV_ISI(CV_ISI==0) = nan;
% tmp = fillmissing(CV_ISI,'movmean',3);
% for i = 1:1000
%     tmp = fillmissing(tmp,'movmean',3);
% end
% CV_ISI2 = smooth2a(tmp,3);
[C1, ~] = contour(a_range,rin_range,CV_ISI,[0.8,0.8]);
%C1 = [C1(:,2:31),[5,5,;0.2,0.01]];
C1 = [C1(:,2:end)];
figure, fill(C1(1,:),C1(2,:),[0 1 0.3],'linestyle','none')
xlim([0,100])
ylim([0.01,0.2])
axis square
alpha(0.5)
set(gca, 'Color', 'none');
%set(gca,'xtick',[],'ytick',[])
export_fig 5.png -transparent 

% COR_I(COR_I==0) = nan;
% tmp = fillmissing(COR_I,'movmean',3);
% for i = 1:1000
%     tmp = fillmissing(tmp,'movmean',3);
% end
% COR_I2 = smooth2a(tmp,1);
[C1, ~] = contour(a_range,rin_range,COR_I,[0.4,0.4]);
C1 = [C1(:,2:end),[0,0,C1(1,2);0.2,0.005,0.005]];
figure, fill(C1(1,:),C1(2,:),[0.5 0.2 0.3],'linestyle','none')
xlim([0,100])
ylim([0.01,0.2])
axis square
alpha(0.5)
set(gca, 'Color', 'none');
%set(gca,'xtick',[],'ytick',[])
export_fig 6.png -transparent 

% SPKS_COR(SPKS_COR==0) = nan;
% tmp = fillmissing(SPKS_COR,'movmean',3);
% for i = 1:1000
%     tmp = fillmissing(tmp,'movmean',3);
% end
% SPKS_COR2 = smooth2a(tmp,2);
% figure,[C1, ~] = contour(a_range,rin_range,SPKS_COR,[0.04,0.04]);
% C1 = [C1(:,2:end)];
% [C2, ~] = contour(a_range,rin_range,SPKS_COR,[0.15,0.15]);
% C2 = [C1,fliplr(C2(:,2:end))];
% figure, fill(C2(1,:),C2(2,:),[0.5 0.2 0.7],'linestyle','none')
% xlim([5,100])
% ylim([0.01,0.2])
% axis square
% alpha(0.5)
% set(gca, 'Color', 'none');
% set(gca,'xtick',[],'ytick',[])
% export_fig 7.png -transparent 

close all