%% Requires ploterr function: https://www.mathworks.com/matlabcentral/fileexchange/22216-ploterr

clear
clc
clf

folder = "/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Analysis_Spreadsheets/";
cd(folder)
D = readtable('Confocal_DT_Table.csv');

%% - Normalize
% Normalize HC count to the length of each turn; scale the count to 212.5um

D.Apex_IHC_Norm     = (D.Apex_IHC ./ D.Apex_length) .* 212.5;
D.Apex_OHC_Norm     = (D.Apex_OHC ./ D.Apex_length) .* 212.5;
D.Middle_IHC_Norm   = (D.Middle_IHC ./ D.Middle_length) .* 212.5;
D.Middle_OHC_Norm   = (D.Middle_OHC ./ D.Middle_length) .* 212.5;
D.Base_IHC_Norm     = (D.Base_IHC ./ D.Base_length) .* 212.5;
D.Base_OHC_Norm     = (D.Base_OHC ./ D.Base_length) .* 212.5;

%% - Divide into 4 tables
% Single_DT -> singly-injected mice, DT genotype
% Single_WT -> singly-injected mice, WT genotype
% Double_DT -> doubly-injected mice, DT genotype
% Double_WT -> doubly-injected mice, WT genotype

% In the "Genotype_DT1_WT0" column, a value of 1 means DT, and 0 means WT
% (for genotype)

single_DT_idx   = (D.Genotype_DT1_WT0 == 1 & D.N_injections == 1);
Single_DT       = D(single_DT_idx,:); 

single_WT_idx   = (D.Genotype_DT1_WT0 == 0 & D.N_injections == 1);
Single_WT       = D(single_WT_idx,:); 

double_DT_idx   = (D.Genotype_DT1_WT0 == 1 & D.N_injections == 2);
Double_DT       = D(double_DT_idx,:); 

double_WT_idx   = (D.Genotype_DT1_WT0 == 0 & D.N_injections == 2);
Double_WT       = D(double_WT_idx,:); 

%% - Compile the normalized cell counts, use NaN's to equalize row lengths
% Row   Contents
% 1     Single, WT, IHC, Apex
% 2     Single, WT, IHC, Middle
% 3     Single, WT, IHC, Base
% 4     Single, DT, IHC, Apex
% 5     Single, DT, IHC, Middle
% 6     Single, DT, IHC, Base
% 7     Single, WT, OHC, Apex
% 8     Single, WT, OHC, Middle
% 9     Single, WT, OHC, Base
% 10    Single, DT, OHC, Apex
% 11    Single, DT, OHC, Middle
% 12    Single, DT, OHC, Base
% 13    Double, WT, IHC, Apex
% 14    Double, WT, IHC, Middle
% 15    Double, WT, IHC, Base
% 16    Double, DT, IHC, Apex
% 17    Double, DT, IHC, Middle
% 18    Double, DT, IHC, Base
% 19    Double, WT, OHC, Apex
% 20    Double, WT, OHC, Middle
% 21    Double, WT, OHC, Base
% 22    Double, DT, OHC, Apex
% 23    Double, DT, OHC, Middle
% 24    Double, DT, OHC, Base

clear compiled
compiled(:,1) = Single_WT.Apex_IHC_Norm;
compiled(:,2) = Single_WT.Middle_IHC_Norm;
compiled(:,3) = Single_WT.Base_IHC_Norm;

compiled(:,4) = Single_DT.Apex_IHC_Norm;
compiled(:,5) = Single_DT.Middle_IHC_Norm;
compiled(:,6) = Single_DT.Base_IHC_Norm;

compiled(:,7) = Single_WT.Apex_OHC_Norm;
compiled(:,8) = Single_WT.Middle_OHC_Norm;
compiled(:,9) = Single_WT.Base_OHC_Norm;

compiled(:,10) = Single_DT.Apex_OHC_Norm;
compiled(:,11) = Single_DT.Middle_OHC_Norm;
compiled(:,12) = Single_DT.Base_OHC_Norm;

compiled(end+1:7,:) = missing;

compiled(:,13) = Double_WT.Apex_IHC_Norm;
compiled(:,14) = Double_WT.Middle_IHC_Norm;
compiled(:,15) = Double_WT.Base_IHC_Norm;

compiled(end+1:12,:) = missing;

compiled(:,16) = Double_DT.Apex_IHC_Norm;
compiled(:,17) = Double_DT.Middle_IHC_Norm;
compiled(:,18) = Double_DT.Base_IHC_Norm;

% Make a temporary array in which I can add NaN's before adding it to the
% main array
clear temp_dat
temp_dat(:,1) = Double_WT.Apex_OHC_Norm;
temp_dat(:,2) = Double_WT.Middle_OHC_Norm;
temp_dat(:,3) = Double_WT.Base_OHC_Norm;

temp_dat(end+1:12,:) = missing;

compiled(:,19) = temp_dat(:,1);
compiled(:,20) = temp_dat(:,2);
compiled(:,21) = temp_dat(:,3);

compiled(:,22) = Double_DT.Apex_OHC_Norm;
compiled(:,23) = Double_DT.Middle_OHC_Norm;
compiled(:,24) = Double_DT.Base_OHC_Norm;

%% - Plot bar, swarm, error bars
clf
subplot(2,2,1)
bar(1,mean(Single_WT.Apex_IHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)
hold on
bar(2,mean(Single_WT.Middle_IHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)
bar(3,mean(Single_WT.Base_IHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)
bar(4,mean(Single_DT.Apex_IHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)
bar(5,mean(Single_DT.Middle_IHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)
bar(6,mean(Single_DT.Base_IHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)

s = swarmchart([1 2 3 4 5 6], compiled(:,1:6),20,'k');
s(1).XJitter = 'rand';
s(1).XJitterWidth = 1.0;

ploterr(1, mean(Single_WT.Apex_IHC_Norm,'omitnan'), [], std(Single_WT.Apex_IHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);
ploterr(2, mean(Single_WT.Middle_IHC_Norm,'omitnan'), [], std(Single_WT.Middle_IHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);
ploterr(3, mean(Single_WT.Base_IHC_Norm,'omitnan'), [], std(Single_WT.Base_IHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);
ploterr(4, mean(Single_DT.Apex_IHC_Norm,'omitnan'), [], std(Single_DT.Apex_IHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);
ploterr(5, mean(Single_DT.Middle_IHC_Norm,'omitnan'), [], std(Single_DT.Middle_IHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);
ploterr(6, mean(Single_DT.Base_IHC_Norm,'omitnan'), [], std(Single_DT.Base_IHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);

text(0.9,-2,'A','fontsize',16); text(1.9,-2,'M','fontsize',16); text(2.9,-2,'B','fontsize',16)
text(3.9,-2,'A','fontsize',16); text(4.9,-2,'M','fontsize',16); text(5.9,-2,'B','fontsize',16)
ylabel('IHCs','fontsize',18)
title('Single-Injected','fontsize',16)
ylim([0 40])
set(gca,'xticklabel',{[]},'fontsize',16)


subplot(2,2,3)
bar(1,mean(Single_WT.Apex_OHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)
hold on
bar(2,mean(Single_WT.Middle_OHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)
bar(3,mean(Single_WT.Base_OHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)
bar(4,mean(Single_DT.Apex_OHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)
bar(5,mean(Single_DT.Middle_OHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)
bar(6,mean(Single_DT.Base_OHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)

s = swarmchart([1 2 3 4 5 6], compiled(:,7:12),20,'k');
s(1).XJitter = 'rand';
s(1).XJitterWidth = 1;

ploterr(1, mean(Single_WT.Apex_OHC_Norm,'omitnan'), [], std(Single_WT.Apex_OHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);
ploterr(2, mean(Single_WT.Middle_OHC_Norm,'omitnan'), [], std(Single_WT.Middle_OHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);
ploterr(3, mean(Single_WT.Base_OHC_Norm,'omitnan'), [], std(Single_WT.Base_OHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);
ploterr(4, mean(Single_DT.Apex_OHC_Norm,'omitnan'), [], std(Single_DT.Apex_OHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);
ploterr(5, mean(Single_DT.Middle_OHC_Norm,'omitnan'), [], std(Single_DT.Middle_OHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);
ploterr(6, mean(Single_DT.Base_OHC_Norm,'omitnan'), [], std(Single_DT.Base_OHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);

text(0.9,-5,'A','fontsize',16); text(1.9,-5,'M','fontsize',16); text(2.9,-5,'B','fontsize',16)
text(3.9,-5,'A','fontsize',16); text(4.9,-5,'M','fontsize',16); text(5.9,-5,'B','fontsize',16)
text(1.8,-14,'WT','fontsize',18)
text(4.8,-14,'DT','fontsize',18)
ylabel('OHCs','fontsize',18)
set(gca,'xticklabel',{[]},'fontsize',16)
ylim([0 100])


subplot(2,2,2)
bar(1,mean(Double_WT.Apex_IHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)
hold on
bar(2,mean(Double_WT.Middle_IHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)
bar(3,mean(Double_WT.Base_IHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)
bar(4,mean(Double_DT.Apex_IHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)
bar(5,mean(Double_DT.Middle_IHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)
bar(6,mean(Double_DT.Base_IHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)

s = swarmchart([1 2 3 4 5 6], compiled(:,13:18),20,'k');
s(1).XJitter = 'rand';
s(1).XJitterWidth = 1;

ploterr(1, mean(Double_WT.Apex_IHC_Norm,'omitnan'), [], std(Double_WT.Apex_IHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);
ploterr(2, mean(Double_WT.Middle_IHC_Norm,'omitnan'), [], std(Double_WT.Middle_IHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);
ploterr(3, mean(Double_WT.Base_IHC_Norm,'omitnan'), [], std(Double_WT.Base_IHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);
ploterr(4, mean(Double_DT.Apex_IHC_Norm,'omitnan'), [], std(Double_DT.Apex_IHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);
ploterr(5, mean(Double_DT.Middle_IHC_Norm,'omitnan'), [], std(Double_DT.Middle_IHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);
ploterr(6, mean(Double_DT.Base_IHC_Norm,'omitnan'), [], std(Double_DT.Base_IHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);

text(0.9,-2,'A','fontsize',16); text(1.9,-2,'M','fontsize',16); text(2.9,-2,'B','fontsize',16)
text(3.9,-2,'A','fontsize',16); text(4.9,-2,'M','fontsize',16); text(5.9,-2,'B','fontsize',16)
title('Double-Injected','fontsize',16)
ylim([0 40])
set(gca,'xticklabel',{[]})


subplot(2,2,4)
bar(1,mean(Double_WT.Apex_OHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)
hold on
bar(2,mean(Double_WT.Middle_OHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)
bar(3,mean(Double_WT.Base_OHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)
bar(4,mean(Double_DT.Apex_OHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)
bar(5,mean(Double_DT.Middle_OHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)
bar(6,mean(Double_DT.Base_OHC_Norm,'omitnan'),'facecolor',[1 1 1],'linestyle','-','facealpha',0.8)

s = swarmchart([1 2 3 4 5 6], compiled(:,19:24),20,'k');
s(1).XJitter = 'rand';
s(1).XJitterWidth = 1;

ploterr(1, mean(Double_WT.Apex_OHC_Norm,'omitnan'), [], std(Double_WT.Apex_OHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);
ploterr(2, mean(Double_WT.Middle_OHC_Norm,'omitnan'), [], std(Double_WT.Middle_OHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);
ploterr(3, mean(Double_WT.Base_OHC_Norm,'omitnan'), [], std(Double_WT.Base_OHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);
ploterr(4, mean(Double_DT.Apex_OHC_Norm,'omitnan'), [], std(Double_DT.Apex_OHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);
ploterr(5, mean(Double_DT.Middle_OHC_Norm,'omitnan'), [], std(Double_DT.Middle_OHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);
ploterr(6, mean(Double_DT.Base_OHC_Norm,'omitnan'), [], std(Double_DT.Base_OHC_Norm,'omitnan'), 'k', 'abshhxy', 0.3);

text(0.9,-5,'A','fontsize',16); text(1.9,-5,'M','fontsize',16); text(2.9,-5,'B','fontsize',16)
text(3.9,-5,'A','fontsize',16); text(4.9,-5,'M','fontsize',16); text(5.9,-5,'B','fontsize',16)
text(1.8,-14,'WT','fontsize',18)
text(4.8,-14,'DT','fontsize',18)
set(gca,'xticklabel',{[]})
ylim([0 100])



%% - Run Kruskal-Wallis and add *, ** or NS to plots
% *  = p < 0.05
% ** = p < 0.001
% NS = p > 0.05≤

% Test for normality
kstest([compiled(:,1);compiled(:,2);compiled(:,3)])
kstest([compiled(:,4);compiled(:,5);compiled(:,6)])
% Not normal

subplot(2,2,1)
%[p, tbl, stats] = anova1(compiled(:,1:6),{'WT', 'WT', 'WT','DT', 'DT', 'DT'},'off');
[p] = KruskalWallisPlot([compiled(:,1);compiled(:,2);compiled(:,3)],[compiled(:,4);compiled(:,5);compiled(:,6)], 1, 35, [2 4])
% hold on
% line([2 5],[35 35],'color','k','linewidth',2)
% %line([4 5],[35 35],'color','k','linewidth',2)
% if p >= 0.05
%     text(3.3,34,'NS','fontsize',30)
% end
% if p < 0.05 && p>= 0.01
%     text(3.3,35,'*','fontsize',30)
% end
% if p < 0.01 && p>=0.001
%     text(3.3,35,'**','fontsize',30)
% end
% if p < 0.001 
%     text(3.3,35,'***','fontsize',30)
% end


subplot(2,2,2)
%[p, tbl, stats] = anova1(compiled(:,13:18),{'WT', 'WT', 'WT','DT', 'DT', 'DT'},'off');
[p] = KruskalWallisPlot([compiled(:,13);compiled(:,14);compiled(:,15)],...
    [compiled(:,16);compiled(:,17);compiled(:,18)], 1, 35, [2 4])
% hold on
% line([2 3],[35 35],'color','k','linewidth',2)
% line([4 5],[35 35],'color','k','linewidth',2)
% if p >= 0.05
%     text(3.3,34,'NS','fontsize',30)
% end
% if p < 0.05 && p>= 0.01
%     text(3.3,35,'*','fontsize',30)
% end
% if p < 0.01 && p>=0.001
%     text(3.3,35,'**','fontsize',30)
% end
% if p < 0.001 
%     text(3.3,35,'***','fontsize',30)
% end


subplot(2,2,3)
%[p, tbl, stats] = anova1(compiled(:,7:12),{'WT', 'WT', 'WT','DT', 'DT', 'DT'},'off');
[p] = KruskalWallisPlot([compiled(:,7);compiled(:,8);compiled(:,9)],...
    [compiled(:,10);compiled(:,11);compiled(:,12)], 1, 100, [2 4])
% hold on
% line([2 3],[95 95],'color','k','linewidth',2)
% line([4 5],[95 95],'color','k','linewidth',2)
% if p >= 0.05
%     text(3.3,90,'NS','fontsize',30)
% end
% if p < 0.05 && p>= 0.01
%     text(3.3,90,'*','fontsize',30)
% end
% if p < 0.01 && p>=0.001
%     text(3.3,90,'**','fontsize',30)
% end
% if p < 0.001 
%     text(3.3,90,'***','fontsize',30)
% end


subplot(2,2,4)
%[p, tbl, stats] = anova1(compiled(:,19:24),{'WT', 'WT', 'WT','DT', 'DT', 'DT'},'off');
[p] = KruskalWallisPlot([compiled(:,19);compiled(:,20);compiled(:,21)],...
    [compiled(:,22);compiled(:,23);compiled(:,24)], 1, 100, [2 4])

% hold on
% line([2 3],[95 95],'color','k','linewidth',2)
% line([4 5],[95 95],'color','k','linewidth',2)
% if p >= 0.05
%     text(3.3,90,'NS','fontsize',30)
% end
% if p < 0.05 && p>= 0.01
%     text(3.3,90,'*','fontsize',30)
% end
% if p < 0.01 && p>=0.001
%     text(3.3,90,'**','fontsize',30)
% end
% if p < 0.001 
%     text(3.3,90,'***','fontsize',30)
% end

cd("/Users/julien/Library/CloudStorage/Box-Box/Julien Azimzadeh's Files/Research/Ricci Lab/Data/Results")
exportgraphics(gcf, 'Confocal_HairCellCounts.pdf', 'ContentType', 'vector');

%% - T test
%std([compiled(:,1); compiled(:,2); compiled(:,3)],'omitnan');
%std([compiled(:,4); compiled(:,5); compiled(:,6)],'omitnan');

% Make a variable that contains the mean count from each cochlea
% Single, IHCs
Count_means(1:6,1) = mean(compiled(1:6,1:3),2,'omitnan');
Count_means(1:6,2) = 1; % WT
Count_means(7:12,1) = mean(compiled(1:6,4:6),2,'omitnan');
Count_means(7:12,2) = 2; % DT

p1 = vartestn(Count_means(1:12,1),Count_means(1:12,2),'TestType','LeveneAbsolute');

%Double, IHCs
Count_means(13:18,1) = mean(compiled(1:6,7:9),2,'omitnan');
Count_means(13:18,2) = 3; % WT
Count_means(19:24,1) = mean(compiled(1:6,10:12),2,'omitnan');
Count_means(19:24,2) = 4; % DT

p2 = vartestn(Count_means(13:24,1),Count_means(13:24,2),'TestType','LeveneAbsolute');

%%
Mean_Single_WT_IHC = mean([mean(Single_WT.Apex_IHC_Norm), mean(Single_WT.Middle_IHC_Norm), nanmean(Single_WT.Base_IHC_Norm)]);
Mean_Single_WT_OHC = mean([mean(Single_WT.Apex_OHC_Norm), mean(Single_WT.Middle_OHC_Norm), nanmean(Single_WT.Base_OHC_Norm)]);
Mean_Double_WT_IHC = mean([mean(Double_WT.Apex_IHC_Norm), mean(Double_WT.Middle_IHC_Norm), nanmean(Double_WT.Base_IHC_Norm)]);
Mean_Double_WT_OHC = mean([mean(Double_WT.Apex_OHC_Norm), mean(Double_WT.Middle_OHC_Norm), nanmean(Double_WT.Base_OHC_Norm)]);

Mean_Single_DT_IHC = mean([mean(Single_DT.Apex_IHC_Norm), mean(Single_DT.Middle_IHC_Norm), nanmean(Single_DT.Base_IHC_Norm)]);
Mean_Single_DT_OHC = mean([mean(Single_DT.Apex_OHC_Norm), mean(Single_DT.Middle_OHC_Norm), nanmean(Single_DT.Base_OHC_Norm)]);
Mean_Double_DT_IHC = mean([mean(Double_DT.Apex_IHC_Norm), mean(Double_DT.Middle_IHC_Norm), nanmean(Double_DT.Base_IHC_Norm)]);
Mean_Double_DT_OHC = mean([mean(Double_DT.Apex_OHC_Norm), mean(Double_DT.Middle_OHC_Norm), nanmean(Double_DT.Base_OHC_Norm)]);

