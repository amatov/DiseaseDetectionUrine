function [ newMirs ] = miRnormalize

% miR = xlsread('healthy_miRlog2.xlsx');
%miR = xlsread('healthy_mircount.xlsx');

% convert miR to log2 and save to xls
%data1 = xlsread('healthy_mircount_full947noLabel.xlsx');

data1 = xlsread('lung_oncomirsQuaNor_743H2.xlsx');
aux=0;
for i= 1:13
    aux = aux + var(data1(i,:));
end
aux = aux/13
aux1=0;
for i =14:28
    aux1 = aux1 + var(data1(i,:));
end
aux1 = aux1/15
%-------panel 13 markers common between stage iv lung cancer
miR891a = data1(719,:)
PmiR891a = miR891a(1:13)
PmiR891a(14) = NaN
PmiR891a(15) = NaN
HmiR891a = miR891a(14:28)
p891a = ranksum(PmiR891a,HmiR891a)
figure, boxplot([PmiR891a',HmiR891a'],'Notch','off','Labels',{'Patient ','Control'},'Whisker',1,'PlotStyle','compact')
lgd891a = sprintf('p = %.5f', p891a);
ylabel('miR-891a-5p')
title(lgd891a)

miR196a = data1(194,:)
PmiR196a = miR196a(1:13)
PmiR196a(14) = NaN
PmiR196a(15) = NaN
HmiR196a = miR196a(14:28)
p196a = ranksum(PmiR196a,HmiR196a)
figure, boxplot([PmiR196a',HmiR196a'],'Notch','off','Labels',{'Patient ','Control'},'Whisker',1,'PlotStyle','compact')
lgd196a = sprintf('p = %.5f', p196a);
ylabel('miR-196a-5p')
title(lgd196a)

miR141 = data1(114,:)
PmiR141 = miR141(1:13)
PmiR141(14) = NaN
PmiR141(15) = NaN
HmiR141 = miR141(14:28)
p141 = ranksum(PmiR141,HmiR141)
figure, boxplot([PmiR141',HmiR141'],'Notch','off','Labels',{'Patient ','Control'},'Whisker',1,'PlotStyle','compact')
lgd141 = sprintf('p = %.5f', p141);
ylabel('miR-141-3p')
title(lgd141)

miR577 = data1(590,:)
PmiR577 = miR577(1:13)
PmiR577(14) = NaN
PmiR577(15) = NaN
HmiR577 = miR577(14:28)
p577 = ranksum(PmiR577,HmiR577)
figure, boxplot([PmiR577',HmiR577'],'Notch','off','Labels',{'Patient ','Control'},'Whisker',1,'PlotStyle','compact')
lgd577 = sprintf('p = %.5f', p577);
ylabel('miR-577')
title(lgd577)

miR200a = data1(205,:)
PmiR200a = miR200a(1:13)
PmiR200a(14) = NaN
PmiR200a(15) = NaN
HmiR200a = miR200a(14:28)
p200a = ranksum(PmiR200a,HmiR200a)
figure, boxplot([PmiR200a',HmiR200a'],'Notch','off','Labels',{'Patient ','Control'},'Whisker',1,'PlotStyle','compact')
lgd200a = sprintf('p = %.5f', p200a);
ylabel('miR-200a-5p')
title(lgd200a)

miR1246 = data1(44,:)
PmiR1246 = miR150(1:13)
PmiR1246(14) = NaN
PmiR1246(15) = NaN
HmiR1246 = miR1246(14:28)
p1246 = ranksum(PmiR150,HmiR150)
figure, boxplot([PmiR1246',HmiR1246'],'Notch','off','Labels',{'Patient ','Control'},'Whisker',1,'PlotStyle','compact')
lgd1246 = sprintf('p = %.5f', p1246);
ylabel('miR-1246')
title(lgd1246)

miR150 = data1(135,:)
PmiR150 = miR150(1:13)
PmiR150(14) = NaN
PmiR150(15) = NaN
HmiR150 = miR150 (14:28)
p150 = ranksum(PmiR150,HmiR150)
figure, boxplot([PmiR150',HmiR150'],'Notch','off','Labels',{'Patient ','Control'},'Whisker',1,'PlotStyle','compact')
lgd150 = sprintf('p = %.5f', p150);
ylabel('miR-150-5p')
title(lgd150)

miR320 = data1(322,:)
PmiR320 = miR320(1:13)
PmiR320(14) = NaN
PmiR320(15) = NaN
HmiR320 = miR320 (14:28)
p320 = ranksum(PmiR320,HmiR320)
figure, boxplot([PmiR320',HmiR320'],'Notch','off','Labels',{'Patient ','Control'},'Whisker',1,'PlotStyle','compact')
lgd320 = sprintf('p = %.5f', p320);
ylabel('miR-320c')
title(lgd320)

miR584 = data1(592,:)
PmiR584 = miR584(1:13)
PmiR584(14) = NaN
PmiR584(15) = NaN
HmiR584 = miR584 (14:28)
p584 = ranksum(PmiR584,HmiR584)
figure, boxplot([PmiR584',HmiR584'],'Notch','off','Labels',{'Patient ','Control'},'Whisker',1,'PlotStyle','compact')
lgd584 = sprintf('p = %.5f', p584);
ylabel('miR-584-5p')
title(lgd584)

miR615 = data1(630,:)
PmiR615 = miR615(1:13)
PmiR615(14) = NaN
PmiR615(15) = NaN
HmiR615 = miR615 (14:28)
p615 = ranksum(PmiR615,HmiR615)
figure, boxplot([PmiR615',HmiR615'],'Notch','off','Labels',{'Patient ','Control'},'Whisker',1,'PlotStyle','compact')
lgd615 = sprintf('p = %.5f', p615);
ylabel('miR-615-3p')
title(lgd615)

miR744 = data1(694,:)
PmiR744 = miR744(1:13)
PmiR744(14) = NaN
PmiR744(15) = NaN
HmiR744 = miR744 (14:28)
p744 = ranksum(PmiR744,HmiR744)
figure, boxplot([PmiR744',HmiR744'],'Notch','off','Labels',{'Patient ','Control'},'Whisker',1,'PlotStyle','compact')
lgd744 = sprintf('p = %.5f', p744);
ylabel('miR-744-5p')
title(lgd744)

%-------panel eight markers healhty vs stage iv lung cancer
miR16 = data1(150,:)
PmiR16 = miR16(1:13)
PmiR16(14) = NaN
PmiR16(15) = NaN
HmiR16 = miR16 (14:28)
p16 = ranksum(PmiR16,HmiR16)
figure, boxplot([PmiR16',HmiR16'],'Notch','off','Labels',{'Patient ','Control'},'Whisker',1,'PlotStyle','compact')
lgd16 = sprintf('p = %.5f', p16);
ylabel('miR-16-5p')
title(lgd16)

miR25 = data1(254,:)
PmiR25 = miR25(1:13)
PmiR25(14) = NaN
PmiR25(15) = NaN
HmiR25 = miR25 (14:28)
p25 = ranksum(PmiR25,HmiR25)
figure, boxplot([PmiR25',HmiR25'],'Notch','off','Labels',{'Patient ','Control'},'Whisker',1,'PlotStyle','compact')
lgd25 = sprintf('p = %.5f', p25);
ylabel('miR-25-3p')
title(lgd25)

miR122 = data1(35,:)
PmiR122 = miR122(1:13)
PmiR122(14) = NaN
PmiR122(15) = NaN
HmiR122 = miR122 (14:28)
p122 = ranksum(PmiR122,HmiR122)
figure, boxplot([PmiR122',HmiR122'],'Notch','off','Labels',{'Patient ','Control'},'Whisker',1,'PlotStyle','compact')
lgd122 = sprintf('p = %.5f', p122);
ylabel('miR-122-5p')
title(lgd122)
figure,bar(PmiR122')
figure,bar (HmiR122')
figure,hist(PmiR122)
figure,hist(HmiR122)

miR143 = data1(117,:)
PmiR143 = miR143(1:13)
PmiR143(14) = NaN
PmiR143(15) = NaN
HmiR143 = miR143 (14:28)
p143 = ranksum(PmiR143,HmiR143)
figure, boxplot([PmiR143',HmiR143'],'Notch','off','Labels',{'Patient ','Control'},'Whisker',1,'PlotStyle','compact')
lgd143 = sprintf('p = %.5f', p143);
ylabel('miR-143-3p')
title(lgd143)

miR200 = data1(206,:)
PmiR200 = miR200(1:13)
PmiR200(14) = NaN
PmiR200(15) = NaN
HmiR200 = miR200 (14:28)
p200 = ranksum(PmiR200,HmiR200)
figure, boxplot([PmiR200',HmiR200'],'Notch','off','Labels',{'Patient','Control'},'Whisker',1,'PlotStyle','compact')
ylabel('miR-200b-5p')
lgd200 = sprintf('p = %.5f', p200);
title(lgd200)

miR223 = data1(240,:)
PmiR223 = miR223(1:13)
PmiR223(14) = NaN
PmiR223(15) = NaN
HmiR223 = miR223 (14:28)
p223 = ranksum(PmiR223,HmiR223)
figure, boxplot([PmiR223',HmiR223'],'Notch','off','Labels',{'Patient','Control'},'Whisker',1,'PlotStyle','compact')
ylabel('miR-233-5p')
lgd223 = sprintf('p = %.5f', p223);
title(lgd223)

miR486 = data1(506,:)
PmiR486 = miR486(1:13)
PmiR486(14) = NaN
PmiR486(15) = NaN
HmiR486 = miR486 (14:28)
p486 = ranksum(PmiR486,HmiR486)
figure, boxplot([PmiR486',HmiR486'],'Notch','off','Labels',{'Patient ','Control'},'Whisker',1,'PlotStyle','compact')
lgd486 = sprintf('p = %.5f', p486);
ylabel('miR-486-5p')
title(lgd486)
 
miR708 = data1(690,:)
PmiR708 = miR708(1:13)
PmiR708(14) = NaN
PmiR708(15) = NaN
HmiR708 = miR708 (14:28)
p708 = ranksum(PmiR708,HmiR708)
figure, boxplot([PmiR708',HmiR708'],'Notch','off','Labels',{'Patient ','Control'},'Whisker',1,'PlotStyle','compact')
lgd708 = sprintf('p = %.5f', p708);
ylabel('miR-708-5p')
title(lgd708)

% data1 = xlsread('lung_oncomirs_full947noLabelH3.xlsx');
% data1(isnan(data1))=0;
% xlswrite('lung_oncomirs_full947noLabelnoZeroH3.xlsx',data1);

%miRlog = log2(miR);
%aux = find(miRlog == -Inf);
%miRlog(aux)=0;
%xlswrite('healthy_miRlog2.xlsx',miRlog);

%QUANTILE NORMALIZATION BY MIRS (ROWS/HORIZONTAL)
miRm = [];
for i = 1:223
    aux = median(miR(i,:));
    for j = 1:45
        miRm(i,j) = miR(i,j) - aux;
    end
end
xlswrite('healthy_miRmedNor.xlsx',miRm);

%LOESS NORMALIZATION
% miR1 = miR(:);%aux = miR';%miR2 = aux(:);
% y=[1 2 3];
% z= [];y1 = [];
% for i = 1:223
%     z(i)  =mean (miR(i,:));% this has to be computed over 15 people but for every time point %y1 = [y1; y'];
% end
% z1 = [];
% for i = 1:45
%     z1 = [z1;z'];
% end
% for i = 1:3345
%             y1 = [y1; y'];
% end
% f = fit([miR1 y1],z1,'lowess')%load franke
% plot(f,[miR1 y1],z1,'Exclude',out)
% xlabel( 'miR [abundance]' );
% ylabel( 'time points 1 2 3' );
% zlabel( 'mean miR per every time point over the 15 people' );
% miR

