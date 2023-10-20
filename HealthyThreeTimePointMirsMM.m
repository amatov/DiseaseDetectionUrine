function  health = HealthyThreePointMirsMM



mirs = xlsread('D:\Desktop10oct2018\JOHAN\UriHealthy1.xlsx');

mir155=mirs(164,:);

mir21_3p=mirs(244,:);
mir21_5p=mirs(245,:);
mir_106b_3p= mirs(31,:);
%mir_106b_5p= mirs(32,:);
mir_181a_2_3p=mirs(174,:);
%mir_181a_3p=mirs(175,:);
mir_181a_5p=mirs(176,:);%mir_181b_3p=mirs(177,:);
mir_181b_5p=mirs(178,:);
mir_221_3p=mirs(266,:);
mir_221_5p=mirs(267,:);
mir_222_3p=mirs(268,:);
mir_16_5p=mirs(171,:);

mirsMM(1,:)=mir21_3p;
mirsMM(2,:)=mir21_5p;
mirsMM(3,:)=mir_106b_3p;
%mirsMM(4,:)=mir_106b_5p;
mirsMM(4,:)=mir_181a_2_3p;
%mirsMM(6,:)=mir_181a_3p;
mirsMM(5,:)=mir_181a_5p;
mirsMM(6,:)=mir_181b_5p;
mirsMM(7,:)=mir_221_3p;
mirsMM(8,:)=mir_221_5p;
mirsMM(9,:)=mir_222_3p;
mirsMM(10,:)=mir_16_5p;

% M = csvread('Data_1_12_2017.csv');
% M = importdata('Data_1_12_2017.txt');
figure
x1 = mirsMM(:,2)';
xx1 = log2(x1);
xx1= xx1-mean(xx1);
y= linspace(0,1,10);
s1 = scatter(y,xx1,'red');
hold on
x2 = mirsMM(:,17)';
xx2 = log2(x2);
xx2= xx2-mean(xx2);
s2 = scatter(y,xx2,'green','d');
hold on
x3 = mirsMM(:,32)';
xx3 = log2(x3);
xx3= xx3-mean(xx3);
s3 = scatter(y,xx3,'blue','filled');

% figure
% x1 = mirsMM(:,2)';
% xx1 = log2(x1);
% xx1= xx1-median(xx1);
% y= linspace(0,1,10);
% s1 = scatter(y,xx1,'red');
% hold on
% x2 = mirsMM(:,17)';
% xx2 = log2(x2);
% xx2= xx2-median(xx2);
% s2 = scatter(y,xx2,'blue');
% hold on
% x3 = mirsMM(:,32)';
% xx3 = log2(x3);
% xx3= xx3-median(xx3);
% s3 = scatter(y,xx3,'green');
% scatter3(x,y,z,'filled')
% view(-30,10)
subplot(3,1,1);
scatter(y,xx1,'red');
%ylim([-1 16])
subplot(3,1,2); 
scatter(y,xx2,'blue');
%ylim([-1 16])
subplot(3,1,3); 
scatter(y,xx3,'green');
%ylim([-1 16])


figure
h = scatter3(xx1,xx2,xx3,'filled')
h.MarkerFaceColor = [0 0.5 0.5];
view(-30,10) %view(40,35)

figure
hs(1) = subplot(2,1,1);
hs(2) = subplot(2,1,2);
scatter3(hs(1),xx1,xx2,xx3,'MarkerFaceColor',[0 .75 .75])
scatter3(hs(2),xx1,xx2,xx3,'*')
M