function  health = HealthyThreePointMirsLung



%mirs = xlsread('D:\Desktop10oct2018\JOHAN\UriHealthy1.xlsx');
mirs = xlsread('C:\Users\dsa\Documents\MATLAB\healthy_miRmedNor_full947.xlsx');



mir_21_3p=mirs(244,:);
mir_21_5p=mirs(245,:);
mir_140_3p= mirs(130,:);
mir_140_5p= mirs(131,:);
mir_155=mirs(164,:);
mir_200b_3p=mirs(231,:);
mir_200b_5p=mirs(232,:);
mir_223_3p=mirs(270,:);
mir_223_5p=mirs(271,:);
mir_221_3p=mirs(266,:);
mir_221_5p=mirs(267,:);
mir_145_3p= mirs(140,:);
mir_145_5p= mirs(141,:);
mir_150_3p=mirs(154,:);
mir_150_5p=mirs(155,:);
mir_200a_3p=mirs(229,:);
mir_200a_5p=mirs(230,:);
mir_205_3p=mirs(238,:);
mir_205_5p=mirs(239,:);
mir_210_3p=mirs(246,:);
mir_210_5p=mirs(247,:);
mir_339_3p=mirs(394,:);
mir_339_5p=mirs(395,:);
mir_93_3p=mirs(929,:);
mir_93_5p=mirs(930,:);

mirsMM(1,:)=mir_21_3p;
mirsMM(2,:)=mir_21_5p;
mirsMM(3,:)=mir_140_3p;
mirsMM(4,:)=mir_140_5p;
mirsMM(5,:)=mir_155;
mirsMM(6,:)=mir_200b_3p;
mirsMM(7,:)=mir_200b_5p;
mirsMM(8,:)=mir_223_3p;
mirsMM(9,:)=mir_223_5p;
mirsMM(10,:)=mir_221_3p;
mirsMM(11,:)=mir_221_5p;
mirsMM(12,:)=mir_145_3p;
mirsMM(13,:)=mir_145_5p;
mirsMM(14,:)=mir_150_3p;
mirsMM(15,:)=mir_150_5p;
mirsMM(16,:)=mir_200a_3p;
mirsMM(17,:)=mir_200a_5p;
mirsMM(18,:)=mir_205_3p;
mirsMM(19,:)=mir_205_5p;
misMM(20,:)=mir_210_3p;
mirsMM(21,:)=mir_210_5p;
mirsMM(22,:)=mir_339_3p;
mirsMM(23,:)=mir_339_5p;
mirsMM(24,:)=mir_93_3p;
mirsMM(25,:)=mir_93_5p;

% M = csvread('Data_1_12_2017.csv');
% M = importdata('Data_1_12_2017.txt');
figure
x1 = mirsMM(:,9)'; %PICK JUST ONE OF THE 15 HEALTHY in TIME POINT 1
xx1 = x1;
%xx1= xx1-mean(xx1);
y= linspace(0,1,25);
s1 = scatter(y,xx1,'red' ,'LineWidth',1.5);
hold on
x2 = mirsMM(:,24)'; %PICK the SAME HEALTHY in TIME POINT 2
xx2 = x2;
%xx2= xx2-mean(xx2);
s2 = scatter(y,xx2,'green','d', 'LineWidth',1.5);
hold on
x3 = mirsMM(:,39)'; %PICK the SAME HEALTHY in TIME POINT 3
xx3 = x3;
%xx3= xx3-mean(xx3);
s3 = scatter(y,xx3,'blue','s', 'LineWidth',1.4);
legend('time point 1','time point 2','time point 3')
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