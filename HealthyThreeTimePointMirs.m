function  health = HealthyThreePointMirs



mirs = xlsread('D:\Desktop10oct2018\JOHAN\UriHealthy1.xlsx');

% M = csvread('Data_1_12_2017.csv');
% M = importdata('Data_1_12_2017.txt');
figure
x1 = mirs(:,2)';
xx1 = log2(x1);
y= linspace(0,1,947);
s1 = scatter(y,xx1,'red');
hold on
x2 = mirs(:,17)';
xx2 = log2(x2);

s2 = scatter(y,xx2,'yellow');
hold on
x3 = mirs(:,32)';
xx3 = log2(x3);
y= linspace(0,1,947);
s3 = scatter(y,xx3,'green');

% scatter3(x,y,z,'filled')
% view(-30,10)
subplot(3,1,1);
scatter(y,xx1,'red');
ylim([-1 16])
subplot(3,1,2); 
scatter(y,xx2,'blue');
ylim([-1 16])
subplot(3,1,3); 
scatter(y,xx3,'green');
ylim([-1 16])


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