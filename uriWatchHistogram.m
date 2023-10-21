function  mirs = uriWatchTest1

mirs = xlsread('E:\Company_Paperwork\JOHAN\Data master.xlsx');
% mirs = xlsread('D:\Desktop10oct2018\JOHAN\Data master.xlsx');

%run mirs=uriWatchTest1;

% figure,plot(mirs(1,:))
% hold on,plot(mirs(2,:),'red')
% hold on,plot(mirs(3,:),'green')
% hold on,plot(mirs(4,:),'yellow')

mirNew = mirs(:,[1,6,11,16, 21, 26, 31, 36, 41, 46, 51, 56, 61, 66, 71, 76]);
% mirNew = mirs(:,[2,7,12,17, 22, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 77]);

mPark = mirNew([79, 547],:)
x = mPark(1,:);
xx = log2(x);
figure
y= linspace(0,1,16);
s = scatter(y,xx,'red');

% mir146a-5p, mir-106b-3p,miR-195-5p, miR-20b-5p, 
% miR-455-3p, mir-29c-3p, mir-93-5p, miR-19b-3p,miR-501-3p, 
% mir-486-5p, mir-483-5p, mir-502-3p, mir-200a-3p, miR-29a
% mir-378a-3p, mir-1291, mir-143-3p, miR-142-3p, miR-328-3p,
% miR-193a-5p, miR-30a-3p, miR-19b-3p, miR-30d-5p, miR-340-5p
% miR-140-5p, hsa-miR-125b-5p, miR-26b-5p, miR-16-5p
% mir-146a-5p, miR-29a-3p, miR-15b-5p, mir-223-3p
mAlzh = mirNew([46 82 535 256 253 197 89 196 31 4 326 129 52 88 43 324 16 603 323 121 32 196 15 125 462 45 66 47 46 88 73 104],:) 

figure,plot(mAlzh(1,:))
hold on,plot(mAlzh(2,:),'red')
hold on,plot(mAlzh(3,:),'green')
hold on,plot(mAlzh(4,:),'yellow')

% Cyl468 = ismember(Cylinders,[4 6 8]);
% parallelcoords(X(Cyl468,:), 'group',Cylinders(Cyl468), ...
%                'standardize','on', 'labels',varNames)
%            
% parall  elcoords(X(Cyl468,:), 'group',Cylinders(Cyl468), ...
%                'standardize','on', 'labels',varNames, 'quantile',.25)
mAlzh