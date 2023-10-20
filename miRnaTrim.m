function miRnaTrim

% select only the unique rows(mirs) of the data

data1 = xlsread('lung_oncomirsQuaNor_full947H2.xlsx');
% tot = []; totV = [];
% for i = 1:28
%     aux = data1(:,i);
%     totV(i) = var(data1(:,i));
%     tot(i) = sum(aux);
% end

[aux1] = unique(data1, 'rows');
xlswrite('lung_oncomirsQuaNor_full618H2unique.xlsx',aux1);

data = data1(aux1, :);

length(aux1)

length(aux2)

figure,hist(data)

data
