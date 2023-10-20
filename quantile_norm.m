function [quant]=quantile_norm(data,wiederholung)
%Usage
%Example
%normalized_data=quantile_norm(data,replicates)
%%data: is an 2 dimensional array that contains data that need to be
%normalized
%replicates: in case each datapoint represents and independent measurement
%and no replicate measurement than set this set this parameter to 1. 
%In case you have replicate measurements and you want them to be averaged
%before the normalization as a data reduction procedure than supply the
%number of replicates in the dataset. 
if nargin == 0
    %data = xlsread('healthy_mircount_full947noZero.xlsx');
 %[data,txt] = xlsread('lung_oncomirs_743H2.xlsx');
       %[data,txt] = xlsread('full28participants_full947.xlsx');
   %[data,txt] = xlsread('lung_oncomirs_742H1.xlsx');
      %[data,txt] = xlsread('lung_oncomirs_722H3.xlsx');
%      [data,txt] = xlsread('Alzheimer.csv');
           % [data,txt] = xlsread('Parkinson.csv');
            [data,txt] = xlsread('Control.csv');


%     m=0;% new matrix miR counter
%     indx = [];
%     ndata = [];% trimmed matrix (less miRs)
%     for i = 1:size(data1,1)% for all patients
%         
%         for j = 1:size(data1,2)% for all miRs
%             if data1(i,j)>=30 & all(data1(i,:)<300)
%                 % if one of the miRs has a value betweem 2^7 and 2^12 -
%                 % worked well for the initial range
%                 m=m+1;
%                 indx = [indx,txt(i+1,1)];
%                 ndata(m,:)=data1(j,:);
%                 break
%             end
%         end
%     end
    %m
    %indx
%     aux1 = find(data1 == NaN);
%     data1(aux1)=0;

%     look for first two local minima
%     data2 = log2(ndata);
%     aux = find(data2 == -Inf);
%     data2(aux)=0;
% %     figure, plot(data2');
%     
%     data = data2';
    %data = data1';
    wiederholung = 1;
end
total1=data';
%Getting the dimensions of the array
[lange breite]=size(total1);
%Preallocating memory saves time
average=ones(size(lange/wiederholung,breite));
quantile_rank=ones(size(lange));
quant=ones(size(lange/wiederholung,breite));
%Starting the clock
%tic
%hier werden erst mal die Mittelwerte fur die Rohdaten berechnet damit das
%normalisieren etwas schneller lauft.
h = waitbar(0,'Please wait  Averaging proteins in progress...');
for i=1 : breite
    waitbar(100/breite)
    dummy=0;
    for ii=1 : lange/wiederholung
        sum=0;
       for iii=1 : wiederholung
            dummy=dummy+1;
            sum=sum+total1(dummy,i);
       end
       average(ii,i)=sum/wiederholung;
    end
end
close(h);
% hier wird die quantile Normalisierung durch gefuhrt und dann an das
% programm ubergeben
total1=average(:,1);
totalrfu=total1;
rank=sort(total1);
lange=length(total1);
for x=1:lange
        for xx=1 : lange
            if rank(x)==total1(xx);
                rank(x)=xx;
            end
        end
end
total1=rank;
for i=2 : breite
    total=average(:,i);
    rank=sort(total);
    for x=1:lange
        for xx=1 : lange
            if total(x)==rank(xx);
                total(x)=xx;
            end
        end
    end
    t=[total1 total];
    total1=t;
end 
rank=t;
%Hier wird der quantile Wert fur jeden Rank im normalsierten Daten feld
%berechnet.
testy=sort(average);
%for i=1 : lange
%    dummy=0;
%    for ii=1 : breite
%        dummy=dummy+testy(i,ii);
%    end 
%    quantile_rank(i)=dummy/breite;
%end     
%quantile_rank;
quantile_rank=mean(testy,2);
for x=1 :lange
       for y=1 : breite
           quant(x,y)=quantile_rank(rank(x,y));
       end
end  
%h(1) = axes('Position',[0.1 0.6 0.8 0.4]);
%mesh(data)
%h(2) = axes('Position',[0.1 0.1 0.8 0.4]);
%mesh(quant)
%stop the watch and display
%toc
%xlswrite('healthy_miRmedNor_full947.xlsx',quant');
%xlswrite('lung_oncomirsLog2QuaNor_93H2trim.xlsx',quant');
%xlswrite('full28participants_medNor_full947.xlsx');
 %xlswrite('lung_oncomirs_medNor_722H3.xlsx',quant');
 %xlswrite('medNor_Alzheimer.csv',quant');
 %xlswrite('medNor_Parkinson.csv',quant');
 xlswrite('medNor_Control.csv',quant');

%filePh = fopen('lung_oncomirsLog2QuaNor_93H2trim.txt','w');
%fprintf(filePh,'%s\n',indx{:});
%fclose(filePh);
%txtwrite('lung_oncomirsQuaNor_634H2trim.txt',indx');
quant