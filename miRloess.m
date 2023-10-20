function miRloess


%miR = xlsread('healthy_mircount.xlsx');

maStruct = gprread('mouse_a1wt.gpr');
cy5data = magetfield(maStruct, 'F635 Median');
cy3data = magetfield(maStruct, 'F532 Median');

[x,y] = mairplot(cy5data, cy3data);
drawnow
ysmooth = malowess(x,y);
hold on;
plot(x, ysmooth, 'rx')
ynorm = y - ysmooth;