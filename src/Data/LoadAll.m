load('umat.csv');
load('watershed.csv');
load('watershed_filter.csv');
load('umatfilter.csv');
load('pmatrix.csv');
load('codebook.csv');
load('ustar.csv');
load('data.csv');
load('codebook_init.csv');
load('umultsch.csv');
load('ustarw.csv');
load('immersion.csv');
load('simulation.csv');

figure;surface(umat);
title('Umat');
figure;surface(pmatrix);
title('PMatrix');
%figure;imagesc(watershed);
%title('Watershed');
%figure;imagesc(watershed_filter);
%title('Watershed filter');
figure;surface(umatfilter);
title('Umat filter');
figure;surface(ustar);
title('UStar');
figure;surface(umultsch);
title('UMatrix Ultsch');
figure;surface(ustarw);
title('UStar Watershed');
figure;surface(immersion);
title('Immersion');

figure;plot(data(:,1), data(:,2), '*b');
hold on;
plot(codebook_init(:,1), codebook_init(:,2), '*r');
title('Codebook Initial vs Input');

figure;plot(data(:,1), data(:,2), '*b');
hold on;
plot(codebook(:,1), codebook(:,2),'*r');
title('Codebook vs Input');

figure;gscatter(data(:,1), data(:,2), simulation);

% 
% figure;plot3(data(:,1), data(:,2), data(:,3), '*b');
% hold on;
% plot3(codebook(:,1), codebook(:,2), codebook(:,3), '*r');
% title('Codebook vs Input');
