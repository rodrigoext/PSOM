load('umat.csv');
load('watershed.csv');
load('watershed_filter.csv');
load('umatfilter.csv');
load('pmatrix.csv');
load('codebook.csv');

figure;surface(umat);
title('Umat');
figure;surface(pmatrix);
title('PMatrix');
figure;imagesc(watershed);
title('Watershed');
figure;imagesc(watershed_filter);
title('Watershed filter');
figure;surface(umatfilter);
title('Umat filter');