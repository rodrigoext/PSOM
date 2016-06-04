load('umat.csv');
load('watershed.csv');
load('watershed_pmatrix.csv');
load('pmatrix.csv');
load('codebook.csv');

figure;surface(umat);
figure;surface(pmatrix);
figure;imagesc(watershed);
figure;imagesc(watershed_pmatrix);