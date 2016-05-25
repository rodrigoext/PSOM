codebook = load('codebook.csv');
uc = load('umat.csv');
sm.codebook = codebook;
umat = som_umat(sm);
figure;imagesc(uc);
figure;imagesc(umat);