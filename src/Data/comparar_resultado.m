% data_t = dataseis';
% data_tn = mapminmax(data_t);
% D = data_tn';
sm = som_make(D);
try
    figure;som_show(sm);
catch
    %
end
figure;plot(D(:,1), D(:,2), '.b');
hold on;
plot(sm.codebook(:,1), sm.codebook(:,2), '*r');
sm.codebook = rc;
try
    figure;som_show(sm);
catch
    %
end
figure;plot(D(:,1), D(:,2), '.b');
hold on;
plot(sm.codebook(:,1), sm.codebook(:,2), '*r');