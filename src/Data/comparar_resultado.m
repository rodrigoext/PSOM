% data_t = dataseis';
% data_tn = mapminmax(data_t);
% D = data_tn';
sd = (mapminmax(simple'))';
sm = som_make(sd);
try
    figure;som_show(sm);
catch
    %
end
figure;plot(sd(:,1), sd(:,2), '.b');
hold on;
plot(sm.codebook(:,1), sm.codebook(:,2), '*r');
sm.codebook = rc;
try
    figure;som_show(sm);
catch
    %
end
figure;plot(sd(:,1), sd(:,2), '.b');
hold on;
plot(sm.codebook(:,1), sm.codebook(:,2), '*r');