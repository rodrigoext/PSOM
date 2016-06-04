function [ pmatrix ] = CalculatePMatrix( data, codebook, raio )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
pmatrix = zeros(1,length(codebook));
for i = 1 :  length(codebook)
    temp = codebook(i,:);
    for j = 1 :  length(data)
        tempData = data(j,:);
        distancia = sqrt(sum((temp - tempData) .^ 2));
        if distancia <= raio
            pmatrix(1,i) = pmatrix(1,i) + 1;
        end
    end
end

end

function raio = CalculateRaio(data)

contaDistancia = 0;
percents = zeros(1,length(data));
distancias = zeros(length(data), length(data));
numRaios = (length(data) - 1) * length(data) / 2;
raios = zeros(1,length(data));

for h = 1 : length(data)
    temp1 = data(h,:);
    for i = 1 : h
        temp = data(i,:);
        contaDistancia = contaDistancia + 1;
        raio =  sqrt(sum((temp - temp1) .^ 2));
        distancias(h,i) = raio;
        raios(contaDistancia) = raio;
    end
end

rt = sortrows(raios');
raios = rt';

raio = 1.0;

end

