function ret = cosine(vectorfile, qg, K)
 
%%% E.g. path
%filepath = '/Users/mxd074/Documents/Research/generank/SimGeneRIF-tests/model_ovl_vs30_cw2_generif10-50K.2-gram/model_ovl_vs30_cw2_g10_generif10-50K.2-gram.txt'

filepath = vectorfile;

%Specify the starting row and column
R1=0;
C1=1;

%Get matrix
M = csvread(filepath,R1,C1);

size(M);

%Get the cosine distance matrix
D = pdist(M,'cosine');
size(D);
%Get original angular cosine values
Dt = bsxfun(@minus, D, 1);
cos = bsxfun(@times, Dt, -1);

%Make square distance matrix
Z = squareform(cos);
n=size(Z);
for i=1:1:n
    Z(i,i) = 1;
end

%Write out the matrix file
matpath = strrep(filepath, 'gram.txt', 'gram.cosine.csv');
dlmwrite(matpath,Z);

%%%%
%Get top 10-rankings for query generif
%%%%

%query generif
query_generif_idx = qg; %e.g. 321

%Get column vector of distances for this generif
qv = Z(:,query_generif_idx);

[sqv,I] = sort(qv,'descend');

%Get top-K indices
ret = cellstr(strcat(num2str(I,'%d='),num2str(sqv,'%6.5f')));

end
