R1=0
C1=1
R2=10
C2=5

M = csvread('model_s50_w2_generif10-50K.2-gram/model_s50_w2_g25_generif10-50K.2-gram.txt',R1,C1);

M(1:10)
size(M)

k=10

idx = kmeans(M,k,'Start','plus','Display', 'final','Distance','cosine', 'EmptyAction', 'drop','Replicates',5);

dlmwrite('model_s50_w2_g25_generif10-50K.2-gram-cluster.txt',idx)
size(idx)
%hist(idx)

