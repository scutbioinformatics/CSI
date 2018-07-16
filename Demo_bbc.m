clear;
addpath('cluster_measures')
addpath('tools')
load('bbc.mat'); 

gnd = truth;
c = length(unique(gnd));
m=3;
X{1} = normalize(X{1});
X{2} = normalize(X{2});
X{3} = normalize(X{3});

lambda = 0.3;
beta = 0.2;
alpha = 0.3;

repeat = 1;

for t=1:repeat
    
[A,Z,err]=CSI(X,c,lambda,alpha,beta);

X_concatFea = [X{1};X{2};X{3}];
X_concatPCA = DataProjection(X_concatFea, 100);



res_mv   = SpectralClustering(A,c);
for i = 1:m
    res_z{i} = SpectralClustering(Z{i},c);
    res_v{i} = SpectralClustering3(X{i},c);
end


res_fea  = SpectralClustering3(X_concatFea,c);
res_pca  = SpectralClustering3(X_concatPCA,c);

for i = 1:m
    nmi_z{i} = nmiMeasure(gnd,res_z{i});
    nmi_v{i} = nmiMeasure(gnd,res_v{i});
    acc_z{i} = accuracyMeasure(gnd,res_z{i});
    acc_v{i} = accuracyMeasure(gnd,res_v{i});
    ari_z{i} = adjrandMeasure(gnd,res_z{i});
    ari_v{i} = adjrandMeasure(gnd,res_v{i});
    purity_z{i} = purityMeasure(gnd,res_z{i});
    purity_v{i} = purityMeasure(gnd,res_v{i});
    [f_z{i}, p_z{i}, r_z{i}] = fprMeasure(gnd,res_z{i});
    [f_v{i}, p_v{i}, r_v{i}] = fprMeasure(gnd,res_v{i});
end
nmi_mv = nmiMeasure(gnd,res_mv);
nmi_fea = nmiMeasure(gnd,res_fea);
nmi_pca = nmiMeasure(gnd,res_pca);

acc_mv = accuracyMeasure(gnd,res_mv);
acc_fea = accuracyMeasure(gnd,res_fea);
acc_pca = accuracyMeasure(gnd,res_pca);


ari_mv = adjrandMeasure(gnd,res_mv);
ari_fea = adjrandMeasure(gnd,res_fea);
ari_pca = adjrandMeasure(gnd,res_pca);

purity_mv = purityMeasure(gnd,res_mv);
purity_fea = purityMeasure(gnd,res_fea);
purity_pca = purityMeasure(gnd,res_pca);

[f_mv, p_mv, r_mv] = fprMeasure(gnd,res_mv);
[f_fea, p_fea, r_fea] = fprMeasure(gnd,res_fea);
[f_pca, p_pca, r_pca] = fprMeasure(gnd,res_pca);


nmi_csi(t)=nmi_mv;
acc_csi(t)=acc_mv;
ari_csi(t)=ari_mv;
pur_csi(t)=purity_mv;
f_csi(t)  =f_mv;
p_csi(t)  =p_mv;
r_csi(t)  =r_mv;

nmi_v1(t)=nmi_v{1};
acc_v1(t)=acc_v{1};
ari_v1(t)=ari_v{1};
pur_v1(t)=purity_v{1};
f_v1(t)  =f_v{1};
p_v1(t)  =p_v{1};
r_v1(t)  =r_v{1};


nmi_v2(t)=nmi_v{2};
acc_v2(t)=acc_v{2};
ari_v2(t)=ari_v{2};
pur_v2(t)=purity_v{2};
f_v2(t)  =f_v{2};
p_v2(t)  =p_v{2};
r_v2(t)  =r_v{2};

nmi_v3(t)=nmi_v{3};
acc_v3(t)=acc_v{3};
ari_v3(t)=ari_v{3};
pur_v3(t)=purity_v{3};
f_v3(t)  =f_v{3};
p_v3(t)  =p_v{3};
r_v3(t)  =r_v{3};


nmi_cf(t)=nmi_fea;
acc_cf(t)=acc_fea;
ari_cf(t)=ari_fea;
pur_cf(t)=purity_fea;
f_cf(t)  =f_fea;
p_cf(t)  =p_fea;
r_cf(t)  =r_fea;

nmi_cp(t)=nmi_pca;
acc_cp(t)=acc_pca;
ari_cp(t)=ari_pca;
pur_cp(t)=purity_pca;
f_cp(t)=f_pca;
p_cp(t)=p_pca;
r_cp(t)=r_pca;

end
fprintf('Method    NMI     ACC      ARI        F         P        R     Purity\n');
fprintf('MV       %.3f    %.3f    %.3f    %.3f    %.3f    %.3f    %.3f \n',mean(nmi_csi), mean(acc_csi), mean(ari_csi),mean(f_csi),mean(p_csi),mean(r_csi),mean(pur_csi));
fprintf('std      %.3f    %.3f    %.3f    %.3f    %.3f    %.3f    %.3f \n',std(nmi_csi), std(acc_csi), std(ari_csi),std(f_csi),std(p_csi),std(r_csi),std(pur_csi));
fprintf('V1       %.3f    %.3f    %.3f    %.3f    %.3f    %.3f    %.3f \n',mean(nmi_v1), mean(acc_v1), mean(ari_v1),mean(f_v1),mean(p_v1),mean(r_v1),mean(pur_v1));
fprintf('std      %.3f    %.3f    %.3f    %.3f    %.3f    %.3f    %.3f \n',std(nmi_v1), std(acc_v1), std(ari_v1),std(f_v1),std(p_v1),std(r_v1),std(pur_v1));
fprintf('V2       %.3f    %.3f    %.3f    %.3f    %.3f    %.3f    %.3f \n',mean(nmi_v2), mean(acc_v2), mean(ari_v2),mean(f_v2),mean(p_v2),mean(r_v2),mean(pur_v2));
fprintf('std      %.3f    %.3f    %.3f    %.3f    %.3f    %.3f    %.3f \n',std(nmi_v2), std(acc_v2), std(ari_v2),std(f_v2),std(p_v2),std(r_v2),std(pur_v2));
fprintf('V3       %.3f    %.3f    %.3f    %.3f    %.3f    %.3f    %.3f \n',mean(nmi_v3), mean(acc_v3), mean(ari_v3),mean(f_v3),mean(p_v3),mean(r_v3),mean(pur_v3));
fprintf('std      %.3f    %.3f    %.3f    %.3f    %.3f    %.3f    %.3f \n',std(nmi_v3), std(acc_v3), std(ari_v3),std(f_v3),std(p_v3),std(r_v3),std(pur_v3));
fprintf('fea      %.3f    %.3f    %.3f    %.3f    %.3f    %.3f    %.3f \n',mean(nmi_cf), mean(acc_cf), mean(ari_cf),mean(f_cf),mean(p_cf),mean(r_cf),mean(pur_cf));
fprintf('std      %.3f    %.3f    %.3f    %.3f    %.3f    %.3f    %.3f \n',std(nmi_cf), std(acc_cf), std(ari_cf),std(f_cf),std(p_cf),std(r_cf),std(pur_cf));
fprintf('pca      %.3f    %.3f    %.3f    %.3f    %.3f    %.3f    %.3f \n',mean(nmi_cp), mean(acc_cp), mean(ari_cp),mean(f_cp),mean(p_cp),mean(r_cp),mean(pur_cp));
fprintf('std      %.3f    %.3f    %.3f    %.3f    %.3f    %.3f    %.3f \n',std(nmi_cp), std(acc_cp), std(ari_cp),std(f_cp),std(p_cp),std(r_cp),std(pur_cp));