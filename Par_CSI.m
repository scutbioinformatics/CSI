function [A,Z,err]=Par_CSI(X,c,lambda,alpha,beta)
% X multi-view data, each column is a smaple
% c number of clusters
m = length(X);
n = size(X{1},2);
W0 = eye(n);
Z = cell(1,m);
for v = 1:m
    Z{v}=eye(n);
    pW =constructW(X{v}');
    pD = diag(sum(pW));
    pL{v} = pD-pW;
    W0 = W0 + pW/3;
end


A = W0;
D0 = diag(sum(W0));
L0 = D0 - W0;
P = eig1(L0,c,0);

maxIter = 10;
options.Display = 'off';

H = m*eye(n);
H = (H'+ H)/2;
d = zeros(n,n);


numCores = 4;
try
    fprintf('Closing any pools...\n');
    delete(gcp('nocreate')); 
catch
    %ignore any errors 
end
parpool('local',numCores);


for iter=1:maxIter
    Ak = A;
    % update A
    sum_Z = zeros(n,n);
    for v = 1:m
        sum_Z = sum_Z + Z{v};
    end
    
    % pre calculate d for pararell
    for i=1:n
        for j=1:n
            d(i,j)=(norm(P(i,:)-P(j,:)))^2;
        end
    end
    parfor i = 1:n
        f = -2*sum_Z(:,i)+beta*d(i,:)';
        A(:,i) = quadprog(H,f,[],[],ones(1,n),1,zeros(n,1),ones(n,1),[],options);
    end
    % update Zv
    for v = 1:m
        temp = (X{v}'*X{v}+lambda*eye(n)+alpha*pL{v})\(X{v}'*X{v}+lambda*A);
        Z{v} = temp - temp.*eye(n);
    end
    % update P
    
    A = (A'+A)/2;
    D = diag(sum(A));
    L = D-A;
    [P, ~, ~]=eig1(L, c, 0);
    
    
    error = norm(Ak-A);
    sum1 = 0; sum2 = 0; sum3 = 0; sum4 = trace(P'*L*P);
    for i = 1:m
        sum1 = sum1+norm(X{v}-X{v}*Z{v},'fro')^2;
        sum2 = sum2+norm(A-Z{v},'fro')^2;
        sum3 = sum3+trace(Z{v}*pL{v}*Z{v}');
    end
    err(iter) = sum1+lambda*sum2+alpha*sum3+beta*sum4;
    if(iter>1)
        fprintf('iter = %d, error of A = %f\n', iter, error);
    end
    if error < 1e-4 && iter > 2
        break;
    end
end

fprintf('Closing the pool\n');
delete(gcp('nocreate'))

end

