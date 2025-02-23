function Y = prox_group_lasso1(X,groups,tau)

nbg =size(groups,1);

[P,N] = size(X);

Y = zeros(P,N);
    for p = 1:nbg
        Y(:,groups{p,1}) = vector_soft_row(X(:,groups{p,1}),tau);
    end
