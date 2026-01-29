function residual = fosters_residual(N0,X_c,S,phantom_data,phantom_times,numdipole)
%% do fosters inverse estimate of X_f, then calculate residual with X_c
N=diag(N0,0);
X_f = cell(1, numdipole);
for k=(1:numdipole)
    alpha_cov = cov(X_c{1,k}');
    for i=(1:size(S,2))
        for j=(1:size(S,2))
            alpha_cov_new(i,j)=alpha_cov(i,j)*norm(S(:,i))*norm(S(:,j));
        end
    end 
    S_star = conj(S)'; 
    first = pinv(S*alpha_cov_new*S_star+N);
    B = alpha_cov_new*S_star*first;
    phi_0=phantom_data{1,k};
    for t=(1:size(phantom_times,2))
        X_f{k}(:,t) = B*phi_0(:,t);
    end
end
for k=(1:numdipole)
    resid(:,k) = norm(X_c{1,k}-X_f{1,k});
end
residual = mean(resid);
end