function data_rec_fosters = sim_fosters_inverse(R,EX,EY,EZ,ch_types,phi_0,timestep,alpha,N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[Sin,SNin] = Sin_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,8);
% pS=pinv(S); %95x306
% XN=pS*phi_0; %multiple moments 95x500
alpha_cov = cov(alpha');
for i=(1:size(Sin,2))
    for j=(1:size(Sin,2))
        alpha_cov_norm(i,j)=alpha_cov(i,j)*norm(Sin(:,i))*norm(Sin(:,j));
    end
end
%find inverse matrix B
S_star = conj(SNin)'; 
first = pinv(SNin*alpha_cov_norm*S_star+N);
B = alpha_cov_norm*S_star*first;
m_alpha = mean(alpha,2);
b = m_alpha - B*SNin*m_alpha;
%better estimate for multipole moments
for i=(1:size(timestep:timestep:timestep*size(phi_0,2),2))
    x_bar(:,i) = B*phi_0(:,i) +b;
end
data_rec_fosters= real(SNin*x_bar(1:size(SNin,2),:));
end