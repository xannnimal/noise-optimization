%% Load UCL OPM dataset from MNE-Python, then process with Foster's Inverse
% reconstructed data is then saved to ".fif" format to be loaded back into
% Python for topo/evoked plotting using MNE-Python functions
%% constant variables 
clear
Lin = 6; % Truncation order of the internal VSH basis
Lout = 3; % Truncation order of the external VSH basis
dim_in = (Lin+1)^2 - 1; % Dimension of the internal SSS basis, should be 80
ymin=-1e-10;
ymax=1e-10;

%% read in data
coordsys = 'device'; 
filename= 'C:/Users/xanmc/RESEARCH/UCL_OPM_example/UCL_OPM_audtitory_evoked.fif';
fileraw= 'C:/Users/xanmc/RESEARCH/UCL_OPM_example/UCL_OPM_audtitory_fix_raw.fif';
info = fiff_read_meas_info(fileraw);
nchan=info.nchan;
[raw] = fiff_setup_read_raw(fileraw);
[data,time] = fiff_read_raw_segment(raw);
for i=1:nchan
    opm_matrix(:,i)=info.chs(i).loc(1:3,:);
    theta_hat(:,i)=info.chs(i).loc(4:6,:);
    phi_hat(:,i)=info.chs(i).loc(7:9,:);
    R_hat(:,i)=info.chs(i).loc(10:12,:);
end

keep_37 = data(37,:);
keep_95 = data(95,:);

%% remove all bad NaN channels
%mark bad chans
bad_chans = [37:40,43,44,55,56,95];
raw.info.bads =[]; %list of bad channel names
k=1;
for i =1:info.nchan
    if ismember(i, bad_chans) 
        raw.info.bads{1, k} = info.ch_names{1,i};
        k=k+1;
    end
end

data(bad_chans,:)=[];
opm_matrix(:,bad_chans)=[];
theta_hat(:,bad_chans)=[];
phi_hat(:,bad_chans)=[];
R_hat(:,bad_chans)=[];

ch_types=ones(size(R_hat,2),1);
phi_0=data;

%% calculate SSS basis
% %calculate single in single out
% [Sin,SNin]= Sin_vsh_vv([0,0,0]',opm_matrix,theta_hat,phi_hat,R_hat,ch_types,Lin);
% [Sout,SNout] = Sout_vsh_vv([0,0,0]',opm_matrix,theta_hat,phi_hat,R_hat,ch_types,Lout);
% S = [SNin]; %change to [SNin SNout] for full basis
% pS=pinv(SNin); %ncomp x nchan
% XN=pS*phi_0; %multiple moments ncomp x time
% data_rec=real(SNin*XN(1:size(SNin,2),:));

%% mSSS bsais
center1= [-0.00350699, 0.01138051, 0.05947857]; 
center2= [-0.00433911, 0.04081329, 0.05194245]; 
%adjuct to device coordinate system
center1 = center1 - [0,0,0.05];
center2 = center2 - [0,0,0.05];
Lin =6; %reduce for lower channel count
thresh = 0.05;
[Sin1,SNin1] = Sin_vsh_vv(center1',opm_matrix,theta_hat,phi_hat,R_hat,ch_types,Lin);
[Sin2,SNin2] = Sin_vsh_vv(center2',opm_matrix,theta_hat,phi_hat,R_hat,ch_types,Lin); 
[UN,sigmaN,~] = svd([SNin1 SNin2]);
[U,sigma,~] = svd([SNin1 SNin2]);
sig_num = diag(sigma)';
for i=1:size(sig_num,2)
    ratio(i) = sig_num(i)/sig_num(1);
    if ratio(i) >= thresh
        SNin_tot(:,i) = UN(:,i);
        Sin(:,i)= U(:,i);
    end
end
S=SNin_tot;
pS=pinv(S);   
XN=pS*phi_0;

%% fosters
covariance = load("UCL_OPM_OTP_cov.mat","covar");
N = covariance.covar;
%N = eye(size(R_hat,2),size(R_hat,2))*mean(mean(Nt));
alpha_cov_new = cov(XN'); 
alpha = XN';
%find inverse matrix B
S_star = conj(S)'; 
first = pinv(S*alpha_cov_new*S_star+N);
B = alpha_cov_new*S_star*first;
m_alpha = mean(alpha,1)'; 
b = m_alpha - B*S*m_alpha;
%better estimate for multipole moments
for i=(1:size(time,2))
    x_bar(:,i) = B*phi_0(:,i) +b;
end
data_rec_fosters= real(S*x_bar(1:size(S,2),:));

% check_data_fos = subspace(phi_0(:,2), data_rec_fosters)*180/pi;
% check_data = subspace(phi_0(:,2), data_rec)*180/pi;

%% iterative recon
% ni=10;
% XN_it = xi([SNin,SNout],phi_0,6,2,ni);
% data_rec_fosters = real(SNin*XN_it(1:size(SNin,2),:));

% %% check subspace
% for i=(100000:300000)
%     check_data_fos(i) = subspace(phi_0(:,i), data_rec_fosters)*180/pi;
%     % check_data(i) = subspace(phi_0(:,i), data_rec)*180/pi;
% end
% 
% check_data_min_f = min(check_data_fos);
% check_data_max_f = max(check_data_fos);
% check_data_av_f = mean(check_data_fos);

% check_data_min = min(check_data);
% check_data_max = max(check_data);
% check_data_av = mean(check_data);



%% save files
% data rec is now 86xntimes because we removed bad channels, add them back
%data_rec_vsh = zeros(raw.info.nchan,size(data_rec,2));
data_rec_fos = zeros(raw.info.nchan,size(data_rec_fosters,2));
k=1;
for i=1:95
    if ismember(i, bad_chans)
        %data_rec_vsh(i,:) = data_rec_vsh(i,:);
        data_rec_fos(i,:) = data_rec_fos(i,:);
    else
        %data_rec_vsh(i,:)=data_rec(k,:);
        data_rec_fos(i,:)=data_rec_fosters(k,:);
        k=k+1;
    end

end
%add back in STIM channels
data_rec_fos(37,:)=keep_37;
data_rec_fos(95,:)=keep_95;
%save processed file
outfile = 'C:/Users/xanmc/RESEARCH/UCL_OPM_example/matlab_processed/UCL_OPM_fosters_mSSS_raw.fif';
[outfid,cals] = fiff_start_writing_raw(outfile,raw.info);
fiff_write_raw_buffer(outfid,data_rec_fos,cals);
fiff_finish_writing_raw(outfid);



