%% load phantom data/matrix and perform fosters, then save

%% constant variables 
Lin = 8; % Truncation order of the internal VSH basis
Lout = 3; % Truncation order of the external VSH basis
dim_in = (Lin+1)^2 - 1; % Dimension of the internal SSS basis, should be 80
origin = [0,0,0];
coordsys = 'device'; 
% specify location of infile war data and where you want to save the
% processed data
%infile = 'C:/Users/xanmc/RESEARCH/Data/ILABS/Phantom/231207/phantom_1000nam_default_IASoff_raw.fif';
infile = 'C:/Users/xanmc/RESEARCH/Data/ILABS/Phantom/240708/phantom_32_200nam_20240708_raw.fif';
%covariance
%[covariance] = load('1000nam_default_IASoff_raw_cov.mat');
covariance = load('covariance.mat');

%setup info and load channel positions
info = fiff_read_meas_info(infile);
nchan=info.nchan;
for i=1:nchan
    R(:,i)=info.chs(i).loc(1:3,:);
    EX(:,i)=info.chs(i).loc(4:6,:);
    EY(:,i)=info.chs(i).loc(7:9,:);
    EZ(:,i)=info.chs(i).loc(10:12,:);
end

%% read full data
[raw] = fiff_setup_read_raw(infile,1);
[data,times] = fiff_read_raw_segment(raw);


%% read cropped data
% load raw data (from matrix for now)
% [raw] = fiff_setup_read_raw(infile);
% [data,time] = fiff_read_raw_segment(raw);
% [raw] = load("phantom_32_200nam_20240708.mat");
% data = raw.raw_data;
% times = raw.raw_times;

%% Deal with Bad channels
% we only want channels beginning with 'MEG...'
bad_chans = [];
k=1;
for i = 1:nchan
    if contains(info.ch_names{1,i},'MEG') 
        k=k;
    else
        bad_chans(k) = i;
        k=k+1;
    end
end 
% mark and record stim/hpi/events data for later
keeps = data(bad_chans,:);
% remove bad channels from data for now
data(bad_chans,:)=[];
R(:,bad_chans)=[];
EX(:,bad_chans)=[];
EY(:,bad_chans)=[];
EZ(:,bad_chans)=[];

ch_types=ones(size(EZ,2),1);
phi_0=data;


%% tradition reconstruction
[Sin,SNin] = Sin_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,Lin);
[Sout,SNout] = Sout_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,Lout);
S = [SNin]; 
pS=pinv(S); 
XN=pS*phi_0; 
data_rec=real(SNin*XN(1:size(SNin,2),:));

%% iterative recon
% ni=10;
% XN_it = xi([SNin,SNout],phi_0,Lin,Lout-1,ni);
% data_rec_it = real(SNin*XN_it(1:size(SNin,2),:));


%% reconstruct data fosters
% load in covariance 
%N = covariance.covar;
N = covariance.cov;
%calculate matrix "Sij" for Foster's inverse:  x(bar) = B*phi_0 + b
% alpha is now the estimate from SSS for multipole moments - double check
alpha_cov = cov(XN'); 
alpha = XN';
% for i=(1:size(times/2,2))
%     alpha(:,i) = alpha_dipole(r0',q_t(:,i),Lin)'; %not normalized, in SI units
% end

%normalize alpha_cov
for i=(1:size(Sin,2))
    for j=(1:size(Sin,2))
        alpha_cov_new(i,j)=alpha_cov(i,j)*norm(Sin(:,i))*norm(Sin(:,j));
    end
end
for i=(1:size(Sin,2))
    alpha_norm(:,i) = alpha(:,i)*norm(Sin(:,i));
end

%find inverse matrix B
S_star = conj(S)'; 
first = pinv(S*alpha_cov_new*S_star+N);
B = alpha_cov_new*S_star*first;
m_alpha = mean(alpha,1)'; %should be 80x1
b = m_alpha - B*S*m_alpha;
%better estimate for multipole moments
for i=(1:size(times,2))
    x_bar(:,i) = B*phi_0(:,i) +b;
end
data_rec_fosters= real(SNin*x_bar(1:size(SNin,2),:));


%% sub angles
% for i=(1:size(times,2))
%     check_data(i) = subspace(phi_0(:,i), data_rec)*180/pi;
%     check_data_it(i) = subspace(phi_0(:,i), data_rec_it)*180/pi;
%     check_data_fosters(i) = subspace(phi_0(:,i), data_rec_fosters)*180/pi;
% end
% check_data_min = min(check_data);
% check_data_max = max(check_data);
% check_data_mean = mean(check_data);
% 
% check_data_it_min = min(check_data_it);
% check_data_it_max = max(check_data_it);
% check_data_it_mean = mean(check_data_it);
% 
% check_data_fosters_min = min(check_data_fosters);
% check_data_fosters_max = max(check_data_fosters);
% check_data_fosters_mean = mean(check_data_fosters);


%% plot raw and recon data
% chan=3  ;
% figure(1)
% hold on
% times2 = timestep:timestep:timestep*size(phi_0,2);
% plot(times2,phi_0(chan,:),"green") %raw
% plot(times2,data_rec(chan,:),"red") %sss
% %plot(times2,data_rec_it(chan,:))
% plot(times2,data_rec_fosters(chan,:),"blue") %foster
% title('306 SQUID, Currrent Dipole [5cm,0,0] Reconstruction')
% xlabel('Time (sec)')
% ylabel('Dipole Signal, Chan 3 (T)')
% legend({'raw','SSS','iter','Fosters'},'location','northwest')
% legend({'raw','SSS','Fosters'},'location','northwest')
% %legend({'raw','Fosters'},'location','northwest')

%% save files
% this creates new data matricies that are the original size nchan with
% zeros for the entries of "bad" channels and adds the recon data for not
% bad entries
% preallocate arrays of zeros
data_rec_vsh = zeros(nchan,size(data_rec,2));
data_rec_fos = zeros(nchan,size(data_rec_fosters,2));
%data_rec_its = zeros(nchan,size(data_rec_it,2));
k=1; l=1;
for i=1:nchan
    if ismember(i, bad_chans) %add back in "keeps"
        data_rec_vsh(i,:) = keeps(l,:);
        data_rec_fos(i,:) = keeps(l,:);
        %data_rec_its(i,:) = keeps(l,:);
        l=l+1;
    else
        data_rec_vsh(i,:)=data_rec(k,:);
        data_rec_fos(i,:)=data_rec_fosters(k,:);
        %data_rec_its(i,:)=data_rec_it(k,:);
        k=k+1;
    end
end

%% save data matrix
save("phantom_32_1000nam_20240708_max1_sss.mat", 'data_rec_vsh','-v7.3')
save("phantom_32_1000nam_20240708_max1_fos.mat", 'data_rec_fos','-v7.3')
save("phantom_32_1000nam_20240708_max1_it.mat", 'data_rec_its','-v7.3')

return 
%% save fif files
% make sure to change this to "iterative" or "sss" instead of "fost"
% depending on which data matrix you are saving
outfile = 'C:/Users/xanmc/RESEARCH/Data/fosters_processed/phantom_32_200nam_20240708_test_raw.fif';
[outfid,cals] = fiff_start_writing_raw(outfile,info);
fiff_write_raw_buffer(outfid,data_rec_fos,cals);
fiff_finish_writing_raw(outfid);

%% results - message sepp
% localization
% psd spectra - should reduce white noise after processing 
% check amplitude/localization bias 
% visualization of evoked response - before and after 
% covar from unprocessed data in MNE-python, if it contains signal from
% external environment will probably give bias since SSS already gets rid
% of the external interference - diagonal matrix will probably not cause
% bias
% ask eric about all the different ways to calculate covariance
% different ways to estimate covariance. Idea what time response is,
% foster's inverse to optimize the covariance using simulated annealing to
% match where the phantom dipole set should be- minimize variance in
% baseline period while keeping peak amplitude intact 
