%% load phantom data/matrix and perform fosters, then save

%% constant variables 
Lin = 8; % Truncation order of the internal VSH basis
Lout = 3; % Truncation order of the external VSH basis
dim_in = (Lin+1)^2 - 1; % Dimension of the internal SSS basis, should be 80
origin = [0,0,0];
coordsys = 'device'; 
% specify location of infile war data and where you want to save the
% processed data
infile = 'C:/Users/xanmc/RESEARCH/Data/ILABS/Phantom/240708/phantom_32_200nam_20240708_raw.fif';
%covariance
[covariance] = load('phantom_32_200nam_20240708_cov.mat');

%setup info and load channel positions
info = fiff_read_meas_info(infile);
nchan=info.nchan;
for i=1:nchan
    R(:,i)=info.chs(i).loc(1:3,:);
    EX(:,i)=info.chs(i).loc(4:6,:);
    EY(:,i)=info.chs(i).loc(7:9,:);
    EZ(:,i)=info.chs(i).loc(10:12,:);
end

% load raw data (from matrix for now)
% [raw] = fiff_setup_read_raw(infile);
% [data,time] = fiff_read_raw_segment(raw);
[raw] = load("phantom_32_200nam_20240708.mat");
data = raw.raw_data;
times = raw.raw_times;
% mark and record stim/events data for later
stim_events = data(323,:);

%% remove all bad NaN channels
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
ni=10;
XN_it = xi([SNin,SNout],phi_0,Lin,Lout-1,ni);
data_rec_it = real(SNin*XN_it(1:size(SNin,2),:));


%% reconstruct data fosters
% load in covariance 
N = covariance.covar;
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

var_raw = var(data_rec);
var_f = var(data_rec_fosters);

%% save files
% this creates new data matricies that are the original size nchan with
% zeros for the entries of "bad" channels and adds the recon data for not
% bad entries
% preallocate arrays of zeros
data_rec_vsh = zeros(nchan,size(data_rec,2));
data_rec_fos = zeros(nchan,size(data_rec_fosters,2));
data_rec_its = zeros(nchan,size(data_rec_it,2));
% add back in good data
k=1;
for i=1:nchan
    if ismember(i, bad_chans)
        data_rec_vsh(i,:) = data_rec_vsh(i,:);
        data_rec_fos(i,:) = data_rec_fos(i,:);
        data_rec_its(i,:) = data_rec_its(i,:);
    else
        data_rec_vsh(i,:)=data_rec(k,:);
        data_rec_fos(i,:)=data_rec_fosters(k,:);
        data_rec_its(i,:)=data_rec_it(k,:);
        k=k+1;
    end
end
% add back in the "stim_events" channel data for python processing
data_rec_vsh(323,:)=stim_events;
data_rec_fos(323,:)=stim_events;
data_rec_its(323,:)=stim_events;

%% save!!
% make sure to change this to "iterative" or "sss" instead of "fost"
% depending on which data matrix you are saving
outfile = 'C:/Users/xanmc/RESEARCH/Data/fosters_processed/phantom_32_200nam_20240708_fost.fif';
[outfid,cals] = fiff_start_writing_raw(outfile,info);
fiff_write_raw_buffer(outfid,data_rec_fos,cals);
fiff_finish_writing_raw(outfid);
