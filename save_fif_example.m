%% trying to get a good .fif file saved to read into mne-python
%% constant variables 
clear
Lin = 8; % Truncation order of the internal VSH basis
Lout = 3; % Truncation order of the external VSH basis
dim_in = (Lin+1)^2 - 1; % Dimension of the internal SSS basis, should be 80
origin = [0,0,0];
coordsys = 'device'; 
% specify location of infile war data and where you want to save the
% processed data
infile = 'C:/Users/xanmc/RESEARCH/Data/ILABS/Phantom/240708/phantom_32_200nam_20240708_raw.fif';
outfile = 'C:/Users/xanmc/RESEARCH/Data/fosters_processed/phantom_32_200nam_20240708_fost.fif';

info = fiff_read_meas_info(infile);
nchan=info.nchan;
[fid, tree, dir] = fiff_open(infile);
[raw] = fiff_setup_read_raw(infile);
[data,time] = fiff_read_raw_segment(raw);
for i=1:nchan
    R(:,i)=info.chs(i).loc(1:3,:);
    EX(:,i)=info.chs(i).loc(4:6,:);
    EY(:,i)=info.chs(i).loc(7:9,:);
    EZ(:,i)=info.chs(i).loc(10:12,:);
end

return

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
phi_0p=data;

%% calculate SSS basis
%calculate single in single out
[Sin,SNin]= Sin_vsh_vv([0,0,0]',opm_matrix,theta_hat,phi_hat,R_hat,ch_types,Lin);
[Sout,SNout] = Sout_vsh_vv([0,0,0]',opm_matrix,theta_hat,phi_hat,R_hat,ch_types,Lout);

%Do fosters inverse

%% reconstrct internal data
%single in, single out
pS_p=pinv([SNin SNout]);
XN_p=pS_p*phi_0p;
data_rec_vsh=real(SNin*XN_p(1:size(SNin,2),:));

%% save files
data_rec_vsh = zeros(raw.info.nchan,size(data_rec_vsh_p,2));
data_rec_mvsh = zeros(raw.info.nchan,size(data_rec_multi_vsh_p,2));
data_rec_oid = zeros(raw.info.nchan,size(data_rec_sph_sph_p,2));
data_rec_oidvsh = zeros(raw.info.nchan,size(data_rec_sph_vsh_p,2));
k=1;
for i=1:95
    if ismember(i, bad_chans)
        data_rec_vsh(i,:) = data_rec_vsh(i,:);
        data_rec_mvsh(i,:) = data_rec_mvsh(i,:);
        data_rec_oid(i,:) = data_rec_oid(i,:);
        data_rec_oidvsh(i,:) = data_rec_oidvsh(i,:);
    else
        data_rec_vsh(i,:)=data_rec_vsh_p(k,:);
        data_rec_mvsh(i,:)=data_rec_multi_vsh_p(k,:);
        data_rec_oid(i,:)=data_rec_sph_sph_p(k,:);
        data_rec_oidvsh(i,:)=data_rec_sph_vsh_p(k,:);
        k=k+1;
    end

end
data_rec_oidvsh(37,:)=keep_37;
data_rec_oidvsh(95,:)=keep_95;

[outfid,cals] = fiff_start_writing_raw(outfile,raw.info);
fiff_write_raw_buffer(outfid,data_rec_oidvsh,cals);
fiff_finish_writing_raw(outfid);
