function data = test_matlab_function_foster()

%% load phantom data/matrix and perform fosters, then save
clear
%% constant variables 
Lin = 8; % Truncation order of the internal VSH basis
Lout = 3; % Truncation order of the external VSH basis
dim_in = (Lin+1)^2 - 1; % Dimension of the internal SSS basis, should be 80
origin = [0,0,0];
coordsys = 'device'; 
% specify location of infile war data and where you want to save the
% processed data
%infile = 'C:/Users/xanmc/RESEARCH/Data/ILABS/Phantom/240708/phantom_32_200nam_20240708_raw.fif';
%infile = '../Phantom/240708/phantom_32_200nam_20240708_raw.fif';
infile = 'c:/Downloads_/lab/lab_taulu/Phantom/240708/phantom_32_200nam_20240708_raw.fif';
%covariance
[covariance] = load('./raw data mats/phantom_32_200nam_20240708_cov.mat');
%[covariance] = load('./noise-optimization/covariance.mat');


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
%[raw] = load("phantom_32_200nam_20240708.mat");
[raw] = load("./raw data mats/phantom_32_200nam_20240708.mat");
data = raw.raw_data;
times = raw.raw_times;
% mark and record stim/events data for later
stim_events = data(323,:);

end