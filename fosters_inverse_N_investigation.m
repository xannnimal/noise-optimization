%% implementing Foster's inverse with simulated data
% noise optimization of basic component extraction
% this is for simulated data
% edits since 6/62025
clear
%% constant variables 
Lin = 8; % Truncation order of the internal VSH basis
Lout = 3; % Truncation order of the external VSH basis
dim_in = (Lin+1)^2 - 1; % Dimension of the internal SSS basis, should be 80

%% generate SQUID magnetometers
% load in positions from file, can be found on MNE-Python
coordsys = 'device'; 
rawfile = 'sample_audvis_raw.fif';
[R,EX,EY,EZ] = fiff_getpos(rawfile, coordsys);
info = fiff_read_meas_info(rawfile);
grad = ft_read_sens(rawfile, 'coordsys', 'dewar', 'senstype', 'meg', 'coilaccuracy', 2); % with coilaccuracy being 0, 1 or 2.
EZ=grad.chanori';
R=grad.chanpos';
for i=(1:size(EZ,2))
    if mod(i,3)==0 %every third is a magnetometer
        ch_types(i)=1;
    else
        ch_types(i)=0;
    end
end
k=1;
for i=(1:306)
    if ch_types(i)==1 %every third is a magnetometer
        mags(k)=i;
        k=k+1;
    else
        k=k;
    end
end

%% generate current dipole
%current dipole using Samu's implementation of Sarvas
rs=[0,0,0];
q=[0,1,0]; %y direction
r0=[0.05,0,0]; %5cm along x axis
f_start = 100; % start frequency
f_end = 50; % end frequency
f_start_out = 50; % start frequency
f_end_out = 30; % end frequency
timestep = 0.0001;
T = 0.05;
rate_of_change = (f_start - f_end)/T;
rate_of_change_out=(f_start_out-f_end_out)/T;
times = timestep:timestep:T;
for i=(1:3)
    q_t(i,:) = q(i)*sin(2*pi*(f_start*times - times.^2*rate_of_change/2));
end

% interior current dipole, calculate multipole moments
scale = 1e-12;
phi_0 = zeros(306,500);
for i=(1:size(times/2,2))
    phi_0(:,i) = dipole_field_sarvas(rs',q_t(:,i),r0',R,EX,EY,EZ,mags)*scale';
    alpha(:,i) = alpha_dipole(r0',q_t(:,i),Lin)*scale'; %not normalized, in SI units
end
%add gaussian noise at 10 percent of max value of phi_0
rng("default") %set seed
noise = 0.5*randn(size(phi_0,1),2*size(phi_0,2))*scale;
% add the dipole times after the baseline of only noise
phi_0=[zeros(size(phi_0,1),size(phi_0,2)),phi_0];
phi_0=phi_0+noise;

for i=(1:size(phi_0,1))
    if mod(i,3)==0 %every third is a magnetometer
        phi_0(i,:)=phi_0(i,:)*100;
    else
        phi_0(i,:)=phi_0(i,:);
    end
end


%% tradition reconstruction
[Sin,SNin] = Sin_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,Lin);
[Sout,SNout] = Sout_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,Lout);
S = [SNin]; %change to [SNin SNout] for full basis
pS=pinv(S); %95x306
XN=pS*phi_0; %multiple moments 95x500
data_rec=real(SNin*XN(1:size(SNin,2),:));

%% iterative recon
% ni=10;
% XN_it = xi([SNin,SNout],phi_0,Lin,Lout-1,ni);
% data_rec_it = real(SNin*XN_it(1:size(SNin,2),:));

%% reconstruct data fosters
% try different N matrices
% load in covariance 
% covariance = load("covariance.mat","cov");
% N = covariance.cov;
N = diag(cov(noise'));
%diagonal elements, average of N
id_15 = diag(eye(size(EZ,2),size(EZ,2))*1e-20);
id_7 = eye(size(EZ,2),size(EZ,2))*1e-20;
id_6 = eye(size(EZ,2),size(EZ,2))*1e-18;
id_5 = eye(size(EZ,2),size(EZ,2))*1e-15;
N_mean = diag(eye(size(EZ,2),size(EZ,2))*mean(mean(N)));
%N_zeros = zeros([size(EZ,2),size(EZ,2)]);

data_rec_15 = sim_fosters_inverse(R,EX,EY,EZ,ch_types,phi_0,timestep,alpha,id_15);
data_rec_7 = sim_fosters_inverse(R,EX,EY,EZ,ch_types,phi_0,timestep,alpha,id_7);
data_rec_6 = sim_fosters_inverse(R,EX,EY,EZ,ch_types,phi_0,timestep,alpha,id_6);
data_rec_5 = sim_fosters_inverse(R,EX,EY,EZ,ch_types,phi_0,timestep,alpha,id_5);
data_rec_fosters = sim_fosters_inverse(R,EX,EY,EZ,ch_types,phi_0,timestep,alpha,N);
data_rec_mean = sim_fosters_inverse(R,EX,EY,EZ,ch_types,phi_0,timestep,alpha,N_mean);

% var_raw = var(data_rec);
% var_f = var(data_rec_fosters);


%% sub angles
% will be much lower if angle is calculated for just times with signal, not
% baseline
% for i=(501:1000) % only do signal times, not baseline noise
%     check_data(i) = subspace(data_rec(:,i), phi_0)*180/pi;
%     %check_data_it(i) = subspace(phi_0(:,i), data_rec_it)*180/pi;
%     check_data_fosters(i) = subspace(data_rec_fosters(:,i), phi_0)*180/pi;
% end
% check_data_min = min(check_data);
% check_data_max = max(check_data);
% check_data_mean = mean(check_data);
% 
% % check_data_it_min = min(check_data_it);
% % check_data_it_max = max(check_data_it);
% % check_data_it_mean = mean(check_data_it);
% 
% check_data_fosters_min = min(check_data_fosters);
% check_data_fosters_max = max(check_data_fosters);
% check_data_fosters_mean = mean(check_data_fosters);


%% plot raw and recon data
chan=3  ;
figure(1)
hold on
times2 = timestep:timestep:timestep*size(phi_0,2);
plot(times2,phi_0(chan,:),'linewidth',2) %raw
%plot(times2,data_rec(chan,:),"c") %sss
%plot(times2,data_rec_fosters(chan,:),'linewidth',2) %foster
plot(times2,data_rec_mean(chan,:),'linewidth',2) 
plot(times2,data_rec_15(chan,:),'linewidth',2) 
%plot(times2,data_rec_7(chan,:),'linewidth',2) 
%plot(times2,data_rec_6(chan,:),'linewidth',2) 
%plot(times2,data_rec_5(chan,:),'linewidth',2) 
title('306 SQUID, Currrent Dipole [5cm,0,0], Varying N Reconstruction')
xlabel('Time (sec)')
ylabel('Dipole Signal, Chan 3 (T)')
%legend({'raw','SSS','iter','Fosters'},'location','northwest')
legend({'raw','N=2.49e-26',"N=1e-20"},'location','northwest')
%legend({'raw','Fosters'},'location','northwest')


%% SNR calculations (signal/noise)
% noise from baseline, peak from data
% if A is a matrix whose columns are random variables and whose rows are observations 
% then S is a row vector containing the standard deviation corresponding to each column
n1=40; n2=200;
s1=640; s2=800;
%peak/noise
SNR_raw = mean(std(phi_0(:,s1:s2)))/mean(std(phi_0(:,n1:n2)));
SNR_vsh = mean(std(data_rec(:,s1:s2)))/mean(std(data_rec(:,n1:n2)));
%SNR_it = mean(std(data_rec_it(:,s1:s2)))/mean(std(data_rec_it(:,n1:n2)));
SNR_fost = mean(std(data_rec_fosters(:,s1:s2)))/mean(std(data_rec_fosters(:,n1:n2)));
SNR_15 = mean(std(data_rec_15(:,s1:s2)))/mean(std(data_rec_15(:,n1:n2)));
SNR_7 = mean(std(data_rec_7(:,s1:s2)))/mean(std(data_rec_7(:,n1:n2)));
SNR_6 = mean(std(data_rec_6(:,s1:s2)))/mean(std(data_rec_6(:,n1:n2)));
SNR_5 = mean(std(data_rec_5(:,s1:s2)))/mean(std(data_rec_5(:,n1:n2)));


%% calculate PSD - something like this
%estimate of PSD. "periofogram" function has lots of variations we should
% investigate, this is the simplest implementation for now
[psd_raw,f]=periodogram(phi_0);
[psd_sss,f]=periodogram(data_rec);
%[psd_it,f]=periodogram(data_rec_it);
[psd_fosters,f]=periodogram(data_rec_fosters);
[psd_z,f]=periodogram(data_rec_mean);
[psd_15,f]=periodogram(data_rec_15);
[psd_7,f]=periodogram(data_rec_7);
[psd_6,f]=periodogram(data_rec_6);
[psd_5,f]=periodogram(data_rec_5);

%plot median PSD
figure(2)
hold on
%semilogy(f,median(psd_raw,2),'linewidth',2)
%semilogy(f,median(psd_sss,2),'linewidth',2)
%semilogy(f,median(psd_it,2),'linewidth',2)
%semilogy(f,median(psd_fosters,2),'linewidth',2)
% semilogy(f,median(psd_z,2),'linewidth',2)
semilogy(f,median(psd_15,2),'linewidth',2)
% semilogy(f,median(psd_7,2),'linewidth',2)
% semilogy(f,median(psd_6,2),'linewidth',2)
% semilogy(f,median(psd_5,2),'linewidth',2)
% grid on
legend('Raw','N=Mean','N=1e-20' ,'Location','best')
xlabel('frequency (Hz)')
ylabel('noise (fT/rtHz)')
title('PSD of Reconstructed Simulated Data')
set(gca,'FontSize',12)
set(gcf,'color','w')
