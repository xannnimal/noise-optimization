%% trying the "minimum infor" case of foster's inverse using noise-to-signal ratio
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
% grad = ft_read_sens(rawfile, 'coordsys', 'dewar', 'senstype', 'meg', 'coilaccuracy', 2); % with coilaccuracy being 0, 1 or 2.
% EZ=grad.chanori';
% R=grad.chanpos';
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
phi_0 = zeros(306,500);
for i=(1:size(times/2,2))
    phi_0(:,i) = dipole_field_sarvas(rs',q_t(:,i),r0',R,EX,EY,EZ,mags)';
    alpha(:,i) = alpha_dipole(r0',q_t(:,i),Lin)'; %not normalized, in SI units
end
%add gaussian noise at 10 percent of max value of phi_0
rng("default")
noise = .5*randn(size(phi_0,1),2*size(phi_0,2));
% Create an amplitude for that noise that is 10% of the noise-free signal at every element.
amplitude = 0.15 * phi_0;
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
S = [SNin,SNout]; %change to [SNin SNout] for full basis
pS=pinv(S); %95x306
XN=pS*phi_0; %multiple moments 95x500
data_rec=real(SNin*XN(1:size(SNin,2),:));

%% min info fosters inverse
% B = S*(SS* + r^2I)^dagger
% B = (S*S + r^2I)^-1 x S* this one!
n1=40; n2=200;
s1=640; s2=800;
%peak/noise
SNR_raw = mean(std(phi_0(:,s1:s2)))/mean(std(phi_0(:,n1:n2)));
ratio = (1/SNR_raw)^2;
S_star = conj(SNin)';
second = ratio*eye(size(SNin,2),size(SNin,2));
inv = pinv(S_star*SNin +second);
B=inv*S_star;
m_alpha = mean(alpha,2); %normalized alpha needed?
b = m_alpha - B*S*m_alpha;
for i=(1:size(timestep:timestep:timestep*size(phi_0,2),2))
    x_bar(:,i) = B*phi_0(:,i) +b;
end
data_rec_fosters= real(SNin*x_bar(1:size(SNin,2),:));

%% sub angles
% will be much lower if angle is calculated for just times with signal, not
% baseline
for i=(501:1000) % only do signal times, not baseline noise
    check_data(i) = subspace(data_rec(:,i), phi_0)*180/pi;
    check_data_fosters(i) = subspace(data_rec_fosters(:,i), phi_0)*180/pi;
end
check_data_min = min(check_data);
check_data_max = max(check_data);
check_data_mean = mean(check_data);

check_data_fosters_min = min(check_data_fosters);
check_data_fosters_max = max(check_data_fosters);
check_data_fosters_mean = mean(check_data_fosters);


%% plot raw and recon data
chan=3  ;
figure(1)
hold on
times2 = timestep:timestep:timestep*size(phi_0,2);
plot(times2,phi_0(chan,:),'g') %raw
plot(times2,data_rec(chan,:),'r') %sss
plot(times2,data_rec_fosters(chan,:),'b') %foster
title('306 SQUID, Currrent Dipole [5cm,0,0] Reconstruction')
xlabel('Time (sec)')
ylabel('Dipole Signal, Chan 3 (T)')
legend({'raw','SSS','Fosters'},'location','northwest')

%% SNR calculations (signal/noise)
% noise from baseline, peak from data
% if A is a matrix whose columns are random variables and whose rows are observations 
% then S is a row vector containing the standard deviation corresponding to each column
n1=40; n2=200;
s1=640; s2=800;
%peak/noise
SNR_raw = mean(std(phi_0(:,s1:s2)))/mean(std(phi_0(:,n1:n2)));
SNR_vsh = mean(std(data_rec(:,s1:s2)))/mean(std(data_rec(:,n1:n2)));
SNR_fost = mean(std(data_rec_fosters(:,s1:s2)))/mean(std(data_rec_fosters(:,n1:n2)));

%% calculate PSD - something like this
%estimate of PSD. "periofogram" function has lots of variations we should
% investigate, this is the simplest implementation for now
[psd_raw,f]=periodogram(phi_0);
[psd_sss,f]=periodogram(data_rec);
[psd_fosters,f]=periodogram(data_rec_fosters);
%plot median PSD
figure(2)
hold on
semilogy(f,median(psd_raw,2),'linewidth',2)
semilogy(f,median(psd_sss,2),'linewidth',2)
semilogy(f,median(psd_fosters,2),'linewidth',2)
grid on
legend('raw','SSS','Fosters','Location','best')
xlabel('frequency (Hz)')
ylabel('noise (fT/rtHz)')
title('PSD of Reconstructed Simulated Data')
set(gca,'FontSize',12)
set(gcf,'color','w')
