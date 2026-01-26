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
return
%add gaussian noise at 10 percent of max value of phi_0
rng("default") %set seed
noise = 0.2*randn(size(phi_0,1),2*size(phi_0,2))*scale;
% add the dipole times after the baseline of only noise
phi_0=[zeros(size(phi_0,1),size(phi_0,2)),phi_0];
%alpha=[zeros(size(alpha,1),size(alpha,2)),alpha];
phi_0=phi_0+noise;
%N = cov(noise');
N = diag(cov(noise'));

for i=(1:size(phi_0,1))
    if mod(i,3)==0 %every third is a magnetometer
        phi_0(i,:)=phi_0(i,:)*100;
    else
        phi_0(i,:)=phi_0(i,:);
    end
end


%% calculate SSS basis
[Sin,SNin] = Sin_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,Lin);
[Sout,SNout] = Sout_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,Lout);
S = [SNin]; %change to [SNin SNout] for full basis, or SNin_tot for mSSS
pS=pinv(S); %95x306
XN=pS*phi_0; %multiple moments 95x500
data_rec=real(SNin*XN(1:size(SNin,2),:));
%diff = XN-alpha;

%% iterative recon
% ni=10;
% XN_it = xi([SNin,SNout],phi_0,Lin,2,ni);
% data_rec_it = real(SNin*XN_it(1:size(SNin,2),:));

%% try different N matrices
% load in covariance from MNE-Python
% covariance = load("covariance.mat","cov");
% N_check = covariance.cov;

%% reconstruct data fosters
%%using calculated multipole moments
alpha_cov = cov(alpha');
for i=(1:size(Sin,2))
    for j=(1:size(Sin,2))
        alpha_cov_new(i,j)=alpha_cov(i,j)*norm(Sin(:,i))*norm(Sin(:,j));
    end
end
m_alpha = mean(alpha,2); 

%%using estimated multipole moments
% clear S; clear pS; clear XN; clear alpha_cov_new;
% S = SNin;
% pS=pinv(S);
% XN=pS*phi_0;
% alpha_cov_new = cov(XN'); 
% alpha = XN';
% m_alpha = mean(alpha,1)'; 

%%same for both multipole moments
S_star = conj(S)'; 
first = pinv(S*alpha_cov_new*S_star+N);
B = alpha_cov_new*S_star*first;
b = m_alpha - B*S*m_alpha;
%better estimate for multipole moments
for i=(1:size(timestep:timestep:timestep*size(phi_0,2),2))
    x_bar(:,i) = B*phi_0(:,i) +b;
end
data_rec_fosters= real(SNin*x_bar(1:size(SNin,2),:));

return
%% try with mSSS basis
% thresh = 0.005;
% center1= [-0.00350699, 0.01138051, 0.05947857] - [0,0,0.05]; 
% center2= [-0.00433911, 0.04081329, 0.05194245] - [0,0,0.05];
% [SNin_tot, SNout] = multi_sss(center1,center2,R,EX,EY,EZ,ch_types,Lin, Lout, thresh);
% clear S; clear pS; clear XN; clear alpha_cov_new;
% S = [SNin_tot]; %change to [SNin SNout] for full basis, or SNin_tot for mSSS
% SNin = SNin_tot;
% pS=pinv(S);
% XN=pS*phi_0;
% alpha_cov_new = cov(XN'); 
% alpha = XN';
% m_alpha = mean(alpha,1)'; 
% S_star = conj(S)'; 
% first = pinv(S*alpha_cov_new*S_star+N);
% B = alpha_cov_new*S_star*first;
% b = m_alpha - B*S*m_alpha;
% %better estimate for multipole moments
% clear x_bar
% for i=(1:size(timestep:timestep:timestep*size(phi_0,2),2))
%     x_bar(:,i) = B*phi_0(:,i) +b;
% end
% data_rec_fos_mSSS= real(SNin*x_bar(1:size(SNin,2),:));


%% sub angles
% will be much lower if angle is calculated for just times with signal, not
% baseline
for i=(501:1000) % only do signal times, not baseline noise
    check_data(i) = subspace(data_rec(:,i), phi_0)*180/pi;
    check_data_it(i) = subspace(phi_0(:,i), data_rec_it)*180/pi;
    check_data_fosters(i) = subspace(data_rec_fosters(:,i), phi_0)*180/pi;
    %check_data_fos_mSSS(i) = subspace(data_rec_fos_mSSS(:,i), phi_0)*180/pi;
end
check_data_min = min(check_data);
check_data_max = max(check_data);
check_data_mean = mean(check_data);

check_data_it_min = min(check_data_it);
check_data_it_max = max(check_data_it);
check_data_it_mean = mean(check_data_it);

check_data_fosters_min = min(check_data_fosters);
check_data_fosters_max = max(check_data_fosters);
check_data_fosters_mean = mean(check_data_fosters);

% check_data_fos_mSSS_min = min(check_data_fos_mSSS);
% check_data_fos_mSSS_max = max(check_data_fos_mSSS);
% check_data_fos_mSSS_mean = mean(check_data_fos_mSSS);


%% plot raw and recon data
chan=3  ;
figure(1)
hold on
times2 = timestep:timestep:timestep*size(phi_0,2);
plot(times2,phi_0(chan,:),'g') %raw
plot(times2,data_rec(chan,:),'r') %sss
plot(times2,data_rec_it(chan,:),'c')
plot(times2,data_rec_fosters(chan,:),'b') %foster
%plot(times2,data_rec_fos_mSSS(chan,:),'m') %foster
title('306 SQUID, Currrent Dipole [5cm,0,0] Reconstruction')
xlabel('Time (sec)')
ylabel('Dipole Signal, Chan 3 (T)')
legend({'raw','SSS','iter','Fosters'},'location','northwest')
%legend({'raw','SSS','Fosters'},'location','northwest')
%legend({'raw','Fosters'},'location','northwest')


%% SNR calculations (signal/noise)
% noise from baseline, peak from data
% if A is a matrix whose columns are random variables and whose rows are observations 
% then S is a row vector containing the standard deviation corresponding to each column
n1=40; n2=200;
s1=640; s2=800;
%peak/noise
% SNR_raw = mean(std(phi_0(:,s1:s2)))/mean(std(phi_0(:,n1:n2)));
% SNR_vsh = mean(std(data_rec(:,s1:s2)))/mean(std(data_rec(:,n1:n2)));
% SNR_it = mean(std(data_rec_it(:,s1:s2)))/mean(std(data_rec_it(:,n1:n2)));
% SNR_fost = mean(std(data_rec_fosters(:,s1:s2)))/mean(std(data_rec_fosters(:,n1:n2)));
%SNR_fost_mSSS = mean(std(data_rec_fos_mSSS(:,s1:s2)))/mean(std(data_rec_fos_mSSS(:,n1:n2)));

SNR_raw = mean((std(phi_0(:,s1:s2))).^2)/mean((std(phi_0(:,n1:n2))).^2);
SNR_vsh = mean((std(data_rec(:,s1:s2))).^2)/mean((std(data_rec(:,n1:n2))).^2);
SNR_it = mean((std(data_rec_it(:,s1:s2))).^2)/mean((std(data_rec_it(:,n1:n2))).^2);
SNR_fost = mean((std(data_rec_fosters(:,s1:s2))).^2)/mean((std(data_rec_fosters(:,n1:n2))).^2);

%% calculate PSD 
%estimate of PSD. "periofogram" function has lots of variations we should
% investigate, this is the simplest implementation for now
[psd_raw,f]=periodogram(phi_0);
[psd_sss,f]=periodogram(data_rec);
[psd_it,f]=periodogram(data_rec_it);
[psd_fosters,f]=periodogram(data_rec_fosters);
%[psd_fos_mSSS,f]=periodogram(data_rec_fos_mSSS);
%plot median PSD
figure(2)
hold on
semilogy(f,median(psd_raw,2),'g','linewidth',2)
semilogy(f,median(psd_sss,2),'r','linewidth',2)
semilogy(f,median(psd_it,2),'c','linewidth',2)
semilogy(f,median(psd_fosters,2),'b','linewidth',2)
%semilogy(f,median(psd_fos_mSSS,2),'m','linewidth',2)
grid on
legend('raw','SSS','iter','Fosters','Location','best')
xlabel('frequency (Hz)')
ylabel('noise (fT/rtHz)')
title('PSD of Reconstructed Simulated Data')
set(gca,'FontSize',12)
set(gcf,'color','w')
