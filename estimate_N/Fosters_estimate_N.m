%% december 2025
% use foster's inverse to estimate N starting with known phantom dipole
% locations and S matrix. compare estimated x-bar with calculated x "alpha"
% as an optimization problem.

%% constant variables 
Lin = 8; % Truncation order of the internal VSH basis
Lout = 3; % Truncation order of the external VSH basis
dim_in = (Lin+1)^2 - 1; % Dimension of the internal SSS basis, should be 80

%% generate SQUID magnetometers
% load in positions from file, can be found on MNE-Python
coordsys = 'device'; 
[R,EX,EY,EZ] = fiff_getpos('phantom_1000nam_default_IASoff_evoked_1.fif', coordsys);
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

%% load evoked data for 32 dipoles
numdipole = 32; % try with only done dipole
phantom_data = cell(1, numdipole);
for k = 1:numdipole
  i = k-1;
  filename = sprintf('phantom_1000nam_default_IASoff_evoked_%d.fif', i);
  [data] = fiff_read_evoked(filename);
  phantom_data{k}=data.evoked.epochs;
end
phantom_times=data.evoked.times;
save('phantom_times.mat','phantom_times')
save('phantom_data.mat','phantom_data')
clear data

figure(9)
hold on
title('phantom data')
plot(phantom_times,phantom_data{1,1}(:,:))
hold off

%% load phantom dipole positions
phantom_dip = load("Phantom_dipoles.mat");
pos_actual = (phantom_dip.pos)';
ori_actual = (phantom_dip.ori)';

% figure(1)
% hold on
% scatter3(pos_actual(1,:),pos_actual(2,:),pos_actual(3,:),'r')
% scatter3(R(1,:),R(2,:),R(3,:),'g')

%% setup constant knowns
% create simulated phantom dipole activations
rs=[0,0,0];
f_start = 20; % start frequency
f_end = 20; % end frequency
times = phantom_times;
q_t= cell(1, numdipole);
for k=(1:numdipole)
    for i=(1:3)
        q_t{k}(i,:) = ori_actual(i,k)*sin(-2*pi*(f_start*times));
    end
    q_t{k}(:,1:51)=0; %add period of no activations before and after
    q_t{k}(:,151:201)=0;
end
% calculate known multipole moments from phantom dipoles
scale = 1e-12;
X_c = cell(1, numdipole);
phi_0 = cell(1, numdipole);
for k=(1:numdipole)
    for i=(1:size(times,2))
        X_c{k}(:,i) = alpha_dipole(pos_actual(:,k)',q_t{1,k}(:,i),Lin)';
        phi_0{k}(:,i) = dipole_field_sarvas(rs',q_t{1,k}(:,i),pos_actual(:,k),R,EX,EY,EZ,mags)*scale';
    end
end

% figure(10)
% hold on
% title('simulated dipole activations')
% plot(times,real(X_c{1,1}(:,:)))
% hold off

% calculate SSS basis
[Sin,SNin] = Sin_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,Lin);
[Sout,SNout] = Sout_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,Lout);
S = [SNin]; %change to [SNin SNout] for full basis, or SNin_tot for mSSS

%% reconstruct data fosters
%%using calculated multipole moments
%first guess for N, calculated MNE-Python empirical covariance
covar = load('N_1000nam_default_IASoff.mat');
N0= diag(covar.covar);
%resid = fosters_residual(N,X_c,S,phantom_data,phantom_times,numdipole);

%initial point
scale = 1e10;
objFun = @(N) fosters_residual(N,X_c,S,phantom_data,phantom_times,numdipole); %objective function
lb = 1e-30*ones([size(N0,1),size(N0,2)]); %lower bound
ub = 1e-15*ones([size(N0,1),size(N0,2)]); %upper bound
%x0 = N0;


%run for 1000 iterations, plot condition number
options = optimoptions('simulannealbnd','Display','iter','MaxIterations',100, ...
    'PlotFcn',{@saplotbestx,@saplotbestf,@saplotx,@saplotf}); %'OutputFcn',@save_output_iter_cond)
tic
[N,naf,exitflag,output] = simulannealbnd(objFun,N0,lb,ub,options);
toc

%N=N/scale;
%save('N_1050_opt_1000nam_default_IASoff_evoked.mat','N');


%%%% do fosters with optimized N'
%%using estimated multipole moments
clear S; clear pS; clear XN; clear alpha_cov_new;
S = SNin;
pS=pinv(S);
phi_0=phantom_data{1,1}(:,:);
XN=pS*phi_0;
alpha_cov_new = cov(XN'); 
alpha = XN';
m_alpha = mean(alpha,1)'; 

%%same for both multipole moments
S_star = conj(S)'; 
first = pinv(S*alpha_cov_new*S_star+N);
B = alpha_cov_new*S_star*first;
b = m_alpha - B*S*m_alpha;
%better estimate for multipole moments
for i=(1:size(phantom_times,2))
    x_bar(:,i) = B*phi_0(:,i) +b;
end
data_rec_fosters= real(SNin*x_bar(1:size(SNin,2),:));

figure(11)
hold on
title('N-optimized Foster reconstructed data')
plot(phantom_times,data_rec_fosters(:,:))
hold off