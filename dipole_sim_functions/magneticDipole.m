function magneticFlux = magneticDipole(chanpos,EX,EY,EZ,dip_pos, dip_mom, ch_types)
% Iman and Xan 2024
% Calculating the magnetic field from a dipole using Griffiths EM eq. 5.89
% dip_mom: magnetic moment 1x3
% dip_pos: location of dipole
% chanpos: channel position matrix nchanx3
% sensing_dir: the orientation/measuring direction of the sensors
% grad_pos_ex: for gradiometers, need EX positions, for OPM sensors this
%                part of the function is never used
ch_types=ch_types';
mag_size= 21.0e-3;
weights_mag = [16/81;25/324;25/324;25/324;25/324;10/81;10/81;10/81;10/81];
d = sqrt(3/5)*(mag_size/2);
D_mag = [0 0; d d; -d d; -d -d; d -d; 0 d; 0 -d; d 0; -d 0]';
baseline = 16.69e-3;
weights_grad = (1/(4*baseline))*ones(4,1);
weights_grad(5:8) = -(1/(4*baseline))*ones(4,1);
dx1 = 5.89e-3;
dx2 = 10.8e-3;
dy = 6.71e-3;
D_grad = [dx1 dy; dx2 dy; dx1 -dy; dx2 -dy; -dx1 dy; -dx2 dy; -dx1 ...
   -dy; -dx2 -dy]';

%make r vector from (x,y,z) chanpos 
% for i=(1:size(chanpos,2))
%     %r(i)=sqrt(chanpos(i,1)^2+chanpos(i,2)^2+chanpos(i,3)^2);
%     EZ(i,:)=EZ(i,:)/norm(EZ(i,:));
% end
%calculate B for each channel position
nchan=size(chanpos,2);
for k=(1:nchan)
    if ch_types(k)==1 %magnetometer calculations
        %expressing the location vector r as the product of its magnitude times the unit vector in its direction
        %rhat points from dipole location to sensor
        magneticFlux(k) = mag(chanpos(:,k),EX(:,k),EY(:,k),EZ(:,k),dip_pos, dip_mom,weights_mag, D_mag);
    else %gradiometers, ch_types=0
        magneticFlux(k) = grad(chanpos(:,k),EX(:,k),EY(:,k),EZ(:,k),dip_pos, dip_mom,weights_grad, D_grad);
    end
end

%magnetometer function
function mag_flux = mag(chanpos_i,EX_i,EY_i,EZ_i,dip_pos, dip_mom,weights_mag, D_mag)
    scale = 1.0e-7;
    for i=1:9
        r_integral = chanpos_i + D_mag(1,i)*EX_i + D_mag(2,i)*EY_i;
        rvec = r_integral - dip_pos;
        r_hat = norm(rvec);
        Bfield = (3*dot(dip_mom,rvec)/r_hat^5)*rvec;
        Bfield2= dip_mom/r_hat^3;
        B(:,i)= scale*(Bfield-Bfield2);
    end
    mag_flux = dot(B*weights_mag,EZ_i);
end

%gradiometer function
function mag_flux = grad(chanpos_i,EX_i,EY_i,EZ_i,dip_pos, dip_mom,weights_grad, D_grad)
    scale = 1.0e-7;
    for i=1:8
        r_integral = chanpos_i + D_grad(1,i)*EX_i + D_grad(2,i)*EY_i;
        rvec = r_integral - dip_pos;
        r_hat = norm(rvec);
        Bfield = (3*dot(dip_mom,rvec)/r_hat^5)*rvec;
        Bfield2= dip_mom/r_hat^3;
        B(:,i)= scale*(Bfield-Bfield2);
    end
    mag_flux = dot(B*weights_grad,EZ_i);
end
end