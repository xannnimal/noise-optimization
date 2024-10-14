% Flux through pick-up loops (magnetometers and gradiometers) based on the
% 1987 formula by J. Sarvas
% vuovec = dipole_field_sarvas(rs,q,r0,R,NX,NY,NZ,mags) where
% rs: center of the conducting sphere
% q: dipole moment vector
% r0: location of the dipole in the coordinate system of the conductor
% R: sensor locations
% NX,NY,NZ: direction vectors of the sensors (z-direction is the normal)
% mags: indices of the magnetometer channels (others are assumed
% gradiometers)

function vuovec = dipole_field_sarvas(rs,q,r0,R,NX,NY,NZ,mags)

mag_size = 21.0e-3; % Side length of a square magnetometer loop (in meters)
weights_mag = [16/81;25/324;25/324;25/324;25/324;10/81;10/81;10/81;10/81]; % Weights for the discretized surface integral
d = sqrt(3/5)*(mag_size/2);
D_mag = [0 0; d d; -d d; -d -d; d -d; 0 d; 0 -d; d 0; -d 0]';
baseline = 16.69e-3; % Gradiometer baseline (distance between the centers of the oppositely wound loops
weights_grad = (1/(4*baseline))*ones(4,1);
weights_grad(5:8) = -(1/(4*baseline))*ones(4,1);
dx1 = 5.89e-3;
dx2 = 10.8e-3;
dy = 6.71e-3;
D_grad = [dx1 dy; dx2 dy; dx1 -dy; dx2 -dy; -dx1 dy; -dx2 dy; -dx1 ...
   -dy; -dx2 -dy]';

nchan = size(R,2); % Number of channels
for i = 1:nchan
    if isempty(find(mags==i))
        vuovec(i) = grad_vaste(q,r0,R(:,i)-rs,NX(:,i),NY(:,i),NZ(:,i),weights_grad,D_grad);
   else
        vuovec(i) = magn_vaste(q,r0,R(:,i)-rs,NX(:,i),NY(:,i),NZ(:,i),weights_mag,D_mag);
   end
end


function vuo = magn_vaste(q,r0,rm,nx,ny,nz,weights_mag,D_mag)

for i = 1:9
   rint = rm + D_mag(1,i)*nx + D_mag(2,i)*ny;
   a = rint - r0;
   na = norm(a);
   nr = norm(rint);
   F = na*(nr*na + nr^2 - dot(r0,rint));
   NABLA_F = ((1/nr)*na^2 + (1/na)*dot(a,rint) + 2*na + 2*nr)*rint - ...
   (na + 2*nr + (1/na)*dot(a,rint))*r0;
   qr = cross(q,r0);
   B(:,i) = (F*qr - dot(qr,rint)*NABLA_F)/F^2;
end
vuo = dot(B*weights_mag,nz);


function vuo = grad_vaste(q,r0,rm,nx,ny,nz,weights_grad,D_grad)

baseline = 16.69e-3;
for i = 1:8
   rint = rm + D_grad(1,i)*nx + D_grad(2,i)*ny;
   a = rint - r0;
   na = norm(a);
   nr = norm(rint);
   F = na*(nr*na + nr^2 - dot(r0,rint));
   NABLA_F = ((1/nr)*na^2 + (1/na)*dot(a,rint) +  2*na + 2*nr)*rint - ...
   (na + 2*nr + (1/na)*dot(a,rint))*r0;
   qr = cross(q,r0);
   B(:,i) = (F*qr - dot(qr,rint)*NABLA_F)/F^2;
end
vuo = dot(B*weights_grad,nz);



