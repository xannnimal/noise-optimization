%
% [Sin,SNin] = Sin_vsh_vv(r_sphere,R,EX,EY,EZ,ch_types,Lin)
%
% Calculate the internal SSS basis Sin using
% vector spherical harmonics
%
function [Sin,SNin] = Sin_vsh_vv(r_sphere,R,EX,EY,EZ,ch_types,Lin)

MAG = 1;
GRAD = 0;
mu0 = 1.25664e-6; % Permeability of vacuum
%
% For numerical surface integration:
%
mag_size = 21e-3;
baseline = 16.69e-3;
d = sqrt(3/5)*mag_size/2;
dx1 = 5.89e-3;
dx2 = 10.8e-3;
dy = 6.71e-3;
Dmag = [0 0; d d; -d d; -d -d; d -d; 0 d; 0 -d; d 0; -d 0]';
Dgrad = [dx1 dy; dx2 dy; dx1 -dy; dx2 -dy; -dx1 dy; -dx2 dy; -dx1 -dy; -dx2 -dy]';
weights_mag = [16/81 25/324 25/324 25/324 25/324 10/81 10/81 10/81 10/81]';
for j = 1:8
   if j <= 4
      weights_grad(j) = 1/(4*baseline);
   else
      weights_grad(j) = -1/(4*baseline);
   end
end
weights_grad = weights_grad';

nchan = length(ch_types);
for ch = 1:nchan
  %disp(ch)
   count = 1;
   R(:,ch) = R(:,ch) - r_sphere;
   if ch_types(ch) == GRAD
      D = Dgrad;
      weights = weights_grad;
   elseif ch_types(ch) == MAG
      D = Dmag;
      weights = weights_mag;
   else
      error('Unknown sensor type!');
   end
   for l = 1:Lin
      for m = -l:l
	 Sin(ch,count) = -mu0*vsh_response(R(:,ch),EX(:,ch),EY(:,ch),EZ(:,ch),D,weights,l,m);
	 count = count + 1;
      end
   end
end
for i=(1:size(Sin,1))
    if ch_types(i)==1
        Sin(i,:)=Sin(i,:)*100;
    else
        Sin(i,:)=Sin(i,:);
    end
end
for j = 1:size(Sin,2)
  SNin(:,j) = Sin(:,j)/norm(Sin(:,j));
end


function Sin_element = vsh_response(r,ex,ey,ez,D,weights,l,m)

for j = 1:length(weights)
  r_this = r + D(1,j)*ex + D(2,j)*ey;
  rn = norm(r_this);
  theta = acos(r_this(3)/rn);
  phi = atan2(r_this(2),r_this(1));
  sint = sin(theta);
  sinp = sin(phi);
  cost = cos(theta);
  cosp = cos(phi);
  vs = vsh_modified_in(theta,phi,l,m)'/rn^(l+2);
  V(1,j) = vs(1)*sint*cosp + vs(2)*cost*cosp - vs(3)*sinp;
  V(2,j) = vs(1)*sint*sinp + vs(2)*cost*sinp + vs(3)*cosp;
  V(3,j) = vs(1)*cost - vs(2)*sint;
end
Sin_element = dot(V*weights,ez); % Cartesian coordinates
                                 %Sin_element = Sin_element/sqrt((l+1)*(2*l+1));  % Back to orthonormal presentation