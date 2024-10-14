%
% [Sin,SNin] = Sin_vsh_vv(r_sphere,R,EX,EY,EZ,ch_types,Lout)
%
% Calculate the external SSS basis Sin using
% vector spherical harmonics
%
function [Sout,SNout] = Sout_vsh_vv(r_sphere,R,EX,EY,EZ,ch_types,Lout)

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
   for l = 1:Lout
      for m = -l:l
	 Sout(ch,count) = -mu0*vsh_response(R(:,ch),EX(:,ch),EY(:,ch),EZ(:,ch),D,weights,l,m);
	 count = count + 1;
      end
   end
end
for i=(1:size(Sout,1))
    if ch_types(i)==1 %every third is a magnetometer
        Sout(i,:)=Sout(i,:)*100;
    else
        Sout(i,:)=Sout(i,:);
    end
end
for j = 1:size(Sout,2)
  SNout(:,j) = Sout(:,j)/norm(Sout(:,j));
end


function Sout_element = vsh_response(r,ex,ey,ez,D,weights,l,m)
for j = 1:length(weights)
  r_this = r + D(1,j)*ex + D(2,j)*ey;
  rn = norm(r_this);
  theta = acos(r_this(3)/rn);
  phi = atan2(r_this(2),r_this(1));
  sint = sin(theta);
  sinp = sin(phi);
  cost = cos(theta);
  cosp = cos(phi);
  ws = vsh_modified_out(theta,phi,l,m)'*rn^(l-1);
  W(1,j) = ws(1)*sint*cosp + ws(2)*cost*cosp - ws(3)*sinp;
  W(2,j) = ws(1)*sint*sinp + ws(2)*cost*sinp + ws(3)*cosp;
  W(3,j) = ws(1)*cost - ws(2)*sint;
end
Sout_element = dot(W*weights,ez);