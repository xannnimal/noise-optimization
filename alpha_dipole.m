% alpha = alpha_dipole(rq,q,Lin)
%
% Vector of multipole moments up to order l = Lin
% for a current dipole q at location rq
% 
function alpha = alpha_dipole(rq,q,Lin)

rqn = norm(rq);
thetaq = acos(rq(3)/rqn);
phiq = atan2(rq(2),rq(1));
qs(1) = q(1)*cos(thetaq)*cos(phiq) + q(2)*cos(thetaq)*sin(phiq) - q(3)*sin(thetaq);
qs(2) = -q(1)*sin(phiq) + q(2)*cos(phiq);

count = 1;
for l = 1:Lin
   scale = (1/(2*l+1))*sqrt(l/(l+1))*rqn^l;
   for m = -l:l
      Xrq = vspharm(thetaq,phiq,l,m);
      alpha(count) = scale*dot(i*conj(Xrq),qs);
      count = count + 1;
   end
end      
alpha = alpha';
