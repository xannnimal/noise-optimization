function X = xi(S,PHI,nm,nt,ni)
%% from Samu Taulu
%   "S" is the full normalized SSS basis, in and out 
%   "PHI" is the data (a vector or a matrix)
%   "nm" is the same as L_in
%   "nt" equals Lout - 1
%   "ni" is the number of trials or iterations
min_change = 1e-2;
ni_default = 10;

if nargin < 5
   ni = ni_default;
end
nsamp = size(PHI,2);
dim_m = (nm+1)^2-1;
for n = 1:nm
   dim1 = (n-1+1)^2;
   dim2 = (n+1)^2-1;
   dimv(n,:) = [dim1 dim2];
end
for n = 1:nt+1
   dit1 = dim_m + n^2;
   dit2 = dim_m + (n+1)^2-1;
   ditv(n,:) = [dit1 dit2];
end
X = zeros(size(S,2),nsamp);
count = 1;
for n = 1:nm
   if n-1 <= nt
      indices{n} = [dimv(n,1):dimv(n,2) ditv(n,1):ditv(n,2)];
   else
      indices{n} = [dimv(n,1):dimv(n,2)];
   end
   pS{n} = pinv(S(:,indices{n}));
end
condition_Lin=cond(S(:,nm))
condition_tot=cond(S)
if nm >= nt
   while count <= ni 
      for i = 1:nm
	 n = i;
	 inds = indices{n}; 
	 X(inds,:) = zeros(length(inds),nsamp);
	 XN = pS{n}*(PHI-S*X);
	 X(inds,:) = XN;
      end
      count = count + 1;
   end
else
   error('nm should be greater than nt!');
end

