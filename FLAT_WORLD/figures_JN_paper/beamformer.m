function B=beamformer(x,F)
% B=beamformer(x,F)
%   B: spatial filter
%   x: data, time X channels (or C^-1, inverse of covariance)
%   F: filter, channels X 1

if size(x,1)==size(x,2) && all (all(x'==x)); % covariance matrix
    Cinv=x;
else
    Cinv=nt_cov(x)^-1;
end
B=(F'*...
    Cinv...
    *F)^-1 ...
    * F' ...
    *Cinv;
B=B';
end % function beamformer