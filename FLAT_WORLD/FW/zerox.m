function z=zerox2(x, widen)
%z=zerox2(x, widen) - 2D zero crossings
%
%  z: zero crossing matrix (0 if zero crossing, 1 otherwise)
%
%  x: data matrix (2D)
%  widen: number of widening iterations [default: 0]

if nargin<2; widen=0; end

if ndims(x)==3
    z=ones(size(x));
    for iTrial=1:size(x,3);
        z(:,:,iTrial)=zerox(x(:,:,iTrial));
    end
    return
end

x=x(:,:);
mask=ones(size(x));
mask(1:end-1,:) = mask(1:end-1,:) .* (x(1:end-1,:) .* x(2:end,:))>0;
mask(:,1:end-1) = mask(:,1:end-1) .* (x(:,1:end-1) .* x(:,2:end))>0;
mask(2:end,:) = mask(2:end,:) .* (x(1:end-1,:) .* x(2:end,:))>0;
mask(:,2:end) = mask(:,2:end) .* (x(:,1:end-1) .* x(:,2:end))>0;
z=mask;

% widen
for iIter=1:widen
    z(1:end-1,:) = z(1:end-1,:) .* z(2:end,:);
    z(:,1:end-1) = z(:,1:end-1) .* z(:,2:end);
    z(2:end,:) = z(2:end,:) .* z(2:end,:);
    z(:,2:end) = z(:,2:end) .* z(:,2:end);
end
