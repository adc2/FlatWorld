function G=zero_outside(G,value)% set to zerooutside circle
% 
% value: use this instead of zero

if nargin<2; value=0; end

NNN=size(G,1);

for ix=1:NNN
    for iy=1:NNN
        if (ix/NNN*2-1)^2+(iy/NNN*2-1)^2>1; G(ix,iy,:)=value; end
    end
end
