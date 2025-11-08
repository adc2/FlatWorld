function sources=sourcegrid(N)
%sources=sourcegrid(N) - source positions on a grid
%
%  sources: list of x,y pairs
%
%  N: size of grid side (nsources = N^2)


a=repmat(linspace(-1,1,N)',1,N); 
b=repmat(linspace(-1,1,N),N,1); 
sources=[a(:),b(:)];
