function S=sensor_positions(nsensors, delete)
if nargin<2; delete=1; else delete=0; end
% create sensor positions

nsensors=nsensors+1;
S=[sin(2*pi*(1:nsensors)/nsensors)', cos(2*pi*(1:nsensors)/nsensors)'];
[~,idx]=max(S(:,2));
if delete
    S(idx,:)=[]; % delete lowest positions
end
