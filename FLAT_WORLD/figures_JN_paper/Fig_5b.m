clear
addpath ../FW

figure(1); clf
set(gcf,'position',[300   597   900   340])

EXP=1;
if 1
    ramp=linspace(1,0,128)';
    cmap=ones(size(jet)); 
    cmap(129:end,[2 3])=repmat(ramp.^EXP,1,2); cmap(128:-1:1,[1 2])=repmat(ramp.^EXP,1,2); 
else
    ramp=linspace(0,1,128)';
    cmap=zeros(size(jet)); 
    cmap(129:end,1)=repmat(ramp.^EXP,1,1); cmap(128:-1:1,3)=repmat(ramp.^EXP,1,1); 
end

% grid of positions at which to sample gain pattern
NNN=400; 
probes=source_grid(NNN);


 % 4 sources
sourcepositions=0.5*[-1,-1; 1,-1; 0,0;-1 1];
nsensors=5; % 5 sensors
S=sensor_positions(nsensors);
s=randn(1000,size(sourcepositions,1)); % sources time series
s=s*randn(size(s,2)); % introduces some correlation between sources (to show it doesn't matter)
M=1./sqdist(sourcepositions,S);
X=s*M; % observations
topcs=nt_pca0(X); %? PCA
G=source2sensor(probes,S); 
G=reshape(G,NNN*NNN,nsensors);
F=topcs(:,5); % spatial filter
GPF=reshape(G*F,NNN,NNN,1); % gain 

h1=axes('position', [0.65    0.080    0.33    0.84]);

G=GPF(:,:,1);
mask=zerox(G); % set pixels near zero to zero
[G,ticks,ticklabels]=symlog(G*100,[],0.1); % --> signed log 
G=G.*mask; % set zero set to zero
G=max(-3.9,G);
G=zero_outside(G,-4);
imagesc([-1 1], [-1, 1],G');
drawcross(sourcepositions,'k',0.1, [], 'linewidth',2);
%cmap=turbo;
cmap(1,:)=1; 
%cmap(129,:)=0;

set(gca, 'colormap', cmap)
set(gca,'box','off')
set(get(gca,'xaxis'), 'visible','off')
set(get(gca,'yaxis'), 'visible','off')

hold on
N=100;
x=sin(2*pi*(1:N+1)/N);
y=cos(2*pi*(1:N+1)/N);
plot(x,y, 'k', 'linewidth', 2);
a=1.05;
h=scatter(S(:,1)*a,S(:,2)*a, 400, 'k', 'filled');hold on
h=scatter(S(:,1)*a,S(:,2)*a, 300, 'w', 'filled');
xlim([-a,a]); ylim([-a,a]);

set(gca,'clim',[-4 4]);
%    plot_tweak([0 0.06 0 -0.06], h7)
% h{1}=colorbar('southoutside');
% set(h{1},'ticks', [-2 -1 0 1 2], 'ticklabels',{'-100','-10', '0','10', '100'},'fontsize',14, 'limits',[-2.5 2.5], 'fontsize',14)
% set(get(h{1},'label'), 'string', 'gain');

%drawcross(sourcepositions,'g');
%title('PC 5', 'rotation',90, 'position', [-1.2 0.2 0], 'fontsize', 16);


% 3 sources
sourcepositions=0.5*[0,0; 1,-1;-1 1]; 
nsensors=5; % 5 sensors
S=sensor_positions(nsensors);
s=randn(1000,size(sourcepositions,1)); % sources time series
s=s*randn(size(s,2)); % introduces some correlation between sources (to show it doesn't matter)
M=1./sqdist(sourcepositions,S);
X=s*M; % observations
topcs=nt_pca0(X); %? PCA
G=source2sensor(probes,S); 
G=reshape(G,NNN*NNN,nsensors);
F=topcs(:,4:5); % spatial filter
GPF=reshape(G*F,NNN,NNN,2); % gain 


% PC4
axes('position', [0.4    0.55    0.16    0.37])
G=GPF(:,:,1);
mask=zerox(G,2); % set pixels near zero to zero
[G,ticks,ticklabels]=symlog(G*100,[],0.1); % --> signed log 
G=G.*mask; % set zero set to zero
G=max(-3.9,G);
G=zero_outside(G,-4);
imagesc([-1 1], [-1, 1],G');
drawcross(sourcepositions,'k',0.2, [], 'linewidth',2);
%cmap=turbo;
cmap(1,:)=1; 
%cmap(129,:)=0;

set(gca, 'colormap', cmap)
set(gca,'box','off')
set(get(gca,'xaxis'), 'visible','off')
set(get(gca,'yaxis'), 'visible','off')

hold on
N=100;
x=sin(2*pi*(1:N+1)/N);
y=cos(2*pi*(1:N+1)/N);
plot(x,y, 'k', 'linewidth', 2);
a=1.05;
h=scatter(S(:,1)*a,S(:,2)*a, 300, 'k', 'filled');hold on
h=scatter(S(:,1)*a,S(:,2)*a, 200, 'w', 'filled');
xlim([-a,a]); ylim([-a,a]);

set(gca,'clim',[-4 4]);
%title('PC 4', 'rotation',90, 'position', [-1.2 0.2 0], 'fontsize', 16);

% PC5
axes('position', [0.4    0.11    0.16    0.37])
G=GPF(:,:,2);
mask=zerox(G,2); % set pixels near zero to zero
[G,ticks,ticklabels]=symlog(G*10,[],0.1); % --> signed log 
G=G.*mask; % set zero set to zero
G=max(-3.9,G);
G=zero_outside(G,-4);
imagesc([-1 1], [-1, 1],G');
drawcross(sourcepositions,'k',0.2, [], 'linewidth',2);
%cmap=turbo;
cmap(1,:)=1; 
%cmap(129,:)=0;

set(gca, 'colormap', cmap)
set(gca,'box','off')
set(get(gca,'xaxis'), 'visible','off')
set(get(gca,'yaxis'), 'visible','off')

hold on
N=100;
x=sin(2*pi*(1:N+1)/N);
y=cos(2*pi*(1:N+1)/N);
plot(x,y, 'k', 'linewidth', 2);
a=1.05;
h=scatter(S(:,1)*a,S(:,2)*a, 300, 'k', 'filled');hold on
h=scatter(S(:,1)*a,S(:,2)*a, 200, 'w', 'filled');
xlim([-a,a]); ylim([-a,a]);

set(gca,'clim',[-4 4]);
%    plot_tweak([0 0.06 0 -0.06], h7)
% h{1}=colorbar('southoutside');
% set(h{1},'ticks', [-2 -1 0 1 2], 'ticklabels',{'-100','-10', '0','10', '100'},'fontsize',14, 'limits',[-2.5 2.5], 'fontsize',14)
% set(get(h{1},'label'), 'string', 'gain');

%drawcross(sourcepositions,'g');
%title('PC 5', 'rotation',90, 'position', [-1.2 0.2 0], 'fontsize', 16);


% 2 sources
sourcepositions=0.5*[1,-1;-1 1]; 
nsensors=5; % 5 sensors
S=sensor_positions(nsensors);
s=randn(1000,size(sourcepositions,1)); % sources time series
s=s*randn(size(s,2)); % introduces some correlation between sources (to show it doesn't matter)
M=1./sqdist(sourcepositions,S);
X=s*M; % observations
topcs=nt_pca0(X); %? PCA
G=source2sensor(probes,S); 
G=reshape(G,NNN*NNN,nsensors);
F=topcs(:,3:5); % spatial filter
GPF=reshape(G*F,NNN,NNN,3); % gain 


% PC3
h3=axes('position', [0.18   0.67    0.1   0.25])
G=GPF(:,:,1);
mask=zerox(G,3); % set pixels near zero to zero
[G,ticks,ticklabels]=symlog(G*100,[],0.1); % --> signed log 
G=G.*mask; % set zero set to zero
G=max(-3.9,G);
G=zero_outside(G,-4);
imagesc([-1 1], [-1, 1],G');
drawcross(sourcepositions,'k',0.2, [], 'linewidth',2);
%cmap=turbo;
cmap(1,:)=1; 
%cmap(129,:)=0;

set(gca, 'colormap', cmap)
set(gca,'box','off')
set(get(gca,'xaxis'), 'visible','off')
set(get(gca,'yaxis'), 'visible','off')

hold on
N=100;
x=sin(2*pi*(1:N+1)/N);
y=cos(2*pi*(1:N+1)/N);
plot(x,y, 'k', 'linewidth', 2);
a=1.05;
h=scatter(S(:,1)*a,S(:,2)*a, 200, 'k', 'filled');hold on
h=scatter(S(:,1)*a,S(:,2)*a, 100, 'w', 'filled');
xlim([-a,a]); ylim([-a,a]);

set(gca,'clim',[-4 4]);
%title('PC 4', 'rotation',90, 'position', [-1.2 0.2 0], 'fontsize', 16);

% PC4
h5=axes('position', [0.18    0.39   0.1    0.25])
G=GPF(:,:,2);
mask=zerox(G,3); % set pixels near zero to zero
[G,ticks,ticklabels]=symlog(G*100,[],0.1); % --> signed log 
G=G.*mask; % set zero set to zero
G=max(-3.9,G);
G=zero_outside(G,-4);
imagesc([-1 1], [-1, 1],G');
drawcross(sourcepositions,'k',0.2, [], 'linewidth',2);
%cmap=turbo;
cmap(1,:)=1; 
%cmap(129,:)=0;

set(gca, 'colormap', cmap)
set(gca,'box','off')
set(get(gca,'xaxis'), 'visible','off')
set(get(gca,'yaxis'), 'visible','off')

hold on
N=100;
x=sin(2*pi*(1:N+1)/N);
y=cos(2*pi*(1:N+1)/N);
plot(x,y, 'k', 'linewidth', 2);
a=1.05;
h=scatter(S(:,1)*a,S(:,2)*a, 200, 'k', 'filled');hold on
h=scatter(S(:,1)*a,S(:,2)*a, 100, 'w', 'filled');
xlim([-a,a]); ylim([-a,a]);

set(gca,'clim',[-4 4]);
%title('PC 5', 'rotation',90, 'position', [-1.2 0.2 0], 'fontsize', 16);


% PC5
h5=axes('position', [0.18   0.11    0.1   0.25])
G=GPF(:,:,3);
mask=zerox(G,3); % set pixels near zero to zero
[G,ticks,ticklabels]=symlog(G*10,[],0.1); % --> signed log 
G=G.*mask; % set zero set to zero
G=max(-3.9,G);
G=zero_outside(G,-4);
imagesc([-1 1], [-1, 1],G');
drawcross(sourcepositions,'k',0.2, [], 'linewidth',2);
%cmap=turbo;
cmap(1,:)=1; 
%cmap(129,:)=0;

set(gca, 'colormap', cmap)
set(gca,'box','off')
set(get(gca,'xaxis'), 'visible','off')
set(get(gca,'yaxis'), 'visible','off')

hold on
N=100;
x=sin(2*pi*(1:N+1)/N);
y=cos(2*pi*(1:N+1)/N);
plot(x,y, 'k', 'linewidth', 2);
a=1.05;
h=scatter(S(:,1)*a,S(:,2)*a, 200, 'k', 'filled');hold on
h=scatter(S(:,1)*a,S(:,2)*a, 100, 'w', 'filled');
xlim([-a,a]); ylim([-a,a]);

set(gca,'clim',[-4 4]);
%    plot_tweak([0 0.06 0 -0.06], h7)
% h{1}=colorbar('southoutside');
% set(h{1},'ticks', [-2 -1 0 1 2], 'ticklabels',{'-100','-10', '0','10', '100'},'fontsize',14, 'limits',[-2.5 2.5], 'fontsize',14)
% set(get(h{1},'label'), 'string', 'gain');

%drawcross(sourcepositions,'g');
%title('PC 5', 'rotation',90, 'position', [-1.2 0.2 0], 'fontsize', 16);




