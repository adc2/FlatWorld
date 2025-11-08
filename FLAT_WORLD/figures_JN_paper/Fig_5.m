clear
addpath ../FW

rng('default');
rand(1,37); % --> nice picture

figure(1); clf
set(gcf,'position',[300   597   340 340])

% color map
EXP=1; % exponent to tweak gradient
if 1
    ramp=linspace(1,0,128)';
    cmap=ones(size(jet)); 
    cmap(129:end,[2 3])=repmat(ramp.^EXP,1,2); cmap(128:-1:1,[1 2])=repmat(ramp.^EXP,1,2); 
else
    ramp=linspace(0,1,128)';
    cmap=zeros(size(jet)); 
    cmap(129:end,1)=repmat(ramp.^EXP,1,1); cmap(128:-1:1,3)=repmat(ramp.^EXP,1,1); 
end


nsensors=5; 
sensors=sensor_positions(nsensors);
sourcepositions=[0.5,-.4];
S=sensor_positions(nsensors);
M=1./sqdist(sourcepositions,S); % forward matrix
s=randn(1000,size(sourcepositions,1)); % sources time series
X=s*M; % observations
topcs=nt_pca0(X); %? PCA

NNN=400; % defines density of source position grid
source=source_grid(NNN);
G=source2sensor(source,S); 
G=reshape(G,NNN*NNN,nsensors);

F=topcs(:,2:end); % spatial filter
GPF=reshape(G*F,NNN,NNN,4); % gain 

for iNullFilter=1:4

    subplot(2,2,iNullFilter);
    G=GPF(:,:,iNullFilter);
    mask=zerox(G,1); % set pixels near zero to zero
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
    h=scatter(sensors(:,1)*a,sensors(:,2)*a, 300, 'k', 'filled');hold on
    h=scatter(sensors(:,1)*a,sensors(:,2)*a, 200, 'w', 'filled');
    % a=1.05;
    % plot(sensors(:,1)*a,sensors(:,2)*a, '.', 'color', [0,0.2,1],'markersize',40);
    xlim([-a,a]); ylim([-a,a]);
    
    set(gca,'clim',4*[-1 1]);
%    plot_tweak([0 0.06 0 -0.06], h7)
    % h{iNullFilter}=colorbar('southoutside');
    % set(h{iNullFilter},'ticks', [-2 -1 0 1 2], 'ticklabels',{'-100','-10', '0','10', '100'},'fontsize',14, 'limits',[-2.5 2.5], 'fontsize',14)
    % set(get(h{iNullFilter},'label'), 'string', 'gain');

end

