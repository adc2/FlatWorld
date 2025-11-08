clear

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



figure(1); clf
set(gcf,'position',[300   597   1200 280])

%rng(0);

% grid of positions at which to sample gain pattern
NNN=400; 
probes=source_grid(NNN);

% create several configurations
S{1}=sensor_positions(15);
T{1}=0.4*[1,1];
C{1}=0.7*(2*rand(14,2)-1);
S{3}=sensor_positions(5);
T{3}=0.4*[0,0];
T{3}=0.4*[1,1];
C{3}=0.5*[-1,-1; 1,-1; 0 1;-1 1; -1 0; 1 0; 1 1];
C{3}=0.5*[-1,-1; 1,-1; 0,0;-1 1; -1 0];
%C{2}=0.5*source_grid(3); C{2}(C{2}(:,1)>0,:)=[];
%S{3}=sensor_positions(9); S{3}( S{3}(:,1)>0,:)=[];
S{4}=sensor_positions(5);
T{4}=0.4*[0,0];
T{4}=0.4*[1,1];
%C{3}=0.5*source_grid(3);
C{4}=0.5*[-1,-1; 1,-1; ;-1 1; -1 0; 1 1; 0 1; 0 -1];
C{4}=0.5*[-1,-1; 1,-1; 0,0;-1 1; 1 0];
S{2}=sensor_positions(10); S{2}( S{2}(:,1)>0,:)=[];
T{2}=0.4*[-1 0];
C{2}=0.5*source_grid(7); C{2}( C{2}(:,1)<0,:)=[];
S{5}=sensor_positions(10); S{5}( S{5}(:,1)<0,:)=[];
T{5}=0.4*[-1 0];
C{5}=0.5*source_grid(7); C{5}( C{5}(:,1)<0,:)=[];

for iConf=1:5
    % create target and competitor sensor signals (separate to ease SNR
    % estimate)
    
    % 1 target
    Mt=1./sqdist(T{iConf},S{iConf});
    % multiple competitors
    Mc=1./sqdist(C{iConf},S{iConf});
    % source signals
    nsamples=10000;
    st=sin(2*pi*(1:nsamples)*10/nsamples)'; % 10 cycles
    sc=randn(nsamples,size(C{iConf},1));
    % sensor signals
    Xt=st*Mt;
    Xc=sc*Mc;
    % adjust amplitude of targe for SNR=1 on best sensor
    pwrt=mean(Xt.^2);
    pwrc=mean(Xc.^2);
    r=sqrt(max(pwrt./pwrc));
    Xt=Xt/r;
    % target SNR
    SNR=sqrt(1);
    Xt=Xt*SNR;
    % sensor signals
    X=Xc+Xt;
    
    % epoch, apply DSS
    X=nt_fold(X,nsamples/10);

    todss=nt_dss1(X);
    f=todss(:,1); % virtual electrode
    % apply to source and competitor separately
    zt=nt_mmat(Xt,f);
    zc=nt_mmat(Xc,f);

    % z=nt_mmat(X,f);
    % clf; subplot 121; nt_bsplot(X(:,1,:)); subplot 122; nt_bsplot(z); pause
    
    boost= (mean(zt.^2)/mean(zc.^2))/SNR.^2;
    disp(boost)

    subplot (1,5,iConf) % DSS
    
    % display gain pattern
    forward=source2sensor(probes,S{iConf}); 
    forward=reshape(forward,NNN*NNN,size(S{iConf},1));
    F=f;
    GF=reshape(forward*F,NNN,NNN,1); % grid to filteroutput
    gain_pattern=GF(:,:,1);

    if iConf==1; gain_pattern=-gain_pattern; end % aesthetics
    
    mask=zerox(gain_pattern); % set pixels near zero to zero
    [gain_pattern,ticks,ticklabels]=symlog(gain_pattern*100,[],0.1); % --> signed log 
    gain_pattern=gain_pattern.*mask; % set zero set to zero
    gain_pattern=max(-3.9,gain_pattern);
    gain_pattern=zero_outside(gain_pattern,-4);
    imagesc([-1 1], [-1, 1],gain_pattern');
    viscircles(T{iConf},0.07, 'color', 'k','linestyle','-', 'linewidth', 3, 'enhancevisibility',0); 
    drawcross(C{iConf},'k',0.07);
    
    %cmap=turbo;
    cmap(1,:)=1; 
    %cmap(129,:)=0;
    set(gca, 'colormap', cmap)
    set(gca,'box','off')
    set(get(gca,'xaxis'), 'visible','off')
    set(get(gca,'yaxis'), 'visible','off')
    
    hold on
    % circle around disk
    N=100;
    x=sin(2*pi*(1:N+1)/N);
    y=cos(2*pi*(1:N+1)/N);
    plot(x,y, 'k', 'linewidth', 2);
    a=1.05;
    h=scatter(S{iConf}(:,1)*a,S{iConf}(:,2)*a, 200, 'k', 'filled');hold on
    h=scatter(S{iConf}(:,1)*a,S{iConf}(:,2)*a, 100, 'w', 'filled');
    xlim([-a,a]); ylim([-a,a]);
    set(gca,'clim',[-4 4]);
    %    plot_tweak([0 0.06 0 -0.06], h7)
    h=colorbar('southoutside');
    set(h,'ticks', [-2 -1 0 1 2], 'ticklabels',{'-100','-10', '0','10', '100'},'fontsize',14, 'limits',[-2.5 2.5], 'fontsize',14)
    set(get(h,'label'), 'string', 'gain');

    if iConf==1
        hh=title('SNR ratio = \infty', 'interpreter','tex');
    else
        hh=title(['SNR ratio = ', num2str(boost, '%.1f')], 'interpreter','tex');
    end
    set(hh,'units','points', 'position',[74.2510 150 0])
    plot_tweak([0 0.05 0 -.1])
end

