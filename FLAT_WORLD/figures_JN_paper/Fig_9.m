clear
addpath ../FW

if 2~=exist('nt_version');
    error('download NoiseTools, adjust path');
end


figure(1); clf
set(gcf,'position',[300   597   1200 600])

%rng(0);

load tmp/VE_globs
cmap=VE_globs.cmap;

% grid of positions at which to sample gain pattern
NNN=400; 
probes=source_grid(NNN);

S{1}=sensor_positions(5);
T{1}=0.4*[1,-1;1 1];
C{1}=0.5*[-1,-1; -1,1; 0,0];


iConf=1;

% 2 targets
Mt=1./sqdist(T{iConf},S{iConf});
% multiple competitors
Mc=1./sqdist(C{iConf},S{iConf});
% source signals
nsamples=10000;
a=zeros(1000,1);
b=zeros(1000,1);
a(1:200)=sin(pi*(1:200)/200)';
b(101:400)=sin(pi*(1:300)/300)';
st=repmat([a,b],10,1);  % 10 cycles
sc=randn(nsamples,size(C{iConf},1));
% sensor signals
Xt=st*Mt;
Xc=sc*Mc;
% adjust amplitude of target for SNR=1 on best sensor
pwrt=mean(Xt.^2);
pwrc=mean(Xc.^2);
r=sqrt(max(pwrt./pwrc));
Xt=Xt/r;
% target SNR
SNR=sqrt(1);
Xt=Xt*SNR;
% sensor signals
X=Xc+Xt;


X=X+0.0001*randn(size(X)); % to avoid singularity

% epoch, apply DSS
X=nt_fold(X,nsamples/10);
todss=nt_dss1(X);
z=nt_mmat(X,todss);

figure(1); clf;

subplot 251
a=mean(nt_fold(st,nsamples/10),3);
a=a/max(abs(a(:)));
plot(a);
ylim(1.1*[-1 1]);
set(gca,'fontsize',14,'ytick',[]);
xlabel('samples')
title('sources')
legend('source 1', 'source 2', 'location', 'south'); legend boxoff
set(gca,'xtick',[1 1000]);

subplot 256
plot(mean(X,3));
%ylim(1.1*[-1 1]);
set(gca,'fontsize',14,'ytick',[]);
xlabel('samples')
title('mix')
set(gca,'xtick',[1 1000]);

subplot 252
a=mean(nt_fold(z(:,1:2,:),nsamples/10),3);
a=a/max(abs(a(:)));
[~,idx]=max(abs(a));
plot(a(:,1)*sign(a(idx(1),1)), 'b'); hold on
plot(a(:,2)*sign(a(idx(2),2)), 'r');
ylim(1.6*[-1 1]);
set(gca,'fontsize',14,'ytick',[]);
xlabel('samples')
title('JD 1');
legend('component 1', 'component 2', 'location', 'north'); legend boxoff
set(gca,'xtick',[1 1000]);

for iComp=1:2
    f=todss(:,iComp); % virtual electrode    
    
    % display gain pattern
    forward=source2sensor(probes,S{iConf}); 
    forward=reshape(forward,NNN*NNN,size(S{iConf},1));
    F=f;
    GF=reshape(forward*F,NNN,NNN,1); % grid to filteroutput
    gain_pattern=GF(:,:,1);
    
    
    mask=zerox(gain_pattern,1); % set pixels near zero to zero
    [gain_pattern,ticks,ticklabels]=symlog(gain_pattern*100,[],0.1); % --> signed log 
    gain_pattern=gain_pattern.*mask; % set zero set to zero
    gain_pattern=max(-3.9,gain_pattern);
    gain_pattern=zero_outside(gain_pattern,-4);

    subplot(2,5,iComp+6);

    imagesc([-1 1], [-1, 1],gain_pattern');
    viscircles(T{iConf},0.05, 'color', 'k','linestyle','-', 'linewidth', 1.5, 'enhancevisibility',0); 
    drawcross(C{iConf},'b',0.1);
    
    % cmap=turbo;
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
    h=scatter(S{iConf}(:,1)*a,S{iConf}(:,2)*a, 300, 'k', 'filled');hold on
    h=scatter(S{iConf}(:,1)*a,S{iConf}(:,2)*a, 200, 'w', 'filled');
    xlim([-a,a]); ylim([-a,a]);
    set(gca,'clim',[-4 4]);
    set(gca,'fontsize',14);
    %    plot_tweak([0 0.06 0 -0.06], h7)
    h=colorbar('southoutside');
    set(h,'ticks', [-2 -1 0 1 2], 'ticklabels',{'-100','-10', '0','10', '100'},'fontsize',14, 'limits',[-2.5 2.5], 'fontsize',14)
    set(get(h,'label'), 'string', 'gain');
    
    %plot_tweak([0 -.07 0 .05])
    h=title(['filter ',num2str(iComp)]);
    set(h,'position', [1.0573e-05 -1.2 0])
end


% figure(2); clf
% set(gcf,'position',[300   297   1200 300])


subplot 254
c0=nt_cov(z(:,1:2,:))/1000;
c1=nt_cov(z(1:100,1:2,:))/200;
todss2=nt_dss0(c0,c1);
zz=nt_mmat(z(:,1:2,:),todss2);
a=mean(nt_fold(zz,nsamples/10),3);
a=a/max(abs(a(:)));
[~,idx]=max(abs(a));
plot(a(:,1)*sign(a(idx(1),1)), 'b:'); hold on
plot(a(:,2)*sign(a(idx(2),2)), 'r','linewidth',2);
ylim(1.6*[-1 1]);
set(gca,'fontsize',14,'ytick',[]);
xlabel('samples')
title('JD 2')
legend('component 1''', 'component 2''', 'location','north'); legend boxoff
set(gca,'xtick',[1 1000]);


subplot 255
c0=nt_cov(z(:,1:2,:))/1000;
c1=nt_cov(z(201:400,1:2,:))/400;
todss2=nt_dss0(c0,c1);
zz=nt_mmat(z(:,1:2,:),todss2);
a=mean(nt_fold(zz,nsamples/10),3);
a=a/max(abs(a(:)));
[~,idx]=max(abs(a));
plot(a(:,1)*sign(a(idx(1),1)), 'b:'); hold on
plot(a(:,2)*sign(a(idx(2),2)), 'r','linewidth',2);
ylim(1.6*[-1 1]);
set(gca,'fontsize',14,'ytick',[]);
xlabel('samples')
title('JD 3')
legend('component 1"', 'component 2"', 'location','north'); legend boxoff
set(gca,'xtick',[1 1000]);


subplot 259
c0=nt_cov(z(:,1:2,:))/1000;
c1=nt_cov(z(1:100,1:2,:))/200;
todss2=nt_dss0(c0,c1);
f=todss(:,1:2)*todss2(:,2); % virtual electrode

% display gain pattern
forward=source2sensor(probes,S{iConf}); 
forward=reshape(forward,NNN*NNN,size(S{iConf},1));
F=f;
GF=reshape(sum(forward*F,2),NNN,NNN,1); % grid to filteroutput
gain_pattern=GF(:,:,1);

mask=zerox(gain_pattern,1); % set pixels near zero to zero
[gain_pattern,ticks,ticklabels]=symlog(gain_pattern*100,[],0.1); % --> signed log 
gain_pattern=gain_pattern.*mask; % set zero set to zero
gain_pattern=max(-3.9,gain_pattern);
gain_pattern=zero_outside(gain_pattern,-4);

imagesc([-1 1], [-1, 1],gain_pattern');
viscircles(T{iConf},0.05, 'color', 'k','linestyle','-', 'linewidth', 1.5, 'enhancevisibility',0); 
drawcross(C{iConf},'b',0.1);

% cmap=turbo;
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
h=scatter(S{iConf}(:,1)*a,S{iConf}(:,2)*a, 300, 'k', 'filled');hold on
h=scatter(S{iConf}(:,1)*a,S{iConf}(:,2)*a, 200, 'w', 'filled');
xlim([-a,a]); ylim([-a,a]);
set(gca,'clim',[-4 4]);
set(gca,'fontsize',14);
%    plot_tweak([0 0.06 0 -0.06], h7)
h=colorbar('southoutside');
set(h,'ticks', [-2 -1 0 1 2], 'ticklabels',{'-100','-10', '0','10', '100'},'fontsize',14, 'limits',[-2.5 2.5], 'fontsize',14)
set(get(h,'label'), 'string', 'gain');
%plot_tweak([0 -.07 0 .05])
h=title(['filter 3']);
set(h,'position', [1.0573e-05 -1.2 0])


c0=nt_cov(z(:,1:2,:))/1000;
c1=nt_cov(z(201:400,1:2,:))/400;
todss2=nt_dss0(c0,c1);

subplot (2,5,10)
c0=nt_cov(z(:,1:2,:))/1000;
c1=nt_cov(z(201:400,1:2,:))/400;
todss2=nt_dss0(c0,c1);
% zz=nt_mmat(z(:,1:2,:),todss2);
% a=mean(nt_fold(zz,nsamples/10),3);
% a=a/max(abs(a(:)));
% plot(a);
% ylim(1.1*[-1 1]);
% set(gca,'fontsize',14,'ytick',[]);
% xlabel('samples')
% title('filter 4')
% legend('component 1"', 'component 2"', 'location','northeast'); legend boxoff

f=todss(:,1:2)*todss2(:,2); % virtual electrode

% display gain pattern
forward=source2sensor(probes,S{iConf}); 
forward=reshape(forward,NNN*NNN,size(S{iConf},1));
F=f;
GF=reshape(sum(forward*F,2),NNN,NNN,1); % grid to filteroutput
gain_pattern=GF(:,:,1);

mask=zerox(gain_pattern,1); % set pixels near zero to zero
[gain_pattern,ticks,ticklabels]=symlog(gain_pattern*100,[],0.1); % --> signed log 
gain_pattern=gain_pattern.*mask; % set zero set to zero
gain_pattern=max(-3.9,gain_pattern);
gain_pattern=zero_outside(gain_pattern,-4);

imagesc([-1 1], [-1, 1],gain_pattern');
viscircles(T{iConf},0.05, 'color', 'k','linestyle','-', 'linewidth', 1.5, 'enhancevisibility',0); 
drawcross(C{iConf},'b',0.1);

% cmap=turbo;
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
    h=scatter(S{iConf}(:,1)*a,S{iConf}(:,2)*a, 300, 'k', 'filled');hold on
    h=scatter(S{iConf}(:,1)*a,S{iConf}(:,2)*a, 200, 'w', 'filled');
xlim([-a,a]); ylim([-a,a]);
set(gca,'clim',[-4 4]);
set(gca,'fontsize',14);
%    plot_tweak([0 0.06 0 -0.06], h7)
h=colorbar('southoutside');
set(h,'ticks', [-2 -1 0 1 2], 'ticklabels',{'-100','-10', '0','10', '100'},'fontsize',14, 'limits',[-2.5 2.5], 'fontsize',14)
set(get(h,'label'), 'string', 'gain');
%plot_tweak([0 -.07 0 .05])
h=title(['filter 4']);
set(h,'position', [1.0573e-05 -1.2 0])

subplot 252
plot_tweak([0.08 0 0 0]);

return