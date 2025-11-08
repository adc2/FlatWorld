function plotsensors(sensors)

N=100;
x=sin(2*pi*(1:N+1)/N);
y=cos(2*pi*(1:N+1)/N);
plot(x,y, 'k', 'linewidth', 2);
a=1.05;
plot(sensors(:,1)*a,sensors(:,2)*a, '.', 'color', [0,0.2,1],'markersize',40);
xlim([-a,a]); ylim([-a,a]);
