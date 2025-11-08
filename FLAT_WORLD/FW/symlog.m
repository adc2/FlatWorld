function [y, ticks, ticklabels]=symlog(x,base, minabs)
%[y,ticks,ticklabels]=symlog(x,base, minabs) - symmetric comppressive log transform
%  
%  y: compressed data (-inf to +inf)
%
%  x: data to compress (-inf to inf)
%  base: base for log [default: 10]
%  minabs: minimum absolute value
%

if nargin<2||isempty(base); base=10; end
if nargin<3||isempty(minabs); minabs=0;end

y=zeros(size(x));
idx=x>0;
y(idx)=max(minabs,log(x(idx))/log(base));
idx=x<0;
y(idx)=min(-minabs,-log(-x(idx))/log(base));

idx=x==0;
y(idx)=0;

ticks=ceil(min(y(:))) : floor(max(y(:)));
ticklabels=[-base.^-ticks(ticks<0), 0, base.^ticks(ticks>0)];

