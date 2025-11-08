function drawcross(positions,color,sz,open,varargin)
if nargin<3||isempty(sz); sz=0.1; end
if nargin<4||isempty(open); open=0; end

if numel(sz)==1; sz=[sz;sz]; end
for iSource=1:size(positions,1)
    x=positions(iSource,1); y=positions(iSource,2);
    if open
        line([x+sz(1),x+sz(1)/2],[y, y], 'color', color, varargin{:});
        line([x,x],[y+sz(2), y+sz(2)/2], 'color', color, varargin{:});
        line([x-sz(1),x-sz(1)/2],[y, y], 'color', color, varargin{:});
        line([x,x],[y-sz(2), y-sz(2)/2], 'color', color, varargin{:});
    else
        line([x+sz(1),x-sz(1)],[y, y], 'color', color, varargin{:});
        line([x,x],[y+sz(2), y-sz(2)], 'color', color, varargin{:});
    end
    end % function drawcross
end
