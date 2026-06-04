function ha = TightPlots(Nh, Nw, w, AR, gap, marg_h, marg_w, units)


% Default values if no input
if nargin<3; w = 15; end
if nargin<4; AR = [2 1]; end
if nargin<5 || isempty(gap); gap = 0.8; end
if nargin<6 || isempty(marg_h); marg_h = [ 0.8 0.4 ]; end
if nargin<7 || isempty(marg_w); marg_w = [ 0.8 0.4 ]; end
if nargin<8; units = 'centimeters'; end
if numel(gap)==1; gap = [gap gap]; end
if numel(marg_w)==1; marg_w = [marg_w marg_w]; end
if numel(marg_h)==1; marg_h = [marg_h marg_h]; end

% Ensure appropriate unit input
if strcmp(units, 'centimeters') + strcmp(units, 'inch') + ...
        strcmp(units, 'points') + strcmp(units, 'pixels') == 0;
    error('Units must be centimeters ''inch'' ''points'' or ''pixels''');
end

% Unit conversion to points
if strcmp(units, 'centimeters') == 1; 
    con = 72/2.54;
elseif strcmp(units, 'inch') == 1
    con = 72;
elseif strcmp(units, 'pixels') == 1
    con = 0.75;
else
    con = 1;
end

units = 'points';
w = w * con;
gap = gap * con;
marg_h = marg_h * con;
marg_w = marg_w * con;

% Calculation of axis width and height
axw = (w-sum(marg_w)-(Nw-1)*gap(2))/Nw;
axh = AR(2)/AR(1)*axw;

% Calculation of figure height
h = Nh*axh + (Nh-1)*gap(1) + sum(marg_h);

% Obtain screen dimensions to place figure in the centre
set(0, 'units', units);
screensize = get(0, 'screensize');
figSize = [ screensize(3)/2 - w/2  screensize(4)/2 - h/2 w h];

% Set the same figure dimensions on screen and paper.
set(gcf, 'Units', units, 'Resize', 'off')
set(gcf, 'Position', figSize)
set(gcf, 'PaperUnits', units, 'PaperSize', [w h])
set(gcf, 'PaperPositionMode', 'manual', 'PaperPosition', [0 0 w h])

% Build axes in figure space
py = h-marg_h(2)-axh; 

ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Nh
    px = marg_w(1);
    
    for ix = 1:Nw
        ii = ii+1;
        
        ha(ii) = axes;
        set(gca, 'Units', units, 'Position', [px py axw axh])
        set(gca, 'XTickLabel', '', 'YTickLabel', '')
        
        px = px+axw+gap(2);
    end
    
    py = py-axh-gap(1);
end

end
