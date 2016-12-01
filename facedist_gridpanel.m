%
% [spokes,rings] = facedist_gridpanel(ax,assumesymmetry,custommask,scalefactor)
function [spokes,rings] = facedist_gridpanel(ax,assumesymmetry,custommask,scalefactor)

[phi,rad] = facedist_stimcoords(assumesymmetry);

if ~ieNotDefined('scalefactor')
    rad = rad*scalefactor;
end

ax = gca;
hold(ax,'on');
nanmask = [];
if ~assumesymmetry
    % mask out left side
    nanmask = [.5*pi 1.5*pi];
end

if ~ieNotDefined('custommask')
    nanmask = custommask;
end

[spokes,rings] = polargrid(ax,phi,rad,nanmask);
setplotz(spokes,-1);
setplotz(rings,-.25);
arrayfun(@(x)uistack(x,'bottom'),spokes);
colors = facedist_colors('ecc');
for s = 1:numel(rings)
    set(rings(s),'color',colors(s,:),'linewidth',2);
end
set(spokes,'linewidth',1,'linestyle',':','color',[0 0 0]);
