% for sizing to work this should be called with a figure generated from
% figurebetter([],[12 12]);
%
% [ecchand,dirhand] = facedist_gridmds(F,rdm,stimuli,[polarmode=false])
function [ecchand,dirhand] = facedist_gridmds(F,rdm,stimuli,polarmode)

if ieNotDefined('polarmode')
    polarmode = false;
end

rdv = asrdmvec(rdm);
% rescore negative values (below-chance discriminants) to 0.
rdv(rdv<0) = 0;
rdm = asrdmmat(rdv);
ncon = size(rdm,1);

facedist_refxy = facedist_getrefxy(ncon);
% rotate 
facedist_refxy = rotatepoints(facedist_refxy',1/2*pi)';
if ncon == 24
    % offset second view
    facedist_refxy(13:24,:) = bsxfun(@plus,facedist_refxy(13:24,:),...
        [range(facedist_refxy(:,1))*1.1 0]);
end

xy = mdscale_robust(rdm,2,'criterion','metricstress');

% try to align mdscale to reference distances
[deviance,xy] = procrustes(facedist_refxy,xy);

colors = facedist_colors('ecc');
[dirrdm,eccrdm] = facedist_rsapredictor_neighborsonly(ncon);

% setup colour and type of marker
facecolors = colors(repinds(1,3,4),:);
if ncon==12
    facemarkers = repmat({''},[12 1]);
    scalef = 1;
else
    facemarkers = {'L';'R'};
    facemarkers = facemarkers(repinds(1,2,12),:);
    scalef = 2;
end

ax = axes('position',[.05 .05 .9 .9]);
set(ax,'dataaspectratio',[1 1 1]);

args = arrayfun(@(x){'linewidth',2,'marker','o','markersize',10,...
    'markerfacecolor',colors(x,:),'markeredgecolor',colors(x,:)},...
    3:-1:1,'uniformoutput',0);

for v = 1:size(dirrdm,3)
    dirhand{v} = plotmdslines(ax,xy,dirrdm(:,:,v),[0 0 0],false,...
        'linewidth',1,'linestyle',':');
    ecchand{v} = plotmdslines(ax,xy,eccrdm(:,:,v),colors,false,args{:});
    uistack(ecchand{v},'bottom');
    setplotz(ecchand{v},-1);
    uistack(dirhand{v},'bottom');
    setplotz(dirhand{v},-1);
end
axis(ax,'off');

if ieNotDefined('stimuli')
    thand = textscatter(xy,facemarkers,[1 1 1]);
    set(thand,'fontsize',8);
    uistack(thand,'top');
    maxax = max(range(xy));
    % 10% increase in longest dim
    axsc = maxax / 20;
    % need to do our own scaling
    axis(ax,[min(xy(:,1)) max(xy(:,1)) min(xy(:,2)) max(xy(:,2))] .* 1.3);
else
    % kill markers
    cellfun(@(thish)set(thish,'marker','none'),[dirhand,ecchand]);
    h = scalef * (range(xy(:,2))/ncon);
    outh = imageaxes(gca,xy,{stimuli.image},{stimuli.alpha},h);
    uistack(outh,'top');
end

if polarmode
    % kill lines (we plot our own)
    cellfun(@(thish)set(thish,'linestyle','none'),[dirhand,ecchand]);
end

padaxislims('ax',gca,'axdir','x');
padaxislims('ax',gca,'axdir','y');
