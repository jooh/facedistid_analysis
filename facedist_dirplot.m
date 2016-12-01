% handles = facedist_dirplot(groupres,meanres,customfits,vstruct,roi)
function handles = facedist_dirplot(groupres,meanres,customfits,vstruct,roi,chancelevel)

par = facedist_plotpar;
doingblendedshade = false;
glmmode = false;
if chancelevel==50
    % facepairs mode
    par.ylabel = par.ylabel_facepairs;
end
gtarget = 'r';
if any(strcmp(meanres.rows_contrast,'face01'))
    % GLM mode
    par.ylabel = par.ylabel_glm;
    par.dirxlabel = par.dirxlabel_glm;
    % later we may need this to hack around some old hacks
    glmmode = true;
    gtarget = 'b';
end

% set up plot
fsize = par.dirfsize;
% black for mean
par.basecolors(4,:) = [0 0 0];
handles.f = figurebetter(fsize);
handles.ax = gca(handles.f);
% nb we may or may not plot a datapoint at x 1 (depends on whether we are
% within or across views)
% if in the final analysis we don't want to use this analysis, we probably
% want to forget this and plot with the same settings as ecc here.
xl = [.5 4.5];
if numel(vstruct)==1
    viewlabel = vstruct.xlabel2;
else
    viewlabel = '';
end
set(handles.f,'name',['dir ' viewlabel ' ' roi]);
hold(handles.ax,'on');

% find the relevant data
roiind = strcmp(groupres.cols_roi,roi);
assert(sum(roiind)==1);

subind = ~isnan(squeeze(groupres.(gtarget)(1,roiind,:)));
nsub = sum(subind);
subx = linspace(par.subxscale(1),par.subxscale(2),nsub);
nboot = 0;
bootind = bootindices(nsub,nboot);

if glmmode
    % something entirely different - straight line fit through all 4 points
    xvals = linspace(1,4,1000)';
    xdist = xvals;
    valdist = [1:4]';
else
    % this function describes how distance increases with angle (in radians)
    dfun = @(xv)sin(xv/2)*2;

    % these are the (idealized) test angles in radians
    valrad = [pi/3 pi/3*2 pi]';
    % and in distance units (where we fit a straight line)
    valdist = dfun(valrad);

    % now our x values are going to have to be transformed similarly
    xrad = linspace(pi/3,pi,1000)';
    xdist = dfun(xrad);
    % rescale to 2:4 range
    xvals = (xrad/max(xrad)) * 3 + 1;
end

nview = numel(vstruct);

% ecc by nsub by dir
handles.points = NaN([4 nsub 4 nview]);
handles.meanpoints = NaN([4 4 nview]);
handles.meanerrors = NaN([4 4 nview]);

% stylings for the two viewpoints (within, across)
viewstyles = struct('marker',{par.markertype,'s'},'line',{'-','--'},...
    'facecolor',{par.basecolors,ones(size(par.basecolors))},...
    'xoff',{par.viewxoff(1),par.viewxoff(2)},...
    'textcolor',{par.basecolors,par.basecolors});
if numel(vstruct)==1
    [viewstyles.xoff] = deal(0);
end

for v = 1:numel(vstruct)
    % for each eccentricity level (with mean as the hacky 4th entry)
    for ecc = 1:4
        if ecc < 4
            thisecc = [par.absecclabels{ecc} '-' par.absecclabels{ecc}];
            thisxoff = par.eccxoff(ecc);
            conname = [vstruct(v).prefix 'slope_' thisecc];
        else
            thisecc = 'mean';
            thisxoff = 0;
            conname = [vstruct(v).prefix 'slope_no_eccentricity_tuning'];
        end
        customind = strcmp({customfits.name},conname);
        assert(sum(customind)==1);
        if ecc < 4 || glmmode
            connames = customfits(customind).cons;
            if strcmp(vstruct(v).prefix,'view_across_')
                % insert across view - we want to plot this like any other
                % EXCEPT we don't fit a slope.
                connames = [{[vstruct(v).prefix thisecc '_direction_11_m0.0']};...
                    connames];
            end
        else
            connames = arrayfun(@(x)sprintf('%smean_direction_%02d',...
                vstruct(v).prefix,x),[11 17 25 28],'uniformoutput',0);
            if ~strcmp(vstruct(v).prefix,'view_across_')
                connames(1) = [];
            end
        end
        customconind = strcmp(groupres.rows_contrast,conname);
        assert(sum(customconind)==1);
        [~,conind] = intersect(groupres.rows_contrast,connames,'stable');
        if isfield(meanres,'ppara')
            slopep = meanres.ppara(customconind,roiind);
            chancep = meanres.ppara(conind,roiind);
        else
            slopep = NaN;
            chancep = NaN([3 1]);
        end
        % keep track of data - squeeze out roi dim to make nsub by ncon
        % matrix
        fulldatamat = squeeze(groupres.(gtarget)(conind,roiind,subind))';
        fullpmat = [];
        if isfield(groupres,'pperm') && ~glmmode
            fullpmat = squeeze(groupres.pperm(conind,roiind,subind))';
        end
            
        if ~strcmp(vstruct(v).prefix,'view_across_') && ~glmmode
            fulldatamat = [NaN([size(fulldatamat,1),1]) fulldatamat];
            chancep = [NaN; chancep];
            if isfield(groupres,'pperm') && ~glmmode
                fullpmat = [NaN([size(fulldatamat,1),1]) fullpmat];
            end
        end

        % for each direction
        % (note that we plot everything 1 x offset from this to leave room for
        % delta 0 case below)
        xmat = repmat(valdist',[size(fulldatamat,1) 1]);
        for x = 1:4
            datamat = fulldatamat(:,x);
            if ~all(isnan(datamat))
                % data points
                ecolor = par.basecolors(ecc,:);
                handles.points(ecc,:,x,v) = plotarray(...
                    subx+x+thisxoff+viewstyles(v).xoff,datamat,...
                    viewstyles(v).marker,'markersize',par.markersize,...
                    'markerfacecolor',viewstyles(v).facecolor(ecc,:),...
                    'markeredgecolor',ecolor);
                % mean
                handles.meanerrors(ecc,x,v) = errorbar2(...
                    x+thisxoff+viewstyles(v).xoff,mean(datamat),...
                    sterr(datamat),1,'-','color',...
                    min([viewstyles(v).facecolor(ecc,:); ...
                    viewstyles(v).textcolor(ecc,:)]),...
                    'linewidth',par.errorbarlinewidth);
                handles.meanpoints(ecc,x,v) = plot(...
                    x+thisxoff+viewstyles(v).xoff,mean(datamat),...
                    viewstyles(v).marker,'markersize',par.markersize,...
                    'markerfacecolor',viewstyles(v).facecolor(ecc,:),...
                    'markeredgecolor',ecolor);
            end
            % GLM for bootstrap later
            if x~=1
                model(x) = GLM([xmat(:,x-1) ones(numel(datamat),1)],...
                    datamat);
            end
        end % x = 1:4
        xlim(handles.ax,xl);

        % get the bootstrapped parameter estimates
        % bootstrap fit
        model(1) = [];
        % NB we fit with the distance test values, but we plot the fit with
        % the radians
        bootdist = bootstrapsamples(model,bootind,@booter,[],xdist);
        % now we can obtain the bootstrap errors like so
        bootest = prctile(squeeze(bootdist),[5 50 95],2);
        low(:,ecc,v) = bootest(:,1);
        high(:,ecc,v) = bootest(:,3);
        if isempty(bootind)
            % plot LS line instead of median BS estimate
            bootest(:,2) = booter(model,xdist);
        end
        handles.slope(ecc,v) = plot(handles.ax,xvals+thisxoff+viewstyles(v).xoff,...
            bootest(:,2),viewstyles(v).line,...
            'color',morph(par.fitlinecolor,par.basecolors(ecc,:),.5),...
            'linewidth',par.fitlinewidth);
        % add p values for each ecc level, R2 estimate 
        textx = xl(2) + range(xl)*.05;
        if isfield(meanres,'ppara')
            rstr = char(p2str(slopep,3,true));

            [~,ind] = min(abs(xvals-xl(2)));

            textc = min([viewstyles(v).facecolor(ecc,:); ...
                        viewstyles(v).textcolor(ecc,:)]);
            % thresholded
            fontweight = 'normal';
            fontangle = 'italic';
            if slopep<0.05
                fontweight = 'bold';
                fontangle = 'normal';
            end

            ry = bootest(ind,2);
            handles.tr(ecc,v) = text(textx,ry,rstr,...
                'horizontalalignment','left',...
                'verticalalignment','middle','color',textc,...
                'edgecolor','none',...
                'linestyle',viewstyles(v).line,...
                'linewidth',1,'fontweight',fontweight,...
                'fontangle',fontangle);

            bgcolor = viewstyles(v).textcolor(ecc,:);
            if all(bgcolor==[1 1 1])
                bgcolor = 'none';
            end

            handles.tp(ecc,:,v) = addptext([1:4]+thisxoff+viewstyles(v).xoff,...
                repmat(max(ylim(handles.ax)),...
                [1 4]),chancep,3,'p=','rotation',45,...
                'color',...
                viewstyles(v).textcolor(ecc,:));

            % just plot the errorshades here. We will alpha blend offline.
            if ~doingblendedshade
                if all(isnan(low(:,ecc,v))) || nboot<1
                    handles.shade(ecc,v) = NaN;
                else
                    handles.shade(ecc,v) = errorshade(...
                        xvals+thisxoff+viewstyles(v).xoff,low(:,ecc,v),...
                        high(:,ecc,v),par.basecolors(ecc,:));
                end
            end
        else
            % fill in NaNs in handles
            handles.shade(ecc,v) = NaN;
            handles.tp(ecc,:,v) = NaN;
            handles.tr(ecc,v) = NaN;
        end
    end % ecc = 1:4
    % put at bottom (avoid background colour blocking everything)
    arrayfun(@(x)uistack(x,'bottom'),withoutnan(handles.tp(:,:,v)));
    if strcmp(vstruct(v).prefix,'view_within_')
        delete(handles.tp(:,1,v));
        handles.tp(:,1,v) = NaN;
    end
end

% plot the error shades - right, for now we are going to give up on
% blending these
if doingblendedshade
    % NB currently does not shift with the eccentricity groups
    if numel(vstruct)>1 && nboot >0
        % change in x per entry
        xpertick = range(xvals) ./ numel(xvals);
        % total x range (before allowing for shifts)
        xshift = range([viewstyles.xoff]);
        % so before the beginning of the line we need
        xbase = min(xvals);
        xa = [xbase+viewstyles(1).xoff:xpertick:xbase-xpertick];
        % and after the line we need
        xend = max(xvals);
        xb = [xend+xpertick:xpertick:xend+viewstyles(2).xoff];
        % and the full xvals is
        newx = [xa; xvals; xb];
        % number of samples to shift by
        noff = numel(xvals) ./ range(xvals) * xshift;
        ntoadd = ceil(noff);
        low = [NaN([ntoadd 4 2]); low];
        high = [NaN([ntoadd 4 2]); high];
        low(:,:,1) = circshift(low(:,:,1),[-ntoadd 0 0]);
        high(:,:,1) = circshift(high(:,:,1),[-ntoadd 0 0]);
        handles.shade = errorshade(newx,cat(2,low(:,1:3,1),low(:,1:3,2)),...
            cat(2,high(:,1:3,1),high(:,1:3,2)),...
            par.basecolors([1:3 1:3],:));
        handles.mshade = errorshade(newx,squeeze(low(:,4,:)),...
            squeeze(high(:,4,:)),par.shadecolor([1 1],:));
    else
        % nice and simple
        if nboot > 0
            handles.shade = errorshade(xvals,low(:,1:3),high(:,1:3),...
                par.basecolors(1:3,:));
            handles.mshade = errorshade(xvals,low(:,4),high(:,4),[.8 .8 .8]);
        else
            handles.shade = NaN;
            handles.mshade = NaN;
        end
    end
else
    handles.mshade = handles.shade(end,:);
    handles.shade(end,:) = [];
end

arrayfun(@(thish)uistack(thish,'bottom'),handles.slope);
if chancelevel~=50
    handles.chance = line([0 4.5],[chancelevel chancelevel],...
        'linestyle',par.chancelinestyle,...
        'linewidth',par.chancelinewidth,'color',par.chancelinecolor);
    handles.tchance = text(textx,chancelevel,'chance',...
        'horizontalalignment','left',...
        'verticalalignment','middle');
    uistack(handles.chance,'bottom');
end

if ~isnan(handles.mshade)
    arrayfun(@(x)uistack(handles.shade(x),'bottom'),1:numel(handles.shade));
    arrayfun(@(x)uistack(handles.mshade(x),'bottom'),1:numel(handles.mshade));
end
yl = ylim(handles.ax);

% final styling
set(handles.ax,'xtick',1:4,'xticklabel',par.dirlabels);
ylabel(handles.ax,par.ylabel);
xlabel(handles.ax,{par.dirxlabel,viewlabel});
set(handles.ax,'position',[repmat(par.axscale,[1 4])],'tickdir','out',...
    'ticklength',par.ticklength);

% transplant over the other m series for ease of indexing later
handles.mpoints = handles.points(4,:,:,:);
handles.points(4,:,:,:) = [];
handles.mslope = handles.slope(4,:);
set(handles.mslope,'color',par.fitlinecolor);
handles.slope(4,:) = [];
handles.mtr = handles.tr(4,:);
handles.tr(4,:) = [];
handles.mtp = handles.tp(4,:,:);
handles.tp(4,:,:) = [];
handles.mmeanpoints = handles.meanpoints(4,:,:);
handles.meanpoints(4,:,:) = [];
handles.mmeanerrors = handles.meanerrors(4,:,:);
handles.meanerrors(4,:,:) = [];

% chance shade. Here a slight wrinkle is that we need to have two shades
% basically - one for grouped data and one for means.
if ~par.dochanceshade || ~isfield(meanres,'ppara') || all(isnan(meanres.ppara(:)))
    [handles.chanceshade,handles.chancetext,handles.mchanceshade,...
        handles.mchancetext] = deal(NaN);
else
    xm = arrayfun(@(x)getn(x,'xdata'),handles.mpoints);
    ym = arrayfun(@(x)getn(x,'ydata'),handles.mpoints);
    [handles.mchanceshade,handles.mchancetext] = ...
        facedist_chanceshadewrap(squeeze(xm),squeeze(ym),textx);
    xg = shiftdim(arrayfun(@(x)getn(x,'xdata'),handles.points),1);
    yg = shiftdim(arrayfun(@(x)getn(x,'ydata'),handles.points),1);
    [handles.chanceshade,handles.chancetext] = facedist_chanceshadewrap(...
        xg,yg,textx);
end

function x = booter(bootmod,xvals)

est = fit(bootmod);
x = est(2) + xvals * est(1);
