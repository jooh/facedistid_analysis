% handles = facedist_abseccplot(groupres,meanres,customfits,vstruct,roi)
function handles = facedist_abseccplot(groupres,meanres,customfits,vstruct,roi,chancelevel)

par = facedist_plotpar;
if chancelevel==50
    % facepairs mode
    par.ylabel = par.ylabel_facepairs;
end

glmmode = false;
gtarget = 'r';
if any(strcmp(meanres.rows_contrast,'face01'))
    % GLM mode
    par.ylabel = par.ylabel_glm;
    glmmode = true;
    gtarget = 'b';
end

fsize = par.abseccfsize;
handles.f = figurebetter(fsize);
handles.ax = gca(handles.f);

% find the relevant data
roiind = strcmp(groupres.cols_roi,roi);
assert(sum(roiind)==1);

for v = 1:numel(vstruct)
    customind = strcmp({customfits.name},[vstruct(v).prefix 'slope_absecc']);
    assert(sum(customind)==1);
    customconind{v} = strcmp(groupres.rows_contrast,[vstruct(v).prefix ...
        'slope_absecc']);
    assert(sum(customconind{v})==1);
    connames = customfits(customind).cons;
    [~,conind{v}] = intersect(groupres.rows_contrast,connames,'stable');
    assert(numel(conind{v})==3);
    %[~,slopeind{v}] = intersect(groupres.rows_contrast,...
        %customfits(customind).name);
    %assert(numel(slopeind{v})==1);
end
nsub = numel(groupres.z_subject);

% set up plot
subx = linspace(par.subxscale(1),par.subxscale(2),nsub);
xl = [.5 3.5];
textx = xl(2) + range(xl)*.05;
xvals = linspace(1,3,1000);
nboot = 0;

if numel(vstruct)==1
    viewlabel = vstruct.xlabel2;
else
    viewlabel = '';
end
set(handles.f,'name',['absecc ' viewlabel ' ' roi]);
hold(handles.ax,'on');
pmat = [];

vstyles = struct('marker',{par.markertype,'s'},'line',{'-','--'},...
    'facecolor',{par.basecolors,ones(size(par.basecolors))},...
    'xoff',{par.viewxoff(1),par.viewxoff(2)},...
    'textcolor',{par.basecolors,par.basecolors});
if numel(vstruct)==1
    [vstyles.xoff] = deal(0);
end

for v = 1:numel(vstruct)
    xmat = [];
    datamat = [];
    for x = 1:3
        ecolor = par.basecolors(x,:); % keep track of data
        datamat(:,x) = squeeze(groupres.(gtarget)(conind{v}(x),roiind,:));
        if isfield(groupres,'pperm') && ~glmmode
            pmat(:,x,v) = squeeze(groupres.pperm(conind{v}(x),roiind,:));
        end
        nanmat = ~isnan(datamat(:,x));
        xmat(1:size(nanmat,1),x) = x;
        % GLM for bootstrap later
        model{v}(x) = GLM([xmat(nanmat,x) ones(sum(nanmat),1)],...
            datamat(nanmat,x));
        % data points
        handles.points(:,x,v)= plotarray(subx+x+vstyles(v).xoff,...
            datamat(:,x),...
            vstyles(v).marker,'markersize',par.markersize,...
            'markerfacecolor',vstyles(v).facecolor(x,:),...
            'markeredgecolor',ecolor,'linewidth',1);
        % mean
        handles.meanerrors(x,v) = errorbar2(x+vstyles(v).xoff,...
            nanmean(datamat(:,x)),...
            nansterr(datamat(:,x)),1,'-','color',par.basecolors(x,:),...
            'linewidth',par.errorbarlinewidth);
        handles.meanpoints(x,v) = plot(x+vstyles(v).xoff,...
            nanmean(datamat(:,x)),...
            vstyles(v).marker,'markersize',par.markersize,...
            'markerfacecolor',vstyles(v).facecolor(x,:),...
            'markeredgecolor',ecolor,'linewidth',1);
    end
    % get the bootstrapped parameter estimates
    bootind = preparesampleboots(model{v},nboot);
    if nboot > 0 && ~isempty(bootind)
        % bootstrap fit
        bootdist = bootstrapsamples(model{v},bootind,@booter,[],xvals);
        % now we can obtain the bootstrap errors like so
        bootest(:,:,v) = prctile(squeeze(bootdist),[5 50 95],2);
    else
        % plot LS line instead of median BS estimate
        bootest(:,2,v) = booter(model{v},xvals);
    end
    handles.slope(v) = plot(handles.ax,xvals+vstyles(v).xoff,...
        squeeze(bootest(:,2,v)),...
        vstyles(v).line,...
        'color',par.fitlinecolor,'linewidth',par.fitlinewidth);
    if isfield(meanres,'ppara')
        rstr = char(p2str(...
            meanres.ppara(customconind{v},roiind),3,true));
        textcolor = [0 0 0];
        fontweight = 'normal';
        fontangle = 'italic';
        if meanres.ppara(customconind{v},roiind)<0.05
            fontweight = 'bold';
            fontangle = 'normal';
            % textcolor = par.notsigcolor;
        end
        [~,ind] = min(abs(xvals-xl(2)));
        ry = bootest(ind,2,v);
        handles.tr(v) = text(textx,ry,rstr,'horizontalalignment','left',...
            'verticalalignment','middle','color',textcolor,'edgecolor',...
            'none','fontweight',fontweight,'fontangle',fontangle,...
            'linestyle',vstyles(v).line,'linewidth',1);
        yl = ylim(handles.ax);
        for x = 1:3
            handles.tp(x,v) = addptext(x+vstyles(v).xoff,yl(2),...
                meanres.ppara(conind{v}(x),roiind),3,'p=','rotation',45,...
                'color',vstyles(v).textcolor(x,:));
        end
        % put at bottom (avoid background colour blocking everything)
        arrayfun(@(x)uistack(x,'bottom'),withoutnan(handles.tp(:,v)));
    else
        handles.tr(v) = NaN;
        handles.tp(1:3,v) = NaN;
    end
end

% and plot these
if numel(vstruct)>1
    % now add NaNs to offset instead
    
    % change in x per entry
    xpertick = range(xvals) ./ numel(xvals);
    % total x range (before allowing for shifts)
    xshift = range([vstyles.xoff]);
    % so before the beginning of the line we need
    xbase = min(xvals);
    xa = [xbase+vstyles(1).xoff:xpertick:xbase-xpertick];
    % and after the line we need
    xend = max(xvals);
    xb = [xend+xpertick:xpertick:xend+vstyles(2).xoff];
    % and the full xvals is
    newx = [xa xvals xb]';

    % number of samples to shift by
    noff = numel(xvals) ./ range(xvals) * xshift;
    ntoadd = ceil(noff);
    if nboot > 0
        low = [NaN([ntoadd 2]); squeeze(bootest(:,1,:))];
        high = [NaN([ntoadd 2]); squeeze(bootest(:,3,:))];
        low(:,1) = circshift(low(:,1),[-ntoadd 0]);
        high(:,1) = circshift(high(:,1),[-ntoadd 0]);
        handles.shade = errorshade(newx,low,high,par.shadecolor([1 1],:));
    else
        handles.shade = NaN;
    end
else
    if nboot > 0
        handles.shade = errorshade(xvals',...
            bootest(:,1),...
            bootest(:,3),...
            par.shadecolor);
    else
        handles.shade = NaN;
    end
end


arrayfun(@(x)uistack(x,'bottom'),handles.slope);
if chancelevel~=50
    handles.chance = line([0 4],[chancelevel chancelevel],...
        'linestyle',par.chancelinestyle,...
        'linewidth',par.chancelinewidth,'color',par.chancelinecolor);
    % add p values for each ecc level, R2 estimate with CI
    handles.tchance = text(textx,chancelevel,'chance',...
        'horizontalalignment','left',...
        'verticalalignment','middle');
    uistack(handles.chance,'bottom');
end

xlim(handles.ax,xl);

% final styling
set(handles.ax,'xtick',1:3,'xticklabel',par.absecclabels);
ylabel(handles.ax,par.ylabel);
xlabel(handles.ax,{par.abseccxlabel,viewlabel});
set(handles.ax,'position',[repmat(par.axscale,[1 4])],'tickdir','out',...
    'ticklength',par.ticklength);

if nboot > 0
    arrayfun(@(x)uistack(x,'bottom'),handles.shade);
end

% chance shade
if ~par.dochanceshade || ~isfield(meanres,'ppara') || all(isnan(meanres.ppara(:)))
    [handles.chanceshade,handles.chancetext] = deal(NaN);
else
    ydata = arrayfun(@(x)getn(x,'ydata'),handles.points,'uniformoutput',1);
    xdata = arrayfun(@(x)getn(x,'xdata'),handles.points,'uniformoutput',1);
    [handles.chanceshade,handles.chancetext] = facedist_chanceshadewrap(...
        xdata,ydata,textx);
end

function x = booter(bootmod,xvals)

est = fit(bootmod);
x = est(2) + xvals * est(1);
