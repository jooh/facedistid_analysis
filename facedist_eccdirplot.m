% handles = facedist_eccdirplot(groupres,meanres,customfits,vstruct,roi,chancelevel)
% TODO:
% RDM legend with same colors as here. NB we have predictors in the parent
% function...
function handles = facedist_eccdirplot(groupres,meanres,customfits,vstruct,roi,chancelevel)

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

fsize = par.eccdirfsize;
handles.f = figurebetter(fsize);
handles.ax = gca(handles.f);

ffxtest = false;

% find the relevant data
roiind = strcmp(groupres.cols_roi,roi);
assert(sum(roiind)==1);

viewstyles = struct('marker',{par.markertype,'s'},'line',{'-','--'},...
    'facecolor',{par.basecolors,ones(size(par.basecolors))},...
    'xoff',{par.viewxoff(1),par.viewxoff(2)},...
    'textcolor',{ones(size(par.basecolors)),par.basecolors});
if numel(vstruct)==1
    [viewstyles.xoff] = deal(0);
end

ecccolors = facedist_colors('ecc');
eccnames = {'sub','typical','caricature'};
eccstruct = struct('group',{'sub-sub','typical-typical',...
    'caricature-caricature','sub-typical','typical-caricature',...
    'sub-caricature'},'delta',num2cell([0 0 0 1 1 2]),'x',[],'y',[],'p',[],...
    'view',[],'pslope',[]);
necc = numel(eccstruct);
eccstruct = repmat(eccstruct,[numel(vstruct),1]);


mfield = 'mean';
if ~isfield(meanres,mfield)
    mfield = 'r';
    assert(isfield(meanres,mfield),'cannot find correct meanres field');
end

hasp = isfield(meanres,'ppara');

for v = 1:numel(vstruct)
    [eccstruct(v,:).view] = deal(vstruct(v).prefix);
    for ec = 1:necc
        eccstruct(v,ec).color = facedist_colors(eccstruct(v,ec).group);
        pattern = sprintf('%sslope_%s',vstruct(v).prefix,...
            eccstruct(v,ec).group);
        % index into customfits
        customind = strcmp({customfits.name},pattern);
        assert(sum(customind)==1);
        % index into meanres
        customconind{v,ec} = strcmp(meanres.rows_contrast,pattern);
        assert(sum(customconind{v,ec})==1);
        connames = customfits(customind).cons;
        [~,conind{v,ec}] = intersect(meanres.rows_contrast,connames,...
            'stable');
        assert(numel(conind{v,ec})==3 || numel(conind{v,ec})==4);
        % index into mean for this group
        meanind = strfindcell(meanres.rows_contrast,sprintf('%s%s_m',...
            vstruct(v).prefix,eccstruct(v,ec).group));
        assert(numel(meanind)==1);
        % just save the x and y data at this point instead.
        eccstruct(v,ec).x = customfits(customind).x';
        eccstruct(v,ec).xmean = mean(eccstruct(v,ec).x);
        eccstruct(v,ec).y = meanres.(mfield)(conind{v,ec},roiind)';
        eccstruct(v,ec).ymean = meanres.(mfield)(meanind,roiind);
        if hasp
            eccstruct(v,ec).p = meanres.ppara(conind{v,ec},roiind)';
            eccstruct(v,ec).pmean = meanres.ppara(meanind,roiind);
            eccstruct(v,ec).pslope = meanres.ppara(customconind{v,ec},roiind);
        end
        eccstruct(v,ec).ysingles = squeeze(groupres.r(conind{v,ec},...
            roiind,:))';
        % need to hack in the 4th point (delta dir 0) - this is quite
        % messy (only for block off diagonals)
        if strcmp(eccstruct(v,ec).view,'view_across_') && ec<4
            extraind = strfindcell(meanres.rows_contrast,sprintf(...
                '%s%s_direction_11_',vstruct(v).prefix,...
                eccstruct(v,ec).group));
            assert(numel(extraind)==1);
            name = meanres.rows_contrast{extraind};
            ind = strfind(name,'_m');
            ind = ind(end);
            eccstruct(v,ec).x(end+1) = str2double(name(ind+2:end));
            eccstruct(v,ec).y(end+1) = meanres.(mfield)(extraind,roiind);
            eccstruct(v,ec).p(end+1) = meanres.(mfield)(extraind,roiind);
            eccstruct(v,ec).ysingles(:,end+1) = ...
                squeeze(groupres.r(extraind,roiind,:));
        end
    end
end

vstyles = struct('marker',{par.markertype,'s'},'line',{'-','--'},...
    'facecolor',{par.basecolors,ones(size(par.basecolors))},...
    'xoff',{par.viewxoff(1),par.viewxoff(2)},...
    'textcolor',{ones(size(par.basecolors)),par.basecolors});
if numel(vstruct)==1
    [vstyles.xoff] = deal(0);
end
nsub = numel(groupres.z_subject);

subx = linspace(par.subxscale(1),par.subxscale(2),nsub)' * 3;

handles.points = NaN([numel(vstruct) necc 40]);
handles.meanpoints = NaN([numel(vstruct),necc,4]);
handles.textmarker = NaN([numel(vstruct),necc,4]);

% and now we can just plot each group
for v = 1:numel(vstruct)
    for ec = 1:necc
        thisecc = eccstruct(v,ec);
        thissubx = bsxfun(@plus,thisecc.x,subx);
        thisn = numel(thissubx);
        thisc = size(thissubx,2);
        handles.points(v,ec,1:thisn) = plotarray(thissubx(:),...
            thisecc.ysingles(:),vstyles(v).marker,'markersize',3,...
            'markerfacecolor',thisecc.color,'markeredgecolor',...
            thisecc.color);
        % DEBUG - where to the zeros in handles.meanpoints come from?
        handles.meanpoints(v,ec,1:thisc) = plotarray(thisecc.x(:),...
            thisecc.y(:),vstyles(v).marker,'markersize',8,...
            'markerfacecolor',[1 1 1],...
            'markeredgecolor',thisecc.color);
        handles.textmarker(v,ec,1:thisc) = textscatter([thisecc.x(:),thisecc.y(:)],...
            repmat({num2str(thisecc.delta,'%d')},numel(thisecc.x),1));
        %handles.gslope(v,ec) = lslinebetter(handles.meanpoints(v,ec,:),...
            %true,'adaptive','linestyle',viewstyles(v).line,...
            %'color',thisecc.color * .5,...
            %'linewidth',par.fitlinewidth);
        uistackn(handles.meanpoints(v,ec,:),'top');
        uistackn(handles.textmarker(v,ec,:),'top');
        %handles.tp(v,ec) = text(thisecc.xmean,max(ylim),...
            %sprintf('%s %s',thisecc.group,...
            %char(p2str(thisecc.pmean,3,true))),...
            %'rotation',45,'color',thisecc.color);
        %handles.pline(v,ec) = line(thisecc.xmean * [1 1]',...
            %[mean(thisecc.y(:)); 100],'color',thisecc.color,...
            %'linestyle','-','linewidth',1);
        %set(handles.pline(v,ec),'tag','getdatalims=0','clipping','on');
        %if thisecc.pmean < .05
            %set(handles.tp(v,ec),'fontangle','normal','fontweight','bold');
        %else
            %set(handles.tp(v,ec),'fontangle','italic','fontweight',...
                %'normal');
        %end
    end % ec = 1:necc
    % fit slope to ecc(0) cases (no, now testing all again)
    handles.slope(v) = lslinebetter(handles.meanpoints(v,:,:),...
        true,[],'color',par.fitlinecolor,'linewidth',par.fitlinewidth);
end % v = 1:numel(vstruct)


% styling
uistackn(handles.points,'bottom');
% uistackn(handles.gslope,'bottom');
uistackn(handles.slope,'bottom');
ylabel(par.ylabel);
xlabel(facedist_names('full'));
set(gca,'position',par.axscale * [1 1 1 1],'tickdir','out',...
    'ticklength',par.ticklength);
viewlabel = vstruct.xlabel2;
set(handles.f,'name',['eccdir' viewlabel ' ' roi]);

xl = [0 50];
xlim(xl);
textx = xl(2) + range(xl) * .05;
if chancelevel ~= 50
    handles.chance = line(xl,[chancelevel chancelevel],...
        'linestyle',par.chancelinestyle,...
        'linewidth',par.chancelinewidth,...
        'color',par.chancelinecolor);
    uistack(handles.chance,'bottom');
    handles.tchance = text(textx,chancelevel,'chance',...
        'horizontalalignment','left',...
        'verticalalignment','middle');
end

% Now we concatenate all the data we're about to plot and figure out the
% maximal statistic p threshold.
%if ~par.dochanceshade || ~isfield(meanres,'ppara') || all(isnan(meanres.ppara(:)))
if ~isfield(meanres,'ppara') || all(isnan(meanres.ppara(:))) || chancelevel==50
    [handles.chanceshade,handles.chancetext] = deal(NaN);
else
    ally = cat(2,eccstruct.ysingles);
    allx = bsxfun(@plus,cat(2,eccstruct.x),subx);
    [handles.chanceshade,handles.tchancethresh] = ...
        facedist_chanceshadewrap(allx,ally,max(xlim) + range(xlim)*.05);
end
