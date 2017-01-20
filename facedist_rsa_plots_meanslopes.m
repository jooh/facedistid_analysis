% facedist_rsa_plots_meanslopes(figdir,res,groupres,predictors,varargin)
function facedist_rsa_plots_meanslopes(figdir,res,groupres,predictors,varargin)

nroi = numel(res.cols_roi);
mkdirifneeded(figdir);

glmmode = false;
if ~any(strcmp(res.rows_contrast,'face01'))
    % get the customfits (gives us x values)
    customfits = facedist_customfits_mean(res.rows_contrast);
    suffind = [1 2];
else
    % this is a bit hacky but we can actually repurpose the contrasts
    % struct instead
    nface = 12;
    if any(strcmp(res.rows_contrast,'face24'))
        nface = 24;
    end
    customfits = facedist_glm_contrasts(res.rows_contrast(1:nface));
    suffind = [1];
    glmmode = true;
end

% and the rows in the result file to which they belong
[~,custominds] = intersect(...
    groupres.rows_contrast,{customfits.name},'stable');
custominds = num2cell(custominds);
[customfits.ind] = deal(custominds{:});

if glmmode
    if nface == 24
        viewgroups = struct('prefix',{'view_left_','view_right_','all_',...
            'ALL'},...
            'xlabel2',{'(view left)','(view right)',' ','ALL'});
    else
        viewgroups = struct('prefix','all_','xlabel2','');
    end
else
    withinind = strfindcell(groupres.rows_contrast,'view_within_');
    if isempty(withinind)
        viewgroups = struct('prefix','',...
            'xlabel2','');
    else
        % 3rd entry is a placeholder
        viewgroups = struct('prefix',{'view_within_','view_across_','ALL'},...
            'xlabel2',{'(within viewpoint)','(across viewpoint)','ALL'});
    end
end

% detect facepairs mode
chancelevel = 0;
if any(strfindcell(res.cols_roi,'facepairs')) && ~any(strfindcell(res.cols_roi,'EVC'))
    logstr('detected facepairs mode\n');
    chancelevel = 50;
end

% index RDM for MeanRSA legend
if ~glmmode
    names = {'sub-sub','typical-typical',...
        'caricature-caricature','sub-typical','typical-caricature',...
        'sub-caricature'};
    inds = cellfun(@(x)strfindcell({predictors.name},...
        [viewgroups(1).prefix x '_m']),names,...
        'uniformoutput',1);
    thisrdmat = single(asrdmmat(predictors(inds))~=0);
    nmat = size(thisrdmat,3);
    % add offset for each
    urdmat = bsxfun(@times,thisrdmat,reshape(1:nmat,[1 1 nmat]));
    % and sum to obtain final index matrix
    fullmat = sum(urdmat,3);
    colormap = [1 1 1; facedist_colors(names{:})];
    im = asrdmimage(fullmat,colormap);
    % save as image
    fh = figure;
    par = facedist_plotpar;
    imageplot(gca,im,'gridlines',par.gridlines,'gridcolor',[1 1 1]);
    printstandard(fullfile(figdir,'rdm_legend_meanrsa'),'formats=eps');
    close(fh);
end

for v = 1:numel(viewgroups)
    if v>1 && v==numel(viewgroups)
        thisview = viewgroups(1:2);
    else
        thisview = viewgroups(v);
    end
    for r = 1:nroi
        thisroi = groupres.cols_roi{r};
        handles.ecc{v,r} = facedist_abseccplot(groupres,res,...
            customfits,thisview,thisroi,chancelevel);
        handles.dir{v,r} = facedist_dirplot(groupres,res,...
            customfits,thisview,thisroi,chancelevel);
        % new
        if ~glmmode && (numel(viewgroups)==1 || v < numel(viewgroups))
            handles.eccdir{v,r} = facedist_eccdirplot(groupres,res,...
                customfits,thisview,thisroi,chancelevel);
        end
        % initial, crude update of y pos
        settexty(handles.dir{v,r}.tp,max(ylim(handles.dir{v,r}.ax)));
        settexty(handles.ecc{v,r}.tp,max(ylim(handles.ecc{v,r}.ax)));
        % settexty(handles.eccdir{v,r}.tp,max(ylim(handles.eccdir{v,r}.ax)));
    end
end
handles.ecc = cell2mat(handles.ecc);
handles.dir = cell2mat(handles.dir);

% ECCDIR
if ~glmmode
    handles.eccdir = cell2mat(handles.eccdir);
    % both
    outdir = fullfile(figdir,'eccdir_both');
    mkdirifneeded(outdir);
    printbyname([handles.eccdir.f],outdir,'formats=eps');

    % means only
    supersetter(handles.eccdir,{'points'},false,'visible','off');
    updateploty(handles.eccdir,chancelevel);
    outdir = fullfile(figdir,'eccdir_means');
    mkdirifneeded(outdir);
    printbyname([handles.eccdir.f],outdir,'formats=eps');

    % singles only
    supersetter(handles.eccdir,{'points'},false,'visible','on');
    supersetter(handles.eccdir,{'meanpoints'},false,'visible','off');
    supersetter(handles.eccdir,{'textmarker'},false,'visible','off');
    updateploty(handles.eccdir,chancelevel);
    outdir = fullfile(figdir,'eccdir_singles');
    mkdirifneeded(outdir);
    printbyname([handles.eccdir.f],outdir,'formats=eps');

    close([handles.eccdir.f]);
end

% normal mode
gdir = figdir;

% DIR
for suff = suffind
    if suff == 1
        % just carry on
        suff = '';
    else
        % change xlim, make ptext etc invisible
        suff = 'x3';
        set([handles.dir.ax],'xlim',[1.5 4.5]);
        % need to get rid of the tp and mtp, if necessary
        for s = 1:size(handles.dir,1)
            for h = 1:size(handles.dir,2)
                deleten(handles.dir(s,h).tp(:,1,:));
                handles.dir(s,h).tp(:,1,:) = NaN;
                deleten(handles.dir(s,h).mtp(1,1,:));
                handles.dir(s,h).mtp(1,1,:) = NaN;
            end
        end
    end

    mfields = {'mpoints','mslope','mtr','mtp','mmeanpoints','mmeanerrors',...
        'mshade','mchanceshade','mchancetext'};

    % first grouped plots - turn off all means
    supersetter(handles.dir,mfields,false,'visible','off');
    % and turn on all non-means
    supersetter(handles.dir,mfields,true,'visible','on');

    % grouped singles 
    supersetter(handles.dir,{'meanpoints','meanerrors'},false,'visible',...
        'off');
    % adjust y lim and p text
    updateploty(handles.dir,chancelevel);
    outdir = fullfile(gdir,['dir_gsingles' suff]);
    mkdirifneeded(outdir);
    printbyname([handles.dir.f],outdir,'formats=eps');

    % grouped means
    supersetter(handles.dir,{'meanpoints','meanerrors'},false,'visible',...
        'on');
    supersetter(handles.dir,{'points'},false,'visible','off');
    updateploty(handles.dir,chancelevel);
    outdir = fullfile(gdir,['dir_gmeans' suff]);
    mkdirifneeded(outdir);
    printbyname([handles.dir.f],outdir,'formats=eps');

    % now mean plots - turn off everything else
    supersetter(handles.dir,mfields,true,'visible','off');
    % and turn the mean fields back on
    supersetter(handles.dir,mfields,false,'visible','on');

    % mean singles 
    supersetter(handles.dir,{'mmeanpoints','mmeanerrors'},false,...
        'visible','off');
    % adjust y lim and p text
    updateploty(handles.dir,chancelevel);
    outdir = fullfile(gdir,['dir_msingles' suff]);
    mkdirifneeded(outdir);
    printbyname([handles.dir.f],outdir,'formats=eps');

    % mean means
    supersetter(handles.dir,{'mmeanpoints','mmeanerrors'},false,...
        'visible','on');
    supersetter(handles.dir,{'mpoints'},false,'visible','off');
    updateploty(handles.dir,chancelevel);
    outdir = fullfile(gdir,['dir_mmeans' suff]);
    mkdirifneeded(outdir);
    printbyname([handles.dir.f],outdir,'formats=eps');
end
% finished with dir
close(handles.dir.f);

% ECC
% singles only
supersetter(handles.ecc,{'meanpoints','meanerrors'},false,'visible','off');
% adjust y lim and p text
updateploty(handles.ecc,chancelevel);
outdir = fullfile(figdir,'ecc_singles');
mkdirifneeded(outdir);
printbyname([handles.ecc.f],outdir,'formats=eps');

% means only
supersetter(handles.ecc,{'meanpoints','meanerrors'},false,'visible','on');
supersetter(handles.ecc,{'points'},false,'visible','off');
% adjust y lim and p text
updateploty(handles.ecc,chancelevel);
outdir = fullfile(figdir,'ecc_means');
mkdirifneeded(outdir);
printbyname([handles.ecc.f],outdir,'formats=eps');
% finished with ecc
close([handles.ecc.f]);

% Export stats
if ~glmmode
    res.cols_roi = facedist_names(res.cols_roi{:});
    res.rows_contrast = facedist_names(res.rows_contrast{:})';
    roinames = {'early visual cortex',...
        'occipital face area','transverse occipital sulcus',...
        'fusiform face area','parahippocampal place area'};
    if numel(viewgroups)>1
        connames = {'view within: sub-caricature',...
            'view across: sub-caricature',...
        'view within: typical','view across: typical',...
        'view within: caricature','view across: caricature',...
        'view within: eccentricity slope',...
        'view across: eccentricity slope'};
    else
        connames = {'sub-caricature','typical','caricature',...
            'eccentricity slope'};
    end
    roidata2csv(res,fullfile(figdir,'meanrsa_stats.csv'),...
        connames,roinames,3);
end

function settexty(th,y)
settextpos(th,y,'y');

function updateploty(h,chancelevel)
% allow the y lim to float
set([h.ax],'ylimmode','auto','ytickmode','auto','yticklabelmode','auto');

switch chancelevel
    case 0
        p = 1;
        special = chancelevel;
    case 50
        p = 0;
        special = 0;
    otherwise
        error('unknown chance level: %f',chancelevel);
end

% match up the limits pairwise
for r = 1:size(h,2)
    % first find the limits
    for n = 1:size(h,1)
        lims = getdatalims(h(n,r).ax,'y');
        if chancelevel==50 && lims(2)>100
            lims(2) = 100;
        end
        ylim(h(n,r).ax,lims);
        lims = roundplotlimits(h(n,r).ax,'y',special,p);
    end
    % match the limits
    newax = matchaxis([h(:,r).ax]);
    % and update y positions of text labels
    % pad
    newy = newax(4) + range(newax(3:4))*.1;
    if isfield(h,'tp')
        settexty(justcat(h(:,r).tp),newy);
    end
    if isfield(h,'mtp')
        settexty(justcat(h(:,r).mtp),newy);
    end
    % now always set ticks to make sure formatting is consistent across
    % ROIs
    minimalticks([h(:,r).ax],[],special);
    t = get(h(1,r).ax,'ytick');
    %tstr = mat2strcell(t,precision);
    arrayfun(@(thisax)set(thisax,'ytick',t),justcat(h(:,r).ax));
end

function supersetter(st,mfields,mflag,varargin)
if mflag
    mfields = cellfun(@(thism)thism(2:end),mfields,'uniformoutput',0);
end

cellfun(@(thisfield)cellfun(@(thisentry)...
    setn(thisentry,varargin{:}),...
    {st.(thisfield)}),mfields);
