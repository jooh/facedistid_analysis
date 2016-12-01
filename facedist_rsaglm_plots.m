% visualise RSAGLM results.
%
% NB varargin are ignored and only retained for support of calls from
% aamod_pilab_rsa_visualisation_rfx.
%
% facedist_rsaglm_plots(figdir,res,groupres,predictors,varargin)
function facedist_rsaglm_plots(figdir,res,groupres,predictors,varargin)

par = facedist_plotpar;
mkdirifneeded(figdir);

conorder = {'delta_eccentricity','delta_direction'};
ncon = numel(conorder);
conwithin = prefix(conorder,'view_within_');
conacross = prefix(conorder,'view_across_');
[~,~,indswithin] = intersect(conwithin,res.rows_contrast,'stable');
[~,~,indsacross] = intersect(conacross,res.rows_contrast,'stable');
% set up factorial structure for ANOVA
facts = struct('name','metric','levs',[1:ncon]);
fullrdm = false;
if isempty(indswithin) && isempty(indsacross)
    % collapsed view fit
    [~,~,inds] = intersect(conorder,res.rows_contrast,'stable');
    assert(isequal(numel(inds),ncon),...
        'did not find the correct number of conditions in rows_contrast');
else
    fullrdm = true;
    assert(isequal(numel(indswithin),numel(indsacross),ncon),...
        'did not find the correct number of conditions in rows_contrast');
    % add second view
    facts.levs = [ones(1,ncon) ones(1,ncon)*2];
    facts(2) = struct('name','view','levs',[1:ncon 1:ncon]);
    conorder = cellfun(@(x){['view_within_' x],['view_across_' x]},...
        conorder,'uniformoutput',0);
end

res.cols_roi = facedist_names(res.cols_roi{:});
groupres.cols_roi = facedist_names(groupres.cols_roi{:});
res.facecolor = repmat({facedist_colors('rsaglm')},[numel(res.rows_contrast),1]);
res.edgecolor = repmat({facedist_colors('rsaglm')},[numel(res.rows_contrast),1]);
% outline the across view bars if they are present
[res.facecolor{strfindcell(res.rows_contrast,'view_across')}] = deal(...
    [1,1,1]);
res.facecolor = vertcat(res.facecolor{:});
res.edgecolor = vertcat(res.edgecolor{:});
% version with inference
[~,bardir] = mkdirifneeded(fullfile(figdir,'bars_inference'));
[handles,xdata,ydata,zdata] = roidata2figure(res,groupres,...
    'conind',conorder,...
    'xticklabels',facedist_names('delta_eccentricity','delta_direction'),...
    'precision',2,'specialval',0,'ylab',par.ylabel_rsaglm,...
    'xscalefactor',2.5,'grouptarget','b','restarget','mean');
printbyname([handles.figure],bardir);
close([handles.figure]);
% version without inference
[~,bardiralt] = mkdirifneeded(fullfile(figdir,'bars_noinference'));
[handles,xdata,ydata,zdata] = roidata2figure(res,groupres,...
    'conind',conorder,...
    'xticklabels',facedist_names('delta_eccentricity','delta_direction'),...
    'precision',2,'specialval',0,'ylab',par.ylabel_rsaglm,...
    'xscalefactor',2.5,'grouptarget','b','restarget','mean',...
    'dozerotest',false,'docontrasts',false);
printbyname([handles.figure],bardiralt);
close([handles.figure]);

% cook up ANOVA - now by handles.singles I think
% hm. Hard to work with. rather than add even more bells and whistles to
% roidata2figure I think we might just return a data struct. Another idea
% would be to just make all the figures and return so that the user can add
% additional stuff post flight, or just printbyname
nroi = numel(res.cols_roi);
roidata = cell(1,nroi);
for r = 1:nroi
    roidata{r} = vertcat(ydata(r).singles{:})';
    [~,stats(r).table,stats(r).stats] = anovan_rmwrap(roidata{r},facts);
    % export it while we're at it
    exporttable(stats(r).table,fullfile(figdir,sprintf('anova_%s.csv',...
        stripbadcharacters(res.cols_roi{r}))));
end
save(fullfile(figdir,'anovastats.mat'),'stats');

% ROI interaction
nsub = size(groupres.b,3);
if nsub>1
    % metric by ROI interaction
    [~,~,roiind] = intersect(facedist_names('EVC','FFA'),res.cols_roi);
    if numel(roiind)~=2
        logstr('could not find EVC and/or FFA, skipping roi by metric test\n');
    else
        facts(end+1) = struct('name','roi','levs',...
            ones(1,numel(facts(end).levs)));
        facts(2,:) = facts;
        facts(2,end).levs(:) = 2;
        % now roifacts is just
        roifacts = struct('name',{facts(1,:).name},'levs',...
            arrayfun(@(x)horzcat(facts(:,x).levs),1:size(facts,2),...
            'uniformoutput',0));
        [~,roistats.table,roistats.stats] = anovan_rmwrap(...
            horzcat(roidata{roiind}),roifacts);
        exporttable(roistats.table,fullfile(figdir,'roianova_EVC_x_FFA.csv'));
        save(fullfile(figdir,'anovastats_roi.mat'),'roistats');
    end
end

% RDM visualisation
meanrdms = struct('name',{'data','fitted','res'},...
    'RDM',cellfun(@asrdmmat,{nanmean(sqrtsigned(groupres.custom{3}),3),...
    nanmean(sqrtsigned(groupres.custom{1}),3),...
    nanmean(sqrtsigned(groupres.custom{2}),3)},'uniformoutput',0));
F = figure;
Fmds = figurebetter([],[12 12]);

for r = 1:nroi
    lims = [];
    for t = meanrdms
        figure(F);
        if isempty(lims)
            % set lims by first metric (data)
            lims = reduceprecision(max(abs(asrdmvec(t.RDM(:,:,r)))),1,...
                @ceil) .* [-1 1];
        end
        [~,intmap,cmap] = rdmplot(gca,t.RDM(:,:,r),'gridlines',...
            par.gridlines,'gridcolor',[1 1 1],'limits',lims);
        newpos = get(gca,'position')*.75;
        set(gca,'position',newpos);
        cbax = axes('position',[newpos(1:2) + [newpos(3) 0] .2 .2]);
        colorbarbetter(cbax,intmap,cmap,'orientation','horizontal',...
            'label',par.ylabel{1},'scale',.5,'tick',[]);
        printstandard(fullfile(figdir,sprintf('rdm_wcb_%s_%s',...
            stripbadcharacters(res.cols_roi{r}),t.name)),'formats=eps');
        clf(F);
        % MDS plot
        try
            figure(Fmds);
            mdsrdm = t.RDM(:,:,r);
            mdsrdm(mdsrdm<0) = 0;
            if all(mdsrdm(:)==0)
                % go to catch below
                error
            end
            facedist_gridmds(Fmds,mdsrdm,[]);
            printstandard(fullfile(figdir,sprintf(...
                'plotmdslinesnoface_%s_%s',stripbadcharacters(res.cols_roi{r}),t.name)),...
                'formats=eps');
            clf(Fmds);
        catch
            logstr('MDS failed.\n')
        end
    end
end
close(Fmds);
% also visualise the predictors with the same parameters
for p = 1:numel(predictors)
    figure(F);
    if any(strfind(predictors(p).name,'constant'))
        predictors(p).RDM(predictors(p).RDM>0) = 49;
    end
    % keep the same limits as the original
    [~,intmap,cmap] = rdmplot(gca,predictors(p).RDM,'gridlines',...
        par.gridlines,'gridcolor',[1 1 1],'limits',[-49 49]);
    printstandard(fullfile(figdir,['rdm_predictors_' ...
        predictors(p).name]),'formats=eps');
    clf(F);
end
close(F);

function exporttable(exptab,outpath)
exptab(strfindcell(exptab(:,1),'subject'),:) = [];
anovan_exportcsv(outpath,exptab);
