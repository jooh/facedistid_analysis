% facedist_rdmplots(meanres,groupres,stimuli,figdir,ts)
function facedist_rdmplots(meanres,groupres,stimuli,figdir,ts)

nroi = numel(meanres.cols_roi);

ncon = npairs2n(size(meanres.mean,1));

switch ncon
    case 12
        rowlabs = stimuli(1:12);
        iscrossdecode = any(strfind(figdir,'crossdecode'));
        if iscrossdecode
            collabs = stimuli(13:24);
        else
            collabs = stimuli(1:12);
        end
        nrow = 2;
    case 24
        rowlabs = stimuli;
        collabs = stimuli;
        nrow = 4;
    otherwise
        error('unknown ncon: %d',ncon);
end

figdir = fullfile(figdir,'facedist');
mkdirifneeded(figdir);

F = figure;
par = facedist_plotpar;
for r = 1:nroi
    ax = gca;
    [~,intmap,cmap] = rdmplot(ax,meanres.mean(:,r),'gridlines',...
        par.gridlines,'gridcolor',[1 1 1],'limits','symmetrical');
    % plot without cb label
    printstandard(fullfile(figdir,sprintf('rdm_%s',...
        meanres.cols_roi{r})),'formats=eps');
    newpos = get(ax,'position')*.75;
    set(ax,'position',newpos);
    cbax = axes('position',...
        [newpos(1:2) + [newpos(3) 0] .2 .2]);
    colorbarbetter(cbax,intmap,cmap,'orientation','horizontal',...
        'label',par.ylabel{1},'scale',.5,'tick',[]);
    % round the tick labels a bit 
    newlim = roundplotlimits(cbax,'x',[],2);
    set(cbax,'xlim',newlim,'xticklabelmode','auto','xtick',...
        sort([newlim 0]));
    % just kill the main RDM since we never use it anyway
    delete(ax);
    printstandard(fullfile(figdir,sprintf('rdm_cb_%s',...
        meanres.cols_roi{r})),'formats=eps');
    clf(F);
    % plot with labels
    ax = gca;
    [~,intmap,cmap] = rdmplot(ax,meanres.mean(:,r),'gridlines',...
        par.gridlines,'gridcolor',[1 1 1],'limits','symmetrical',...
        'xlabels',{collabs.image},'ylabels',{rowlabs.image},...
        'nrows',nrow);
    printstandard(fullfile(figdir,sprintf('rdm_wlabels_%s',...
        meanres.cols_roi{r})),'formats=eps');
    clf(F);
end
close(F);
