% aamod_pilab_rsa_predictorvisualisation_group and
% aamod_pilab_rdms_visualisation_group.
%
% mds plots with supporting lines using plotmdslines so we have a set of
% colours for eccentricity (purple to pink). So we connect the points
% according to the order specified by direction and colour them according
% to eccentricity. It will hopefully be evident that a) the different
% eccentricity bands are distinct, b) the more eccentric bands are further
% apart, c) the grid folds in on itself for all regions EXCEPT EVC and
% behaviour.
%
% function facedist_rdmplots_view(rdms,stimuli,figdir,ts)
%
% OR
%
% function facedist_rdmplots_view(meanres,groupres,stimuli,figdir,ts)
function facedist_rdmplots_view(varargin)

% handle different calling syntaxes
switch nargin
    case 4
        % model RDMs
        [rdms,stimuli,figdir,ts] = deal(varargin{:});
        names = {rdms.name};
        rdms = asrdmvec(rdms);
    case 5
        % data RDMs
        [meanres,groupres,stimuli,figdir,ts] = deal(varargin{:});
        names = meanres.cols_roi;
        rdms = meanres.mean;
    otherwise
        error('unsupported nargin: %d',nargin);
end

ncon = npairs2n(size(rdms,1));
nroi = size(rdms,2);

basemask = logical(1-eye(ncon));

figdir = fullfile(figdir,'facedist');
mkdirifneeded(figdir);

switch ncon
    case 12
        rowlabels = stimuli(1:12);
        % A bit hacky, but here we go...
        iscrossdecode = any(strfind(figdir,'crossdecode'));
        if iscrossdecode
            collabels = stimuli(13:24);
            fn = 'rdm_across_';
        else
            collabels = stimuli(1:12);
            fn = 'rdm_within_';
        end
        % RDM
        F = figure;
        set(F,'renderer','painters');
        for r = 1:nroi
            thisrdm = asrdmmat(rdms(:,r));
            rdmplot(gca,thisrdm,'ylabels',{rowlabels.image},...
                'xlabels',{collabels.image},...
                'cmap',ts.cmap,'nrows',2,'gridlines',...
                ts.gridlines,'gridcolor',ts.gridcolor,...
                'doranktrans',false,'limits','symmetrical');
            printstandard(fullfile(figdir,[fn ...
                stripbadcharacters(names{r})]));
            clf(F);
        end
        close(F);
        % MDS
        F = figurebetter([],[12 12]);
        for r = 1:nroi
            thisrdvec = rdms(:,r);
            try
                % with faces
                facedist_gridmds(F,thisrdvec,rowlabels);
                printstandard(fullfile(figdir,[fn 'plotmdslines_' ...
                    stripbadcharacters(names{r})]),'formats',{'png'},...
                    'r=1200');
                clf(F);
                facedist_gridmds(F,thisrdvec,[]);
                printstandard(fullfile(figdir,[fn 'plotmdslinesnoface_' ...
                    stripbadcharacters(names{r})]),'formats=eps');
                clf(F);
            catch err
                if strcmp(err.identifier,'stats:mdscale:ColocatedPoints')
                    logstr('full mdscale failed for %s, skipping\n',...
                        names{r});
                    continue
                end
                rethrow(err)
            end
        end
        close(F);
    case 24
        % set up masks
        withinmask = basemask;
        withinmask(13:24,1:12) = false;
        withinmask(1:12,13:24) = false;
        acrossmask = basemask;
        acrossmask(1:12,1:12) = false;
        acrossmask(13:24,13:24) = false;
        acrossmask = zerodiagonal(acrossmask);

        % MDS 
        F = figurebetter([12 12]);
        % full RDM
        fn = 'rdm_full_';
        for r = 1:nroi
            thisrdm = asrdmmat(rdms(:,r));
            if any(isnan(thisrdm))
                continue
            end
            if isequal(thisrdm(1:12,1:12),thisrdm(1:12,13:24))
                % will use 'collapsed' RDM below instead.
                continue
            end
            try
                % with faces
                facedist_gridmds(F,thisrdm,stimuli);
                printstandard(fullfile(figdir,[fn 'plotmdslines_' ...
                    stripbadcharacters(names{r})]),'formats=png','r=1200');
                clf(F);
                % no faces
                facedist_gridmds(F,thisrdm,[]);
                printstandard(fullfile(figdir,[fn 'plotmdslinesnoface_' ...
                    stripbadcharacters(names{r})]),'formats=eps');
                clf(F);
            catch err
                if strcmp(err.identifier,'stats:mdscale:ColocatedPoints')
                    logstr('full mdscale failed for %s, skipping\n',...
                        names{r});
                    continue
                end
                error('non-mdscale error');
            end
        end

        % within view
        fn = 'rdm_withinmean_';
        for r = 1:nroi
            thisrdm = asrdmmat(rdms(:,r));
            % within view - averaged
            thisrdm = matmean(thisrdm(1:12,1:12),thisrdm(13:24,13:24));
            rdv = asrdmvec(thisrdm);
            if any(isnan(rdv)) || all(rdv==rdv(1))
                continue
            end
            try
                % with faces
                facedist_gridmds(F,thisrdm,stimuli(1:12));
                printstandard(fullfile(figdir,[fn 'plotmdslines_' ...
                    stripbadcharacters(names{r})]),'formats',{'png'},...
                    'r=1200');
                clf(F);
                % no faces
                facedist_gridmds(F,thisrdm,[]);
                printstandard(fullfile(figdir,[fn 'plotmdslinesnoface_' ...
                    stripbadcharacters(names{r})]),'formats=eps');
                clf(F);
            catch err
                if strcmp(err.identifier,'stats:mdscale:ColocatedPoints')
                    logstr('within view mdscale failed for %s, skipping\n',...
                        names{r});
                    continue
                end
                error('non-mdscale error');
            end
        end
        close(F);

        % RDMs
        F = figure;
        set(F,'renderer','painters');
        % within view
        for r = 1:nroi
            thisrdm = asrdmmat(rdms(:,r));
            % within view - averaged
            withinrdm = matmean(thisrdm(1:12,1:12),thisrdm(13:24,13:24));
            if all(isnan(asrdmvec(withinrdm)))
                continue
            end
            rdmplot(gca,withinrdm,'labels',{stimuli(1:12).image},...
                'cmap',ts.cmap,'nrows',2,'gridlines',...
                ts.gridlines,'gridcolor',ts.gridcolor,...
                'doranktrans',false,'limits','symmetrical');
            printstandard(fullfile(figdir,['rdm_withinmean_' ...
                stripbadcharacters(names{r})]));
            clf(F);
            % within view - full
            withinrdm = thisrdm;
            withinrdm(acrossmask) = NaN;
            rdmplot(gca,withinrdm,'labels',{stimuli.image},...
                'cmap',ts.cmap,'nrows',4,'gridlines',...
                ts.gridlines,'gridcolor',ts.gridcolor,...
                'doranktrans',false,'limits','symmetrical');
            printstandard(fullfile(figdir,['rdm_withinfull_' ...
                stripbadcharacters(names{r})]));
            clf(F);
        end

        % across view
        for r = 1:nroi
            thisrdm = asrdmmat(rdms(:,r));
            acrossrdm = thisrdm;
            acrossrdm(withinmask) = NaN;
            if all(isnan(asrdmvec(acrossrdm)))
                continue
            end
            assert(issplitdatardm(acrossrdm));
            rdmplot(gca,acrossrdm,'labels',{stimuli.image},...
                'cmap',ts.cmap,'nrows',2,'gridlines',...
                ts.gridlines,'gridcolor',ts.gridcolor,...
                'doranktrans',false,'collapsesplitrdm',true,...
                'limits','symmetrical');
            printstandard(fullfile(figdir,['rdm_across_' ...
                stripbadcharacters(names{r})]));
            clf(F);
        end
        close(F);

    otherwise
        error('unsupported RDM size: %d',ncon);
end
