% Main function for running computational model fits for facedist project.
% If recompute is true we overwrite existing results. 
%
% GWP fit requires code from Kendrick Kay:
% https://github.com/kendrickkay/knkutils
%
% The face-space models require my comp model code:
% https://github.com/jooh/matlab-neuromodeling
%
% (which is included here under code/modeling, and gets added to the path if
% missing)
%
% NAMED INPUTS
% niter: 3125 - number of grid points to search over
% skipntune: true - we used to fit the number of model neurons, but no
% longer do so. This flag is deprecated.
% dirsuffix: '' - suffix to add to output dir (an underscore is prepended)
% collapseview: false - run full view fit, or collapsed view version
%
% EXAMPLE USAGE:
% % reproduce full view model fit in the ms (Figs 3:6)
% facedist_doit_models(false,'dirsuffix','_view_full','collapseview',false);
%
% % reproduce collapsed view model fit in the ms (S7 Fig)
% facedist_doit_models(false,'dirsuffix','_collapseview','collapseview',true);
%
% 2017-06-28 J Carlin
%
% facedist_doit_models(recompute,varargin)
function facedist_doit_models(recompute,varargin)

runstart = clock;
skipntune = true;
getArgs(varargin,{'niter',3125,'skipntune',true,'dirsuffix','',...
    'collapseview',false});

dirname = 'view';
nface = 24;
if collapseview
    dirname = 'collapseview';
    nface = 12;
end

%% PRELIMINARY SETTINGS
% coordinates for the test faces (scaled to standard deviations for the
% 2D space).
scalefactor = (sqrt(2)/sqrt(199));
xy = facedist_getrefxy(nface) * scalefactor;
xy_ecc = sqrt(sum(xy.^2,2));
% sub, typical, caricature for one direction only (in case we have two).
% the face-space models make identical predictions for view two
xyplotind = 1:12;
par = facedist_plotpar;
rois = {'EVC','FFA','OFA','TOS','PPA'};
nroi = numel(rois);
% with facepairs
nall = nroi + 1;
barinds = 1:nall;

%% WORK OUT DIRECTORY STRUCTURE
% 2 dirs up from CWD is the main BIDS directory
rootdir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
derivdir = fullfile(rootdir,'derivatives');
[~,resdir] = mkdirifneeded(fullfile(derivdir,sprintf('modeling_%s',dirsuffix)));

if isempty(which('ExemplarModel'))
    logstr('adding code/modeling to path\n');
    addpath(fullfile(rootdir,'code','modeling'));
end

%% LOAD DATA AND ESTIMATE REPRODUCIBILITY
% we depend on these files a lot for retrieving data RDMs - fMRI data from
% rsaglm_groupres and perceptual data from stimres
rsaglm_groupres = loadbetter(fullfile(dirmostrecent(fullfile(derivdir,...
    sprintf('aafacedist_%s_froi_rsaglm',...
    dirname),'aamod_pilab_rfx_*')),'pilab','result_group.mat'));
% slight hack to handle change in name handling between AA4 and 5
if iscell(rsaglm_groupres.z_subject{1})
    rsaglm_groupres.z_subject = cellfun(@(x)x{1}(1:6),...
        rsaglm_groupres.z_subject,'uniformoutput',0);
end
stimres = loadbetter(fullfile(derivdir,'facepairs','stimanalysis.mat'));
nsub = numel(rsaglm_groupres.z_subject);
ndgwp = 5;

% so here we go
reprofile = fullfile(resdir,'reprores.mat');
if ~exist(reprofile,'file') || recompute
    logstr('loading data, measuring GWP and estimating reproducibility\n')
    splits = logical(1-eye(nsub));
    [~,~,gind] = intersect(rois,rsaglm_groupres.cols_roi,'stable');
    % get the data and residuals
    for field = struct('name',{'data','residual'},'index',{3,2})
        data = arrayfun(...
            @(thisres)sqrtsigned(thisres.custom{field.index}(:,...
            strcmp(thisres.cols_roi,'facepairs'))),...
            stimres.rsaglm,'uniformoutput',0);
        data = zcat(data{:});
        if ~collapseview
            % NaN out to make full shape.
            dmat = zerodiagonal(NaN([nface nface nsub]));
            dmat(1:12,1:12,:) = asrdmmat(reshape(data,...
                [size(data,1) nsub]));
            data = reshape(asrdmvec(dmat),[nchoosek(nface,2),1,nsub]);
        end
        % then get the fROI data
        fdata = sqrtsigned(rsaglm_groupres.custom{field.index}(:,gind,:));
        % [dis x roi x sub] format
        repro.(field.name) = horzcat(data,fdata);
    end
    % now we can cook up a subres and plug in to roidata_rfx!
    for sub = 1:nsub
        substr = rsaglm_groupres.z_subject{sub};
        logstr('predictors for subject %d of %d...\n',sub,nsub);
        % import stimuli
        stimuli = facedist_importstimuli(substr);
        stimuli = stimuli(facedist_resortind(numel(stimuli)));
        [~,imdir] = mkdirifneeded(fullfile(resdir,['stim_' substr]));
        % build GWP predictors
        [repro.gwppredictors{sub},repro.gwpmeanresp{sub},...
            repro.gwprawresp{sub},repro.gwpleftind{sub}] = ...
            facedist_stim2gwppredictors(stimuli,imdir);
        if collapseview
            repro.gwppredictors{sub} = ...
                facedist_rsapredictors_collapseview(...
                repro.gwppredictors{sub});
            repro.gwpmeanresp{sub} = matmean(...
                repro.gwpmeanresp{sub}(1:12,:),...
                repro.gwpmeanresp{sub}(13:24,:));
            repro.gwprawresp{sub} = cellfun(@(thisbank)matmean(...
                thisbank(1:12,:),thisbank(13:24,:)),...
                repro.gwprawresp{sub},'uniformoutput',0);
        end
        rmat = NaN([3,nall]);
        train = splits(sub,:);
        % upper reproducibility estimate
        % [n by roi by sub] format
        rmat(1,:) = corrpairs(repro.data(:,:,sub),...
            nanmean(repro.data,3));
        % lower reproducibility estimate
        rmat(2,:) = corrpairs(repro.data(:,:,sub),...
            nanmean(repro.data(:,:,train),3));
        % lower reproducibility of residuals
        rmat(3,:) = corrpairs(...
            repro.residual(:,:,sub),...
            nanmean(repro.residual(:,:,train),3));
        if isnan(rmat(1,1))
            assert(~collapseview)
            nanind = ~isnan(repro.data(:,1,1));
            rmat(1,1) = corr(repro.data(nanind,1,sub),...
                nanmean(repro.data(nanind,1,:),3));
            rmat(2,1) = corrpairs(repro.data(nanind,1,sub),...
                nanmean(repro.data(nanind,1,train),3));
            rmat(3,1) = corrpairs(...
                repro.residual(nanind,1,sub),...
                nanmean(repro.residual(nanind,1,train),3));
        end
        % Fisher transform
        repro.subres(sub) = struct('rows_contrast',...
            {{'rep_upper','rep_lower','rep_lower_resid'}'},...
            'cols_roi',{[{'facepairs'},rois]},...
            'r',atanh(rmat),'tail',{repmat({'right'},[3 1])});
    end
    [repro.subres.name] = deal(rsaglm_groupres.z_subject{:});
    % (careful - the order of the ROIs here is ['facepairs',rois]
    % subres, but in rfx the order will be different)
    repro.cols_roi = [{'facepairs'},rois];
    save(reprofile,'repro');
else
    repro = loadbetter(reprofile);
end

%% CONFIGURE MODEL FITS
alldimlabels = struct('n','n units (log10)','w','tuning width (FWHM)',...
    's','unit distribution width (proportion)','e','tuning width (FWHM, radial)',...
    'r','tuning width (FWHM,tangential)','f','tuning saturation',...
    'a','population averaging (proportion)','o','horizontal offset',...
    'd','chi^2 degrees of freedom','v','unit distribution shape',...
    'aview','hemifield-specific population averaging (proportion)');
for d = 1:ndgwp
    fn{d} = sprintf('c%d',d);
    fnfull = sprintf('%d cycles',d);
    alldimlabels.(fn{d}) = fnfull;
end

outfile_modelst = fullfile(resdir,'facedist_modelst.mat');
if ~exist(outfile_modelst,'file') || recompute
    logstr('configuring models\n')
    % exemplar model, gaussian unit sample
    modelst(1).field = 'ex_gaussian';
    modelst(1).instance = 'ExemplarModel';
    modelst(1).name = 'exemplar tuning, gaussian distribution';
    % use norminv sampler to get the right spacing
    modelst(1).sampler = '@(x)facedist_norminvsample(x,1)';
    nd = 4;
    nex = ceil(nthroot(niter,nd));
    modelst(1).coord.n = linspace(2,5,nex);
    if skipntune
        % disable n tune
        nd = 3;
        nex = ceil(nthroot(niter,nd));
        modelst(1).coord.n = 3;
    end
    modelst(1).coord.s = linspace(.1,3,nex);
    modelst(1).coord.w = linspace(.2,10,nex);
    modelst(1).coord.a = linspace(0,1,nex);
    [modelst(1).grid.n,modelst(1).grid.s,modelst(1).grid.w,...
        modelst(1).grid.a] = ndgrid(...
        modelst(1).coord.n,modelst(1).coord.s,modelst(1).coord.w,...
        modelst(1).coord.a);
    modelst(1).n = numel(modelst(1).grid.n);
    modelst(1).sz = size(modelst(1).grid.n);
    modelst(1).nd = ndims(modelst(1).grid.n);
    modelst(1).skipntune = skipntune;
    modelst(1).skipaverage = false;

    % exemplar model, inverted gaussian sample
    modelst(end+1) = modelst(end);
    % force to outskirts of space to avoid fits where everything collapses
    % in (and prediction basically reduces to ex_gaussian)
    modelst(end).coord.s = linspace(1,3,nex);
    modelst(end).sampler = '@(x)facedist_norminvsample(x,-1)';
    modelst(end).field = 'ex_negativegaussian';
    modelst(end).name = 'exemplar tuning, negative gaussian distribution';

    % ramp model
    modelst(end+1).field = 'ramp';
    modelst(end).instance = 'RampModel';
    modelst(end).name = 'ramp tuning';
    modelst(end).sampler = 'circsample';
    nd = 4;
    nramp = ceil(nthroot(niter,nd));
    % p1 - n unit
    modelst(end).coord.n = linspace(2,5,nramp);
    if skipntune
        nd = 3;
        nramp = ceil(nthroot(niter,nd));
        modelst(end).coord.n = 3;
    end
    % p2/3 - saturation for sigmoid
    modelst(end).coord.f = linspace(0,3,nramp);
    % p1/2 - offset for sigmoid
    modelst(end).coord.o = linspace(0,5,nramp);
    % p3/4 - averaging proportion for population response
    modelst(end).coord.a = linspace(0,1,nramp);
    [modelst(end).grid.n,modelst(end).grid.f,modelst(end).grid.o...
        modelst(end).grid.a] = ndgrid(modelst(end).coord.n,...
        modelst(end).coord.f,modelst(end).coord.o,modelst(end).coord.a);
    modelst(end).n = numel(modelst(end).grid.n);
    modelst(end).sz = size(modelst(end).grid.n);
    modelst(end).nd = ndims(modelst(end).grid.n);
    modelst(end).skipntune = skipntune;
    modelst(end).skipaverage = false;

    % gwp
    modelst(end+1).field = 'gwpgrid';
    modelst(end).instance = 'GWPModel';
    modelst(end).name = 'gabor wavelet pyramid';
    modelst(end).sampler = NaN;
    % for now skip population averaging. 
    nlev = ceil(nthroot(niter,ndgwp));
    for d = 1:ndgwp
        fn{d} = sprintf('c%d',d);
        fnfull = sprintf('%d cycles',d);
        modelst(end).coord.(fn{d}) = logspace(0,4,nlev);
        coordcell{d} = modelst(end).coord.(fn{d});
    end
    % nb same level of averaging as elsewhere
    modelst(end).coord.a = linspace(0,1,nramp);
    coordcell{end+1} = modelst(end).coord.a; 
    fn{end+1} = 'a';
    if ~collapseview
        modelst(end).coord.aview = linspace(0,1,nramp);
        coordcell{end+1} = modelst(end).coord.aview; 
        fn{end+1} = 'aview';
    end
    % build grids
    [gridcell{1:numel(coordcell)}] = ndgrid(coordcell{:});
    % and assign
    % +~collapseview because we have an additional field-specific averaging
    % parameter if there are multiple views
    for d = 1:(ndgwp+1+~collapseview)
        modelst(end).grid.(fn{d}) = gridcell{d};
    end
    modelst(end).n = numel(modelst(end).grid.c1);
    modelst(end).sz = size(modelst(end).grid.c1);
    modelst(end).nd = ndims(modelst(end).grid.c1);
    % stop first dimension from getting dropped.
    modelst(end).skipntune = false;
    modelst(end).skipaverage = false;

    for n = 1:numel(modelst)
        % search over the same space, but fix averaging to 0
        if ~modelst(n).skipaverage
            modelst(end+1) = modelst(n);
            modelst(end).field = [modelst(end).field '_noav'];
            modelst(end).name = [modelst(end).name ' (a=0)'];
            modelst(end).coord.a = 0;
            modelst(end).grid = structfun(@(x)indexdim(x,1,ndims(x)),...
                modelst(end).grid,'uniformoutput',0);
            if isfield(modelst(end).grid,'aview')
                modelst(end).coord.aview = 0;
                modelst(end).grid = structfun(@(x)indexdim(x,1,ndims(x)),...
                    modelst(end).grid,'uniformoutput',0);
            end
            modelst(end).skipaverage = true;
            % have to re-do these since the dimensionality has changed
            firstfn = asrow(fieldnames(modelst(end).grid),1);
            modelst(end).sz = size(modelst(end).grid.(firstfn{1}));
            modelst(end).n = numel(modelst(end).grid.(firstfn{1}));
            modelst(end).nd = ndims(modelst(end).grid.(firstfn{1}));
        end
    end
    save(outfile_modelst,'modelst');
else
    modelst = loadbetter(outfile_modelst);
end

%% RUN MODEL FIT
outfile = fullfile(resdir,'facedist_modelfits.mat');
if ~exist(outfile,'file') || recompute
    logstr('running model fits.\n');
    datardv = nanmean(repro.data,3);
    for thisst = modelst
        fstr = thisst.field;
        logstr('fitting %s %s...\n',fstr,thisst.instance);
        tstart = clock;
        if ~any(strcmp(fstr,{'gwpgrid','gwpgrid_noav'}))
            thisrdv = cell(thisst.sz);
            for n = 1:thisst.n
                fitpar = structfun(@(thisv)thisv(n),thisst.grid);
                thisrdv{n} = fitit(fitpar,thisst.sampler,xy,...
                    thisst.instance);
            end
            % compare the mean RDMs to targets
            res.all(1).(fstr) = rungridfit(...
                datardv,thisst,thisrdv);
            for sub = 1:nsub
                res.singles(1).(fstr)(sub) = ...
                    rungridfit(repro.data(:,:,sub),thisst,thisrdv);
            end
            % leave-one-out - need to average single subject fitpars
            % and construct models from this for cross-validation.
            for sub = 1:nsub
                % first standard
                trainind = setdiff(1:nsub,sub);
                trainpar = matfun(@nanmean,...
                    res.singles(1).(fstr)(trainind).fitpar);
                for t = 1:nall
                    % model representation for training fit (mean)
                    res.group(1).(fstr)(sub).rdvtrain(:,t) = ...
                        mean(fitit(trainpar(:,t),thisst.sampler,xy,...
                        thisst.instance),2);
                    if any(isnan(res.singles(1).(fstr)(sub).fitpar(:,t)))
                        res.singles(1).(fstr)(sub).meanresponse(:,t) = NaN;
                        continue;
                    end
                    % model representation for test split. The second call
                    % to fitit here is just to get the mean responses out
                    % (rather than storing the mean response for a bunch of
                    % non-winning models).
                    [~,testmean] = fitit(...
                        res.singles(1).(fstr)(sub).fitpar(:,t),...
                        thisst.sampler,xy,...
                        thisst.instance);
                    res.singles(1).(fstr)(sub).meanresponse(:,t) = ...
                        mean(testmean,2);
                end
            end
        else
            % GWP fit
            for sub = 1:nsub
                thisrdv = cell(thisst.sz);
                thismeanresp = cell(thisst.sz);
                % weight the response vectors and recalculate the RDM
                % at this level
                X = repro.gwprawresp{sub};
                leftind = repro.gwpleftind{sub};
                parfor n = 1:thisst.n
                    fitpar = structfun(@(thisv)thisv(n),thisst.grid);
                    responses = arrayfun(@(b)X{b} * fitpar(b),...
                        1:numel(X),'uniformoutput',0);
                    % nface by nfilter
                    responses = cat(2,responses{:});
                    % if we are averaging
                    if numel(fitpar)>ndgwp
                        if numel(fitpar)==ndgwp+2
                            % translate toward population average -
                            % separately for filters left and right of
                            % midline.
                            responses(:,leftind) = morph(mean(...
                                responses(:,leftind),2),...
                                responses(:,leftind),fitpar(end-1));
                            responses(:,~leftind) = morph(mean(...
                                responses(:,~leftind),2),...
                                responses(:,~leftind),fitpar(end-1));
                        end
                        % and global population averaging
                        responses = morph(mean(responses,2),responses,...
                            fitpar(end));
                    end
                    % and build the RDM
                    thisrdv{n} = pdist(responses)';
                    thismeanresp{n} = mean(responses,2);
                end
                res.singles(1).(fstr)(sub) = rungridfit(...
                    repro.data(:,:,sub),thisst,thisrdv);
                allrdv{sub} = thisrdv;
                % now we just index in the relevant mean responses
                for t = 1:nall
                    if any(isnan(res.singles.(fstr)(sub).fitpar(:,t)))
                        res.singles.(fstr)(sub).meanresponse(:,t) = NaN;
                        continue
                    end
                    res.singles.(fstr)(sub).meanresponse(:,t) = ...
                        thismeanresp{res.singles.(fstr)(sub).winner(t)};
                end
            end
            % now we make res.all by just averaging the single subjects
            % so we need to average over the final dimension
            meanrdv = cell(thisst.sz);
            for n = 1:thisst.n
                % get the predictions from GWP
                subdata = cellfun(@(thissub)thissub{n},allrdv,...
                    'uniformoutput',0);
                meanrdv{n} = matfun(@nanmean,subdata{:});
            end
            res.all(1).(fstr) = rungridfit(datardv,thisst,meanrdv);
            % cross-validated estimates 
            % averaging RDMs (we can't generate an average model prediction
            % over subjects since every subject viewed different images and
            % consequently has a distinct model)
            for sub = 1:nsub
                trainind = setdiff(1:nsub,sub);
                res.group(1).(fstr)(sub).rdvtrain = matfun(@nanmean,...
                    res.singles(1).(fstr)(trainind).rdvwinmean);
            end
        end
        logstr('finished in %s\n',seconds2str(etime(clock,tstart)));
    end
    % save mat
    save(outfile,'res');
else
    res = loadbetter(outfile);
end

reprofile_withmodels = fullfile(resdir,'reprores_withmodels.mat');
if ~exist(reprofile_withmodels,'file') || recompute
    logstr('summarising model prediction performance\n')
    predicttargets = {'full','vid_corrdist'};
    for sub = 1:nsub
        % add the fixed model RDMs
        repro.fixedpredictors = loadbetter(fullfile(dirmostrecent(fullfile(...
            derivdir,'aa','aamod_pilab_importrsapredictors*')),...
            rsaglm_groupres.z_subject{sub},'pilab',...
            'pilab_rsapredictors.mat'));
        if collapseview
            repro.fixedpredictors = facedist_rsapredictors_collapseview(...
                repro.fixedpredictors);
        end
        [~,~,pinds] = intersect(predicttargets,...
            {repro.fixedpredictors.name},'stable');
        selectpredictors = repro.fixedpredictors(pinds);
        for thisp = selectpredictors
            pvec = asrdmvec(thisp);
            repro.subres(sub).r(end+1,:) = atanh(pearsonvec(pvec,...
                repro.data(:,:,sub)));
            % and special case for facepairs...
            nanind = ~isnan(repro.data(:,1,sub));
            if any(isnan(repro.data(:,1,sub)))
                repro.subres(sub).r(end,1) = atanh(...
                    pearsonvec(pvec(nanind),repro.data(nanind,1,sub)));
            end
            repro.subres(sub).rows_contrast{end+1} = thisp.name;
            repro.subres(sub).tail{end+1} = 'right';
        end
        % add the CV performance for each model
        for fn = fieldnames(res.group)'
            fnstr = fn{1};
            % Fisher transform
            repro.subres(sub).r(end+1,:) = atanh(corrpairs(...
                repro.data(:,:,sub),res.group.(fnstr)(sub).rdvtrain));
            % special handling of facepairs
            if any(isnan(repro.data(:,1,sub)))
                repro.subres(sub).r(end,1) = atanh(...
                    corrpairs(repro.data(nanind,1,sub),...
                    res.group.(fnstr)(sub).rdvtrain(nanind,1)));
            end
            repro.subres(sub).rows_contrast{end+1} = fnstr;
            repro.subres(sub).tail{end+1} = 'right';
        end
    end
    [repro.rfx,repro.group] = roidata_rfx(repro.subres,'nperm',1024,...
        'nboot',1000,'contrasts','roidata_allpairwisecontrasts');
    save(reprofile_withmodels,'repro');
else
    % NB we over-write the original variable
    repro = loadbetter(reprofile_withmodels);
end

% load up the GLM fMRI effects for later comparisons
groupres_uni_fmri = loadbetter(fullfile(dirmostrecent(fullfile(...
    derivdir,'aafacedist_view_froi_glm','aamod_pilab_rfx_*')),'pilab',...
    'result_group.mat'));
% Have to map the names in the fmri result to full length.
[~,~,coninds_fmri] = intersect({'all_mean_ecc_sub',...
    'all_mean_ecc_typical','all_mean_ecc_caricature'}',...
    groupres_uni_fmri.rows_contrast,'stable');

% run each model through the standard inference tools
[~,fitdir] = mkdirifneeded(fullfile(resdir,'modelfitrsa'));
fitfile = fullfile(fitdir,'fitted.mat');
if ~exist(fitfile,'file') || recompute
    logstr('running model fit through standard RSA inference\n')
    fitted.otherpredictors.meanrsa = facedist_rsapredictors_mean_tuning(...
        repro.fixedpredictors);
    fitted.otherpredictors.rsaglm = facedist_rsapredictors_glm(...
        repro.fixedpredictors);
    fitted.designvol = Volume(eye(nface),'metasamples',struct('chunks',...
        ones(nface,1)),'metafeatures',struct('labels',...
        {mat2strcell([1:nface],'face%02d')}),'frameperiod',1);

    for f = fieldnames(res.singles)'
        fstr = f{1};
        logstr('RSA inference for %s\n',fstr);
        figdir = fullfile(fitdir,fstr);
        mkdirifneeded(figdir);
        thissingle = res.singles.(fstr);

        for s = 1:numel(thissingle)
            valid = ~any(isnan(thissingle(s).rdvwinmean),1);
            disvol = Volume(thissingle(s).rdvwinmean(:,valid),...
                'metafeatures',...
                struct('names',{repro.cols_roi(valid)}));
            fitted.subres.(fstr).meanrsa(s) = roidata_rsa(disvol,...
                fitted.otherpredictors.meanrsa,...
                'rsaclass','MeanRSA','splitrsaclass','MeanRSA',...
                'customfits','facedist_customfits_mean',...
                'contrasts','facedist_rsa_rfx_contrasts_meanslopes');
            fitted.subres.(fstr).meanrsa(s).tail(strfindcell(...
                fitted.subres.(fstr).meanrsa(s).rows_contrast,...
                'mean_direction'))= {'none'};
            fitted.subres.(fstr).rsaglm(s) = roidata_glm(...
                fitted.otherpredictors.rsaglm,disvol,...
                'customfun',{'predictY','residuals','getdata'},...
                'glmclass','MultiRSA','nboot',0,...
                'contrasts','facedist_rsaglm_contrasts');
            datavol = Volume(thissingle(s).meanresponse(:,valid),...
                'metasamples',...
                struct('chunks',ones(nface,1)),'metafeatures',struct(...
                'names',{repro.cols_roi(valid)}));
            fitted.subres.(fstr).glm(s) = roidata_glm(...
                fitted.designvol,datavol,'contrasts',...
                'facedist_glm_contrasts');
        end
        [fitted.subres.(fstr).meanrsa.name] = deal(stimres.meanrsa.name);
        [fitted.subres.(fstr).rsaglm.name] = deal(stimres.meanrsa.name);
        % set all the tails to none to disable inference
        [fitted.subres.(fstr).meanrsa.tail] = deal(repmat({'none'},...
            size(fitted.subres.(fstr).meanrsa(1).tail)));
        [fitted.subres.(fstr).rsaglm.tail] = deal(repmat({'none'},...
            size(fitted.subres.(fstr).rsaglm(1).tail)));
        [fitted.subres.(fstr).glm.name] = deal(stimres.meanrsa.name);
        [fitted.subres.(fstr).glm.tail] = deal(repmat({'none'},...
            size(fitted.subres.(fstr).glm(1).tail)));

        % RFX models and figures
        % RSAGLM
        [fittedrfx.rsaglm,grouprfx.rsaglm] = roidata_rfx(...
            fitted.subres.(fstr).rsaglm,'targetfield','b',...
            'nboot',0,'nperm',0,'assumeregister',0);
        glmdir = fullfile(figdir,'rsaglm');
        mkdirifneeded(glmdir);
        glmargs = {glmdir,fittedrfx.rsaglm,grouprfx.rsaglm,...
            fitted.otherpredictors.rsaglm,'groupmtarget','b','errtarget',...
            'sterr','ptarget',[],'mtarget','mean'};
        facedist_rsaglm_plots(glmargs{:});

        % Mean
        [fittedrfx.meanrsa,grouprfx.meanrsa] = roidata_rfx(...
            fitted.subres.(fstr).meanrsa,'targetfield','r',...
            'nboot',0,'nperm',0,'assumeregister',0);
        meandir = fullfile(figdir,'meanrsa');
        mkdirifneeded(meandir);
        meanarglist = {meandir,fittedrfx.meanrsa,grouprfx.meanrsa,...
            fitted.otherpredictors.meanrsa};
        facedist_rsa_plots_meanslopes(meanarglist{:});

        % GLM
        [rfxres_uni,groupres_uni] = roidata_rfx(fitted.subres.(fstr).glm,...
            'targetfield','b',...
            'nboot',0,'nperm',0,'assumeregister',0);
        % NB we recycle rsa_plots_meanslopes
        glmdir = fullfile(figdir,'meanglm');
        mkdirifneeded(glmdir);
        glmarglist = {glmdir,rfxres_uni,groupres_uni,[],'mtarget','mean',...
            'errtarget',[],'mlabel','mean population response',...
            'errlabel','','groupmtarget','b'};
        facedist_rsa_plots_meanslopes(glmarglist{:});

        % scatter showing univariate response vs model prediction
        % Compare to fMRI - need to 
        % c) intersect the ROIs to match
        % then find the intersection
        [sharedroi,roiinds_model,roiinds_fmri] = intersect(...
            groupres_uni.cols_roi,...
            groupres_uni_fmri.cols_roi);
        % ( work out which rows to plot)
        [~,~,coninds_model] = intersect({'all_mean_ecc_sub',...
            'all_mean_ecc_typical','all_mean_ecc_caricature'}',...
            groupres_uni.rows_contrast,'stable');
        % check that subjects are in register
        assert(isequal(groupres_uni.z_subject,...
            groupres_uni_fmri.z_subject),'mismatched subjects');
        % d) make a scatter figure for each - plot one slope per subject to
        % show reliability within-subjects
        F = figurebetter([],par.fitunisize);
        for n = 1:numel(sharedroi)
            roistr = sharedroi{n};
            modeldata = squeeze(groupres_uni.b(coninds_model,...
                roiinds_model(n),:))';
            fmridata = squeeze(groupres_uni_fmri.b(coninds_fmri,...
                roiinds_fmri(n),:))';
            valid = ~any(isnan(fmridata),2);
            modelcomp_r.(fstr).group.(strrep(roistr,' ','_')) = corrpairs(...
                modeldata(valid,:)',fmridata(valid,:)');
            % use plotarray so we can combine over subjects later
            scph = NaN(size(modeldata));
            colorscheme = facedist_colors('ecc');
            for ecc = 1:3
                scph(:,ecc) = plotarray(modeldata(:,ecc),...
                    fmridata(:,ecc),par.markertype,...
                    'color',colorscheme(ecc,:),...
                    'markersize',par.markersize,...
                    'markerfacecolor',colorscheme(ecc,:),...
                    'markeredgecolor','none');
            end
            % then slopes over the subjects in gray
            lh = NaN([size(scph,1),1]);
            for s = 1:size(scph,1)
                if all(isnan(scph(s,:)))
                    continue;
                end
                lh(s) = lslinebetter(scph(s,:),true,'adaptive','color',...
                    [.5 .5 .5],'linewidth',.5);
                uistack(lh(s),'bottom');
            end
            set(gca,'tickdir','out','ticklength',par.ticklength,...
                'position',par.axscale * [.8 1 .8 .8]);
            lims = getdatalims(scph,'xy');
            axis(gca,lims);
            roundplotlimits(gca,'xy',0,1);
            minimalticks(gca,'xy',0);
            plotorigin(gca);
            lh = legend(scph(find(~isnan(scph(:,1)),1,'first'),:),...
                {'sub','typical','caricature'},'location',...
                'northeastoutside');
            box(lh,'off');
            ylabel(par.ylabel_glm);
            xlabel({'model response',...
                '(population average)'});
            printstandard(fullfile(figdir,sprintf(...
                'scatter_modelcomp_%s',strrep(roistr,' ','_'))));
            clf(F);
        end
    end
    save(fitfile,'fitted')
    save(fullfile(resdir,'modelcomp_r.mat'),'modelcomp_r');
else
    fitted = loadbetter(fitfile);
end

% UNIVARIATE MODEL COMPARISON
% so each subres will contain a LOO r estimate, with one row per model and
% one column per ROI.
unireprofile = fullfile(resdir,'unirepro_withmodels.mat');
% work out correspondence between univariate fMRI result and ROI order 
[~,repind,uniind] = intersect(repro.cols_roi,groupres_uni_fmri.cols_roi,...
    'stable');
if ~exist(unireprofile,'file') || recompute
    unirep = struct;
    unirep.cols_roi = groupres_uni_fmri.cols_roi(uniind);
    % zscore dim 1 and then average
    nanzm = @(x,d)nanmean(zscore(x),d);
    unirep.data = groupres_uni_fmri.b(1:nface,uniind,:);
    for sub = 1:nsub
        trainind = setdiff(1:nsub,sub);
        unirep.subres(sub) = struct('rows_contrast',{{}},...
            'cols_roi',{unirep.cols_roi},'b',[],'tail',{{}},...
            'name',repro.subres(sub).name);
        % upper noise
        unirep.subres(sub).b(end+1,:) = atanh(...
            corrpairs(unirep.data(:,:,sub),nanmean(unirep.data,3)));
        unirep.subres(sub).rows_contrast{end+1,1} = 'rep_upper';
        unirep.subres(sub).tail{end+1,1} = 'none';
        % lower
        unirep.subres(sub).b(end+1,:) = atanh(...
            corrpairs(unirep.data(:,:,sub),...
            nanmean(unirep.data(:,:,trainind),3)));
        unirep.subres(sub).rows_contrast{end+1,1} = 'rep_lower';
        unirep.subres(sub).tail{end+1,1} = 'right';
        % add the CV performance for each model
        for fn = fieldnames(res.group)'
            fnstr = fn{1};
            % train data
            trainresp = matfun(@nanmean,...
                res.singles.(fnstr)(trainind).meanresponse);
            % remove facepairs
            trainresp = trainresp(:,repind);
            % Fisher transform
            unirep.subres(sub).b(end+1,:) = atanh(corrpairs(...
                unirep.data(:,:,sub),trainresp));
            unirep.subres(sub).rows_contrast{end+1,1} = fnstr;
            unirep.subres(sub).tail{end+1,1} = 'right';
        end
    end
    [unirep.rfx,unirep.group] = roidata_rfx(unirep.subres,'nperm',1024,...
        'nboot',1000,'contrasts','roidata_allpairwisecontrasts',...
        'targetfield','b');
    save(unireprofile,'unirep');
else
    unirep = loadbetter(unireprofile);
end

%% FIGURES
rootfigdir = fullfile(resdir,'figures');
mkdirifneeded(rootfigdir);
% global visualisation using repro.rfx and fitted
% ratio plot

%% BAR CHART COMPARING MODELS (using repro and plot_roidata)
barcon = {'ramp','ex_gaussian','ex_negativegaussian','gwpgrid',...
    'vid_corrdist','full'};
[~,barfigdir] = mkdirifneeded(fullfile(rootfigdir,'bars_repro'));
F = plotreprores(repro,barfigdir,barcon,'r','distance-matrix similarity');
close(F);

% and same thing with univariate
barconuni = {'ramp','ex_gaussian','gwpgrid'};
[~,barfigdiruni] = mkdirifneeded(fullfile(rootfigdir,'bars_unirepro'));
F = plotreprores(unirep,barfigdiruni,barconuni,'b','activation-profile similarity');
close(F);

[~,unirespdir] = mkdirifneeded(fullfile(rootfigdir,'bars_responses'));
fhand = figure;
resptargets = {'ramp','ex_gaussian','gwpgrid','ex_negativegaussian'};
respgroups = struct('suff',{'','_noav'},'line1',{[.9 .2 0],[.9 .2 0]},...
    'line2',{[.9 .2 0],[1 1 1]});
for roiind = 1:numel(unirep.cols_roi)
    m = zscore(nanmean(unirep.data(:,roiind,:),3));
    ax = axes('position',(1/3) * [1,1,1,1]);
    faceh = plot(1:nface,m,'marker',par.markertype,'markersize',par.markersize,...
        'markerfacecolor',[1 1 1],'markeredgecolor',[0 0 0],...
        'linestyle','none');
    ylim([-3 3]);
    set(ax,'dataaspectratio',[4 1 1],'xtick',1:nface,'xlim',...
        [.5 nface+.5],'tickdir','out','ticklength',par.ticklength,...
        'xcolor',[1 1 1]);
    % set limits also perhaps
    [outh,lineh] = imageticks(gca,{stimres.stimuli{1}(1:nface).image},'nrows',4,...
        'dim','x','autopos',false);
    % add some supporting lines
    gridx = par.gridlines+.5;
    lh = arrayfun(@(x)line([x,x],[-3,3],'color','k','linestyle',':',...
        'linewidth',.5),gridx,'uniformoutput',1);
    axis(ax,'on');
    box(ax,'off');
    minimalticks(ax,'y',0,0);
    ylabel('region mean response (Z)');
    printstandard(fullfile(unirespdir,sprintf('resp_dataonly_%s',...
        unirep.cols_roi{roiind})));
    hold(ax,'on');
    for target = resptargets
        modelstr = target{1};
        ph = [];
        for thisgroup = respgroups
            % response with averaging
            modelresp = matfun(@nanmean,...
                res.singles.([modelstr thisgroup.suff]).meanresponse);
            % nb extra layer of indexing to get around absence of facepairs
            modelresp = modelresp(:,repind(roiind));
            scaledresp = zscore(modelresp);
            % first line
            ph(end+1) = line(1:nface,scaledresp,'color',thisgroup.line1,...
                'linewidth',2);
            % second line
            ph(end+1) = line(1:nface,scaledresp,'color',thisgroup.line2,...
                'linewidth',.5);
        end
        ylim([-3 3]);
        uistack(faceh,'top');
        printstandard(fullfile(unirespdir,sprintf('resp_model_%s_%s',...
            unirep.cols_roi{roiind},modelstr)));
        delete(ph);
    end
    clf(fhand);
end
close(fhand);

barorder = {'perceptual judgments','early visual cortex',...
    'occipital face area','transverse occipital sulcus',...
    'fusiform face area','parahippocampal place area'};
%% prepare the data RDMs in the desired format
datardms = rdm2struct(asrdmmat(nanmean(repro.data,3)),repro.cols_roi);
printname = facedist_names(datardms.name);
[datardms.printname] = deal(printname{:});
[datardms.edgecolor] = deal([1 1 1]);
[datardms.facecolor] = deal([.7 .7 .7]);

nunit = 10;
unitlims = [-2.1; 2.1] * max(abs(xy));
p = .1;
unitx = unitlims(1,1):p:unitlims(2,1);
unity = unitlims(1,2):p:unitlims(2,2);

% Visualise the model parameters
for thisst = modelst
    fstr = thisst.field;
    [~,figdir] = mkdirifneeded(fullfile(resdir,fstr));
    dimlabels = cellfun(@(x)alldimlabels.(x),fieldnames(thisst.coord),...
        'uniformoutput',0);
    dimcoords = cellfun(@(x)thisst.coord.(x),fieldnames(thisst.coord),...
        'uniformoutput',0);
    dimind = 1:numel(dimlabels);
    if thisst.skipntune
        dimind(1) = [];
    end

    thisres = res.all.(fstr);
    thissingle = res.singles.(fstr);
    % run RFX analysis - nb original res order
    [thissingle.cols_roi] = deal(repro.cols_roi);
    [thissingle.rows_contrast] = deal(dimlabels);
    [thissingle.tail] = deal(repmat({'both'},[numel(dimlabels),1]));
    [thissingle.name] = deal(stimres.meanrsa.name);
    args = {thissingle,'targetfield','fitpar',...
        'roicontrasts','roidata_allpairwisecontrasts','nboot',1000,...
        'nperm',1024,'assumeregister',0};
    rfxres = roidata_rfx(args{:});

    if strcmp(fstr,'ramp')
        responses = struct('final',[],'raw',[]);
        if ~all(rfxres.mean(:,end)==0)
            responses.average = [];
        end
    end

    F = figure;
    for t = 1:nall
        prargs = {fstr,'',datardms(t).name};
        % search grid
        if isfield(thisres,'rmedian')
            est = thisres.rmedian{barinds(t)};
            lims = [0 1];
            sz = thisst.sz;
            if thisst.skipntune || thisst.skipaverage
                est = squeeze(est);
            end
            if thisst.skipntune
                sz(1) = [];
            end
            h = imagend(est,ind2subbetter(sz,thisres.winner(barinds(t))),...
                dimlabels(dimind),dimcoords(dimind),lims);
            colormap(jet(256));
            % the grid sampled each dim equally so this should make the
            % pixels approximately square.
            set(h(1).ax,'plotboxaspectratio',[1 1 1]);
            xlabel(h(end).cb,{'similarity','median r'});
            minimalticks(h(end).cb,'x');
            set([h(1:end-1).ax],'plotboxaspectratio',[1 1 1]);
            printstandard(fullfile(figdir,sprintf(...
                'diagnostic_grid_%s_%s_%s',prargs{:})),'formats=eps');
            clf(F);
        end

        if ~any(strcmp(fstr,{'gwpgrid','gwpgrid_noav'}))
            model = buildmodel(rfxres.mean(:,t),thisst.sampler,thisst.instance);

            if strcmp(fstr,'ramp')
                [responses.final(:,end+1),xvals,responses.raw(:,end+1)] = gettestresponse(model);
            end

            % example units
            Fgrid = figurebetter([],[12 12]);
            plotexampleunits(model,nunit,unitx,unity,scalefactor,xy,figdir,...
                prargs,Fgrid);
            if isa(model,'RampModel') && model.averaging>0
                thisa = model.averaging;
                tempfigdir = fullfile(figdir,'averaging1');
                mkdirifneeded(tempfigdir);
                model.averaging = 1;
                % plot 2 units with full averaging
                plotexampleunits(model,2,unitx,unity,scalefactor,xy,...
                    tempfigdir,prargs,Fgrid);
                % add response for tuning function
                responses.average(:,end+1) = gettestresponse(model);
                model.averaging = 0;
                % plot all units with no averaging
                tempfigdir = fullfile(figdir,'averaging0');
                mkdirifneeded(tempfigdir);
                plotexampleunits(model,nunit,unitx,unity,scalefactor,xy,...
                    tempfigdir,prargs,Fgrid);
                % restore original par
                model.averaging = thisa;
            end

            % example dissimilarity profiles
            for xind = 1:numel(xyplotind)
                [h,ph,intmap,cmap] = plotdissimilarity(model,...
                    xy(xyplotind(xind),:),'xvals',unitx,...
                    'yvals',unity,'ax',gca,'limits','zerobounded');
                plotunitgrid(scalefactor,xy);
                uistack(ph,'top');
                setplotz(ph,2);
                set(ph,'marker','o','markersize',10,...
                    'markerfacecolor',[1 1 1],...
                    'markeredgecolor',[0 0 0]);
                % with CB
                newpos = get(gca,'position')*.75;
                set(gca,'position',newpos);
                cbax = axes('position',...
                    [newpos(1:2) + [newpos(3) 0] .2 .2]);
                blim = [intmap(1) reduceprecision(intmap(end),1)];
                colorbarbetter(cbax,blim,cmap,'orientation',...
                    'horizontal','label','euclidean distance',...
                    'scale',.5,'tick',[]);
                printstandard(fullfile(figdir,sprintf(...
                    'diagnostic_dissimilarityprofilewcb_%s_%s_%s_%02d',...
                    prargs{:},xind)),'formats=eps');
                clf(Fgrid);
            end
            close(Fgrid);
        end

        % mean RDM
        [ax,intmap,cmap] = rdmplot(gca(F),thisres.rdvwinmean(:,barinds(t)),...
            'gridlines',par.gridlines,'cmap',cmap_wr,...
            'gridcolor',[1 1 1],'limits','zerobounded');
        % with colorbar
        newpos = get(gca(F),'position')*.75;
        set(gca,'position',newpos);
        cbax = axes('position',...
            [newpos(1:2) + [newpos(3) 0] .2 .2]);
        colorbarbetter(cbax,intmap,cmap,'orientation','horizontal',...
            'label','euclidean distance','scale',.5,'tick',[]);
        % round the tick labels a bit 
        newlim = roundplotlimits(cbax,'x',[],0);
        set(cbax,'xlim',newlim,'xticklabelmode','auto','xtick',newlim);
        printstandard(fullfile(figdir,sprintf(...
            'rdm_meanfitwcb_%s_%s_%s',prargs{:})),'formats=eps');
        clf(F);
    end % t nall
    close(F);

    prargs = {fstr,''};

    if strcmp(fstr,'ramp')
        % make sigmoid response visualisation
        F = figurebetter([7 7]);
        xt = [-1 * xy_ecc([9 5 1])' 0 xy_ecc([1 5 9])'];
        % bit hacky but fitlinecolor is black and that's what we want the
        % norm to be.
        xl = {'caricature','typical','sub','fitlinecolor','sub','typical',...
            'caricature'};
        linestyles = repmat({'-'},[1 length(xl)]);
        linewidths = ones(1,length(xl));
        linestyles(4) = {':'};
        linewidths(4) = .5;
        xticklab = repmat({''},[1 length(xl)]);
        xticklab{4} = '0';
        % separately for normalised and raw
        for fn = fieldnames(responses)'
            respstr = fn{1};
            set(gcf,'defaultaxescolororder',facedist_colors(...
                datardms(1:nall).name));
            ph = plot(xvals,responses.(respstr));
            L = legend(ph,[facedist_names(datardms(1:nall).name)],...
                'location','southeast');
            set(L,'box','off');
            set(gca,'position',[.25 par.axscale * [1 1 1]],'tickdir',...
                'out','box','off','xtick',xt,'xticklabel',xticklab);
            lh = [];
            ylim([0 1]);
            for xpos = 1:7
                xx = xt(xpos);
                lh(end+1) = line([xx xx],[0 1],...
                    'linewidth',linewidths(xpos),...
                    'color',facedist_colors(xl{xpos}),...
                    'linestyle',linestyles{xpos});
                uistack(lh(end),'bottom');
            end
            ylabel([respstr ' response']);
            xlabel('face space position');
            xlim([min(xvals) max(xvals)]);
            lh = axes('position',[par.axscale+.1 par.axscale * [1 1 1]]);
            axis(lh,'off');
            centerinaxis(L,lh);
            printstandard(fullfile(figdir,sprintf(...
                'tuningfunctions_%s_%s_%s',respstr,prargs{:})));
            clf(F);
        end
    end

    % bar chart with winning fit parameters. The units here don't make
    % sense so need one small bar chart per metric
    F = figurebetter([],par.fitratiofsize,'auto');
    inds = 1:thisst.nd;
    if thisst.skipntune
        inds(1) = [];
    end
    % one pmat per dimension (each row/column is a region)
    allpcon = asrdmmat(rfxres.ppara(:,nall+1:end)');

    % for each dimension
    for d = inds
        % plot mean mean fitpar from singles with standard errors and
        % contrasts instead. 
        m = rfxres.mean(d,barinds);
        err = rfxres.sterr(d,barinds);
        % dump out to CSV
        pcon = allpcon(:,:,d);
        pcon = pcon(barinds,barinds,:);
        % dump out p values for pairwise comparisons
        rdm2csv(pcon,fullfile(figdir,sprintf('bars_fit_%s_%s_%s_pcon.csv',...
            stripbadcharacters(dimlabels{d}),prargs{:})),...
            {datardms(1:nall).printname});
        pcon(pcon>.05) = NaN;
        barchart(m,'labels',...
            {datardms(1:nall).printname},'errors',err,...
            'facecolor',cat(1,datardms(1:nall).facecolor),...
            'edgecolor',[],'fighand',F,'rotatelabels',45);
        set(gca,'position',repmat(par.axscale,[1 4]),...
            'tickdir','out',...
            'ticklength',par.ticklength);
        rotateXLabels(gca,45);
        ylabel(dimlabels{d});
        yl = ylim;
        % lock ylim
        ylim(yl);
        contrastlines(gca,pcon);
        printstandard(fullfile(figdir,sprintf('bars_fit_%s_%s_%s',...
            stripbadcharacters(dimlabels{d}),prargs{:})),'formats=eps');
        clf(F);
    end
end

logstr('finished run in %s\n',seconds2str(etime(clock,runstart)));
logstr('results are at %s\n',resdir);
keyboard;

%%% SUB FUNCTIONS %%%

% sample points on a circle of fixed eccentricity
function xy = circsample(n)

rads = rand(n,1)*2*pi;
% the eccentricity of the caricatures
eccs = ones(n,1) * 2.4042;
[x,y] = pol2cart(rads,eccs);
xy = [x,y];

% [rdv,meanresp] = fitit(fitpar)
function [rdv,meanresp] = fitit(fitpar,sampler,testunits,instanceclass)

niter = 100;
ntest = size(testunits,1);

rdv = cell(1,niter);
meanresp = cell(1,niter);
parfor n = 1:niter
    % index the right number of units for this fitpar
    % create the instance
    model = buildmodel(fitpar,sampler,instanceclass);
    % get the response distance matrix
    resp = NaN([ntest,size(model.exemplars,1)]);
    for t = 1:ntest
        resp(t,:) = populationresponse(model,testunits(t,:));
    end
    rdv{n} = pdist(resp)';
    meanresp{n} = mean(resp,2);
end
rdv = cat(2,rdv{:});
meanresp = cat(2,meanresp{:});

function model = buildmodel(fitpar,sampler,instanceclass)

sampler = str2func(sampler);
% first par is always number of units
fitpar(1) = ceil(10^fitpar(1));
% last par is always population averaging proportion
averaging = fitpar(end);
fitpar(end) = [];
if isequal(sampler,@facedist_chi2sample) || isequal(sampler,@facedist_norminvsample)
    % use first 2 pars to tune sampling scheme
    units = sampler(fitpar(1),fitpar(2));
    fitpar(1:2) = [];
else
    % just one input - n
    units = sampler(fitpar(1));
    fitpar(1) = [];
end

if strcmp(instanceclass,'RampModel')
    % each input is passed to instance
    args = num2cell(fitpar);
else
    % the next par is spacing of points
    units = units * fitpar(1);
    % the rest of the arguments are tuning widths and go into a single
    % vector
    args = {fitpar(2:end)};
end
model = feval(instanceclass,units,args{:},averaging);

function res = addtempres(tempres,bestind)

ntarget = size(tempres.fitpar,2);
for t = 1:ntarget
    res(1).fitpar(:,t) = ...
        tempres.fitpar(:,t,bestind(t));
    res.minrdist(t) = ...
        tempres.minrdist(t,bestind(t));
    res.exitflag(t) = ...
        tempres.exitflag(t,bestind(t));
    res.rdvwinmean(:,t) = ...
        tempres.rdvwinmean(:,t,bestind(t));
    res.rdvwin(:,t) = ...
        tempres.rdvwin(:,t,bestind(t));
    res.iterations(t) = ...
        tempres.output(t,bestind(t)).iterations;
    res.funcCount(t) = ...
        tempres.output(t,bestind(t)).funcCount;
end
res.bestind = bestind;

function plotunitgrid(scalefactor,xy)
% colours and plot pars for gridmds style plot
% this is all adapted from facedist_gridmds, which turns out to be too
% messy to re-use here.
colors = facedist_colors('ecc');
args = arrayfun(@(x){'linewidth',2,'marker','o','markersize',10,...
    'markerfacecolor',colors(x,:),'markeredgecolor',colors(x,:),...
    'linestyle','none'},...
    3:-1:1,'uniformoutput',0);
[dirrdm,eccrdm] = facedist_rsapredictor_neighborsonly(12);
[s,r] = facedist_gridpanel(gca,false,[],scalefactor);
ehand = plotmdslines(gca,xy(1:12,:),eccrdm,colors,false,args{:});
setplotz(s,.5);
setplotz(r,1);
setplotz(ehand,1.5);
uistack(r,'top');

function thisres = rungridfit(datardv,thisopt,thisrdv)
ntarget = size(datardv,2);
onepad = ones(1,thisopt.nd);
sd = thisopt.nd+1;
thisres.rdvwin = NaN([size(thisrdv{1},1),ntarget],'like',thisrdv{1});
thisres.rdvwinmean = NaN([size(thisrdv{1},1),ntarget],'like',thisrdv{1});
thisres.winner = NaN([1 ntarget]);
thisres.fitpar = NaN([numel(fieldnames(thisopt.grid)) ntarget]);
for t = 1:ntarget
    thistarget = datardv(:,t);
    % drop nans from vector and thisrdv.
    nanind = ~isnan(thistarget);
    if ~any(nanind)
        logstr('skipping nan ROI: %d\n',t);
        continue
    end
    thistarget = thistarget(nanind);
    thisr = cell(thisopt.sz);
    % get r for each position
    parfor n = 1:thisopt.n
        % store with singletons to support cell2mat
        % pad in as many singleton dims as we need to put
        % the iterations last
        thisr{n} = reshape(pearsonvec(thistarget,...
            thisrdv{n}(nanind,:)),[onepad size(thisrdv{n},2)]);
    end
    rmat = cell2mat(thisr);
    % placeholder
    thisres.meanresponse = [];
    thisres.rmean{t} = mean(rmat,sd);
    thisres.rstd{t} = std(rmat,[],sd);
    thisres.rz{t} = mean(rmat,sd) ./ ...
        std(rmat,[],sd);
    thisres.rmedian{t} = median(rmat,sd);
    % find the winner (largest median r) - we could also
    % use Z or mean here
    thisres.winner(t) = ...
        find(thisres.rmedian{t} == max(...
        thisres.rmedian{t}(:)),1,'first');
    % and store an example RDV for it (the first)
    % handle potential nan
    thisres.rdvwin(nanind,t) = ...
        thisrdv{thisres.winner(t)}(nanind,1);
    % and the mean (who knows, might work better...)
    thisres.rdvwinmean(nanind,t) = ...
        mean(thisrdv{thisres.winner(t)}(nanind,:),2);
    % fitpars
    thisres.fitpar(:,t) = structfun(...
        @(thisv)thisv(thisres.winner(t)),...
        thisopt.grid);
end % t ntarget

function [resp,xv,rawresp] = gettestresponse(model)
v = 3;
xv = linspace(-v,v,200)';
resp = getresponse(model,[v 0],[xv zeros(numel(xv),1)]);
rawresp = getrawresponse(model,[v 0],[xv zeros(numel(xv),1)]);

function plotexampleunits(model,nunit,unitx,unity,scalefactor,xy,figdir,prargs,Fgrid)

for u = 1:nunit
    % plot unit
    [h,ph,intmap,cmap] = plotunit(model,u,'xvals',unitx,...
        'yvals',unity,'ax',gca,'cmap',cmap_bgr,'limits',[0 1]);
    plotunitgrid(scalefactor,xy);
    uistack(ph,'top');
    setplotz(ph,2);
    otheru = model2data(model,model.exemplars(setdiff(...
        1:nunit,u),:));
    phother = plot(otheru(:,1),otheru(:,2),'o-');
    if isa(model,'RampModel')
        set(ph,'marker','none','color',[1 1 1]);
        set(phother,'marker','none','color',[0 0 0]);
    else
        set(ph,'marker','o','markersize',10,...
            'markerfacecolor',[1 1 1],...
            'markeredgecolor',[0 0 0],'linestyle','none');
        set(phother,'marker','o','markersize',10,...
            'markerfacecolor','none',...
            'markeredgecolor',[0 0 0],'linestyle','none');
    end
    uistack(phother,'top');
    uistack(ph,'top');
    setplotz(ph,3);
    setplotz(phother,3);
    % with CB
    newpos = get(gca,'position')*.75;
    set(gca,'position',newpos);
    cbax = axes('position',...
        [newpos(1:2) + [newpos(3) 0] .2 .2]);
    colorbarbetter(cbax,intmap,cmap,'orientation','horizontal',...
        'label','response','scale',.5,'tick',[]);
    printstandard(fullfile(figdir,sprintf(...
        'diagnostic_unitresponsewcb_%s_%s_%s_%02d',prargs{:},...
        u)),'formats=eps');
    delete(phother);
    printstandard(fullfile(figdir,sprintf(...
        'diagnostic_unitresponsewcbnoother_%s_%s_%s_%02d',prargs{:},...
        u)),'formats=eps');
    clf(Fgrid);
end

function F = plotreprores(repro,barfigdir,barcon,targetfield,ylab)

par = facedist_plotpar;

repro.rfx.cols_roi = facedist_names(repro.rfx.cols_roi{:});
repro.group.cols_roi = facedist_names(repro.group.cols_roi{:});
pmat = asrdmmat(repro.rfx.ppara(strfindcell(repro.rfx.rows_contrast,...
    'contrast_'),:));
% color scheme - we want filled bars for fitted and outlined bars for fixed
% models. And we want distinct colours for models that operate in face
% space and image space.
barfaces = facedist_colors(barcon{:});

barwidth = .4;
suffixes = struct('name',{'av','noav'},'suff',{'','_noav'},...
    'facecolor',{barfaces,ones(size(barfaces))},...
    'edgecolor',{barfaces,barfaces},...
    'xoff',{-barwidth/2,barwidth/2});
nbgroup = numel(barcon);
xpos = 1:nbgroup;
for suff = suffixes
    % handle missing cases
    conind = cellfun(@(b)find(strcmp(...
        [b suff.suff],repro.rfx.rows_contrast)),barcon,'uniformoutput',0);
    missing = cellfun(@isempty,conind);
    conind(missing) = [];
    conind = cat(1,conind{:});
    data.(suff.name).ind = xpos(~missing);
    data.(suff.name).x = xpos(~missing) + suff.xoff;
    data.(suff.name).m = tanh(repro.rfx.mean(conind,:));
    data.(suff.name).ms = tanh(repro.group.(targetfield)(conind,:,:));
    data.(suff.name).p = repro.rfx.ppara(conind,:);
    data.(suff.name).pmat = pmat(conind,conind,:);
    data.(suff.name).errs = tanh(repro.rfx.sterr(conind,:));
    data.(suff.name).conind = conind;
    data.(suff.name).facecolor = suff.facecolor(~missing,:);
    data.(suff.name).edgecolor = suff.edgecolor(~missing,:);
end
% post flight adjustment on bar widths for singles
singles = setdiff(round(data.av.x),round(data.noav.x));
data.av.x(singles) = round(data.av.x(singles));

allx = [data.av.x,data.noav.x];
[allx,sortind] = sort(allx);
allcon = [data.av.conind;data.noav.conind];
allcon = allcon(sortind);
allp = pmat(allcon,allcon,:);
notsingle = round(allx) ~= allx;

% write out one pcon table with roi dim stacked vertically.
rdm2csv(allp,fullfile(barfigdir,'pcon_modelcomp.csv'),...
    facedist_names(repro.rfx.rows_contrast{allcon}),...
    'zlabels',repro.rfx.cols_roi,'precision',4);

F = figurebetter([],par.fitratiofsize .* [1.4 .8],'auto');
xl = [.5 nbgroup+.5];
plotwidth = par.axscale * (nbgroup / 6);
plotpos = [par.axscale par.axscale*1.3 plotwidth par.axscale];
scolor = facedist_colors('singles');
for t = 1:size(repro.data,2)
    % separate plot calls for the two groups
    for suff = suffixes
        [F,handles.(suff.name).bars,handles.(suff.name).errs,...
            handles.(suff.name).ptext] = barchart(...
            data.(suff.name).m(:,t)',...
            'facecolor',data.(suff.name).facecolor,...
            'edgecolor',data.(suff.name).edgecolor,...
            'fighand',F,...
            'pvalues',data.(suff.name).p(:,t)',...
            'errors',data.(suff.name).errs(:,t)','plotbaseline',false,...
            'x',data.(suff.name).x,'width',barwidth);
    end
    set(gca,'xtick',xpos,'xlim',[0.5 nbgroup+.5],'xcolor',[1 1 1]);
    lh = line(xlim,[0,0],'color','k','linestyle',':','linewidth',.5);
    uistack(lh,'bottom');
    ylabel({ylab,'(mean r\pm1 standard error)'});
    % add noise ceiling
    noise_low = tanh(repro.rfx.mean(strcmp(repro.rfx.rows_contrast,...
        'rep_lower'),t));
    noise_high = tanh(repro.rfx.mean(strcmp(repro.rfx.rows_contrast,...
        'rep_upper'),t));
    handles.noiseceil = errorshade(xlim(gca)',noise_low * [1;1],...
        noise_high .* [1; 1],facedist_colors('reproducibility'),0);
    uistack(handles.noiseceil,'bottom');
    th = text(max(xlim(gca)) + range(xlim(gca))*.05,mean([noise_low ...
        noise_high]),'noise ceiling','horizontalalignment',...
        'left','verticalalignment','middle');
    finaly = sety([handles.av.errs,handles.noav.errs],noise_high);
    set(gca,'position',plotpos,'tickdir','out','ticklength',par.ticklength);
    th = text(xpos,repmat(finaly(1)-.1*range(finaly),size(xpos)),...
        facedist_names(barcon{:}),'rotation',45,...
        'horizontalalignment','right','verticalalignment','middle');
    thispcon = allp(:,:,t);
    % add contrast lines
    % restrict to just comparing av vs noav for simplicity.
    thispcon = thispcon(notsingle,notsingle);
    thispcon(~blkinds([sum(notsingle),sum(notsingle)],2)) = NaN;
    thispcon(thispcon>.05) = NaN;
    [handles.con,newy] = contrastlines(gca,thispcon,'ypos',finaly(2) + ...
        range(finaly)*.1,'xpos',allx(notsingle));
    % and by now we probably need to shift the p values to the new
    % limits
    settextpos([handles.av.ptext;handles.noav.ptext],...
        max([finaly(1) newy]) + .1*range(finaly),'y');
    printstandard(fullfile(barfigdir,['bars_modelcomp_',...
        repro.rfx.cols_roi{t}]));
    % singles version - need to iterate over bars instead and plot each
    deleten(handles.con);
    deleten([handles.av.errs,handles.noav.errs]);
    for x = xpos
        ys = asrow(data.av.ms(x,t,:));
        xs = data.av.x(x);
        noavind = find(data.noav.ind == x);
        if any(noavind)
            ys(end+1,:) = data.noav.ms(noavind,t,:);
            xs(end+1,1) = data.noav.x(noavind);
        end
        singhand{x} = plot(xs,ys,'markersize',par.markersize,...
            'marker',par.markertype,'linestyle','-',...
            'linewidth',.5,'markerfacecolor',[1 1 1],...
            'markeredgecolor',scolor,'color',scolor);
    end
    ylabel({ylab,'(r)'});
    finaly = sety(cat(2,singhand{:}),noise_high);
    [handles.con,newy] = contrastlines(gca,thispcon,...
        'ypos',finaly(2) + .1*range(finaly),...
        'xpos',allx(notsingle));
    settextpos([handles.av.ptext;handles.noav.ptext],...
        max([finaly newy]) + .1*range(finaly),'y');
    settextpos(th,finaly(1)-.1*range(finaly),'y');
    printstandard(fullfile(barfigdir,['singles_bars_modelcomp_',...
        repro.rfx.cols_roi{t}]));
    clf(F);
end
% export the basic model comp stats too
exportres = repro.rfx;
exportres.cols_roi = facedist_names(exportres.cols_roi{:});
exportres.rows_contrast = facedist_names(exportres.rows_contrast{:});
exportres = renamefield(exportres,'mean','mean_zr');
exportres = renamefield(exportres,'sterr','sterr_zr');
exportres = renamefield(exportres,'ppara','ppara_zr');
exportres.mean_r = tanh(exportres.mean_zr);
roidata2csv(exportres,fullfile(barfigdir,'modelcomp_stats.csv'),...
    allcon,[],3,'mean_r','mean_zr','sterr_zr','ppara_zr','n');

function finaly = sety(hand,noise_high)

ydata = getdatalims(hand,'y');
ydata(2) = max([ydata(2),noise_high]);
ylim(ydata);
finaly = roundplotlimits(gca,'y',0,1);
minimalticks(gca,'y',0,1);
