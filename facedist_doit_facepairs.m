% master analysis for behavioral similarity judgment task. The basic idea is to
% calculate RDMs, run them through the same analysis pipeline that is used
% for the fMRI analysis (see facedist_aa_frombids), and then plug everything
% into the computational models (facedist_doit_models).
%
% 2017-06-28 J Carlin
%
% function facedist_doit_facepairs([redostats=false],[varargin])
function facedist_doit_facepairs(redostats,varargin)

if ieNotDefined('redostats')
    redostats = false;
end

getArgs(varargin,{'nperm',1024,'nboot',1000});

rootdir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
derivdir = fullfile(rootdir,'derivatives');
[~,outdir] = mkdirifneeded(fullfile(derivdir,'facepairs'));
fmriresdir = fullfile(derivdir,'aa');
stimdir = fullfile(fmriresdir,'aamod_pilab_importstimuli_00001');
rdmdir = fullfile(fmriresdir,'aamod_pilab_importrsapredictors_00001');
outres=  fullfile(outdir,'stimanalysis.mat');

subjects = mat2strcell(1:10,'sub-%02d');
nsub = numel(subjects);

target = 'facepairs';

if ~exist(outres,'file') || redostats
    for s = 1:nsub
        subname = subjects{s};
        stimpath = fullfile(stimdir,subname,'pilab','pilab_stimuli.mat');
        res.stimuli{s} = loadbetter(stimpath);
        rdms = loadbetter(fullfile(rdmdir,subname,'pilab',...
            'pilab_rsapredictors.mat'));
        % get the facepairs predictors - expressed as percent chosen (-50)
        res.singles{s} = 100 * (asrdmvec(rdms(...
            strcmp({rdms.name},target)).RDM(1:12,1:12)) / 65) - 50;
        assert(~any(abs(res.singles{s}) > 50),...
            'we assume 65 trials per pair');
        % collapse the views
        res.rdms_collapse{s} = facedist_rsapredictors_collapseview(rdms);
        disvol{s} = Volume(res.singles{s},'metafeatures',...
            struct('names',{{target}}));
        % single subject mean RSA (NB we overwrite every time because it's
        % the same for all subjects)
        res.rdms_mean = facedist_rsapredictors_mean_tuning(...
            res.rdms_collapse{s});
        res.meanrsa(s) = roidata_rsa(disvol{s},res.rdms_mean,...
            'rsaclass','MeanRSA','splitrsaclass','MeanRSA',...
            'customfits','facedist_customfits_mean',...
            'contrasts','facedist_rsa_rfx_contrasts_meanslopes',...
            'defaulttail','both');
        % make sure the mean contrasts are 2-tailed too
        res.meanrsa(s).tail(strfindcell(res.meanrsa(s).rows_contrast,...
            'mean_direction'))= {'both'};
        % single subject RSA GLM
        % adjust back to 50 so the model fit will have correct units
        facepairind = strcmp(disvol{s}.meta.features.names,'facepairs');
        disvol{s}.data(:,facepairind) = disvol{s}.data(:,facepairind)+50;
        res.rdms_rsaglm = facedist_rsapredictors_glm(res.rdms_collapse{s});
        res.rsaglm(s) = roidata_glm(res.rdms_rsaglm,disvol{s},...
            'contrasts',@roidata_allpairwisecontrasts,...
            'customfun',{'predictY','residuals','getdata'},...
            'glmclass','MultiRSA');
    end
    % RFX meanrsa 
    [res.meanrsa.name] = deal(subjects{:});
    [res.meanrfx,res.groupres_mean] = roidata_rfx(res.meanrsa,...
        'nperm',nperm,'nboot',nboot);
    % need to add 50 here instead (previously we added 50 AFTER running the
    % single subject MeanRSA)
    facepairind = strcmp(res.groupres_mean.cols_roi,'facepairs');
    res.groupres_mean.r(:,facepairind,:) = ...
        res.groupres_mean.r(:,facepairind,:) + 50;

    % RFX RSA GLM
    [res.rsaglm.name] = deal(subjects{:});
    [res.rsaglmrfx,res.groupres_rsaglm] = roidata_rfx(res.rsaglm,...
        'nperm',nperm,'nboot',nboot,'targetfield','b');

    save(outres,'res');
    fprintf('saved res to %s\n',outres);
else
    fprintf('loading existing result from %s\n',outres);
    res = loadbetter(outres);
end 

% Mean RSA visualisation 
meandir = fullfile(outdir,'meanrsa');
mkdirifneeded(meandir);
meanarglist = {meandir,res.meanrfx,res.groupres_mean,res.rdms_mean,...
    'mtarget','mean','errtarget','sterr','ptarget','ppara','pthresh',.05};

facedist_rsa_plots_meanslopes(meanarglist{:});

% visualise full RDM and RSA GLM components on the same scale (face space
% distance) for conceptual sum predictors = full figure (with little
% squares everywhere)
fullrdm = indexstructbyfield(res.rdms_collapse{1},'name','full');
allrdms = [fullrdm res.rdms_rsaglm];
lims = max(ascol(asrdmvec(allrdms))) .* [-1 1];
F = figure;
par = facedist_plotpar;
figdir = fullfile(outdir,'predictors_rsaglm');
mkdirifneeded(figdir);
for p = 1:numel(allrdms)
    [~,intmap,cmap] = rdmplot(gca,allrdms(p),'gridlines',...
        par.gridlines,'gridcolor',[1 1 1],'limits',lims);
    % plot without cb label
    printstandard(fullfile(figdir,['rdm_' allrdms(p).name]),'formats=eps');
    newpos = get(gca,'position')*.75;
    set(gca,'position',newpos);
    cbax = axes('position',...
        [newpos(1:2) + [newpos(3) 0] .2 .2]);
    colorbarbetter(cbax,intmap,cmap,'orientation','horizontal',...
        'label','euclidean face space distance','scale',.5,'tick',[]);
    % round the tick labels a bit 
    newlim = roundplotlimits(cbax,'x',[],0);
    set(cbax,'xlim',newlim,'xticklabelmode','auto','xtick',...
        sort([newlim 0]));
    printstandard(fullfile(figdir,['rdm_wcb_' allrdms(p).name]),...
        'formats=eps');
    clf(F);
    % plot with image labels
    rdmplot(gca,allrdms(p),'gridlines',par.gridlines,...
        'gridcolor',[1 1 1],'limits',lims,'labels',...
        {res.stimuli{1}(1:12).image},'nrows',2);
    printstandard(fullfile(figdir,['rdm_wlabels_' allrdms(p).name]),...
        'formats',{'eps','png'});
    clf(F);
end
close(F);

% RSA GLM visualisation
glmdir = fullfile(outdir,'rsaglm');
mkdirifneeded(glmdir);
glmargs = {glmdir,res.rsaglmrfx,res.groupres_rsaglm,res.rdms_rsaglm,...
    'mtarget','mean','groupmtarget','b','errtarget','sterr',...
    'ptarget','ppara',...
    'mlabel','parameter estimate','errlabel','mean\pm1 standard error'};
facedist_rsaglm_plots(glmargs{:});

% face grid visualisation
figdir = fullfile(outdir,'stimgrid');
mkdirifneeded(figdir);
% get the full RDM
targetrdm = res.rdms_collapse{1}(strcmp(...
    {res.rdms_collapse{1}.name},'full')).RDM;
% use first subject's stimuli for now (we should probably repeat this for
% all subjects actually so we can demonstrate how the stimulus generation
% works out).

F = figurebetter([],[12 12]);
% just plot the empty grid
ecchand = facedist_gridmds(F,targetrdm,[],true);
[spokes,rings] = facedist_gridpanel(gca,false);
printstandard(fullfile(figdir,'empty_plotmds_full'),'F',F,'formats=eps');
set(rings(:),'color',par.notsigcolor);
eccs = ecchand{:}(end:-1:1);
set(eccs,'marker','none');
% visualise absecc
for ecc = 1:3
    set(rings(ecc),'color',par.basecolors(ecc,:));
    set(eccs(ecc),'marker','o');
    uistack(eccs(ecc),'top');
    setplotz(eccs(ecc),1);
    printstandard(fullfile(figdir,...
        ['empty_plotmds_' par.absecclabels{ecc}]),'F',F,...
        'formats',{'eps'});
    % reset to gray
    set(eccs(ecc),'marker','none');
    set(rings(ecc),'color',par.notsigcolor);
end
clf(F);

% This doesn't really fit here but good place perhaps to demonstrate
% dissimilarity variance
F = figurebetter([],[8 8],1/2);
ax = subplot(2,2,1);
yscat = rand(nchoosek(12,2),1);
rdv = asrdmvec(res.rdms_collapse{1}(1).RDM);
plot(rdv,yscat,'.');
[m,ind] = max(rdv);
hold(ax(1),'on');
plot(rdv(ind),yscat(ind),'ok');
ylim([-1 2]);
th = text(rdv(ind)+range(xlim(ax(1)))*.05,yscat(ind)+range(ylim)*.05,...
    {'face-antiface','caricatures'},...
    'horizontalalignment','left','verticalalignment','bottom');
titleh(1) = title(ax(1),'morph vectors');
box(ax(1),'off');
ax(2) = subplot(2,2,3);
plot(pdist(randn(12,199)),rand(nchoosek(12,2),1),'.');
titleh(2) = title(ax(2),'gaussian sampling');
set(titleh,'fontweight','bold');
box(ax(2),'off');
matchaxis(ax);
xlabel({'dissimilarity','(euclidean face space distance)'});
set(ax,'ytick',[],'tickdir','out');
printstandard(fullfile(outdir,'dissimilarityvariance'),'formats=eps');
close(F);

F = figurebetter([],[8 8],1/2);
for s = 1:nsub
    stimuli = res.stimuli{s}(1:12);
    % plot with images
    facedist_gridmds(F,targetrdm,stimuli,true);
    [spokes,rings] = facedist_gridpanel(gca,false);
    printstandard(fullfile(figdir,['plotmdslines_' subjects{s}]),'F',F,...
        'formats',{'png'},'r',1200);
    clf(F);
end
close(F);

% RDM visualisation 
ts = struct('cmap',[],'gridlines',[4 8],'gridcolor',[1 1 1]);
facedist_rdmplots_view(rdm2struct(asrdmmat(matmean(res.singles{:})),...
    res.meanrfx.cols_roi),res.stimuli{1},outdir,ts);

logstr('results are at %s\n',outdir);
