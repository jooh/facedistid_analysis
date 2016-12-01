% Construct RDM predictors for GWP RSAGLM. At the moment we return one RDM
% per filter bank.
%
% [outrdms,meanresp,filterresp,leftind] = facedist_stim2gwppredictors(stimuli)
function [outrdms,meanresp,filterresp,leftind] = facedist_stim2gwppredictors(stimuli,outdir)

if strcmp(stimuli,'DEBUG')
    % just load up from somewhere
    stimuli = loadbetter(['/imaging/jc01/research/idanddistinct/' ...
        'results_fullsample_realign/aamod_pilab_importstimuli_00001/' ...
        'be/pilab/pilab_stimuli.mat']);
end

exportimages = true;
if ieNotDefined('outdir')
    exportimages = false;
end

% so patching from seyed's code (allModelsUntrained_func_J)
% squarify each image
maxdim = max(arrayfun(@(thisx)max(size(thisx.rawimage)),stimuli));
outdim = maxdim;
nstim = numel(stimuli);
outvec = NaN([nstim,outdim*outdim],'single');
for s = 1:nstim
    im = stimuli(s).rawimage;
    squareim = single(expandbackground(im,[maxdim,maxdim]));
    if outdim ~= maxdim
        squareim = imresize(squareim,[outdim outdim],'bicubic');
    end
    % ensure max is exactly 1
    squareim = squareim / max(squareim(:));
    outvec(s,:) = squareim(:);
    if exportimages
        imwrite(gray2ind(squareim,256),fullfile(outdir,sprintf(...
            'raw_im%02d.png',s)),'PNG');
    end
end

cpfovs_deg = [.5,1,2,4,8];
% the longest image dim in the scanner was 8.4 deg. GWP uses image cycles.
% So to convert our degree-based cycle sized to image cycle units, do
cpfovs_im = cpfovs_deg * 8.4;

% we don't expect anisotropies so let's fix this
norient = 8;
nphase = 2;

meanresp = NaN([numel(stimuli),numel(cpfovs_im)],'like',outvec);
for c = 1:numel(cpfovs_im)
    logstr('\nWorking on cpfovs %.2f (%d of %d)\n',cpfovs_deg(c),...
        c,numel(cpfovs_deg));
    cpname = sprintf('cpfovs%.1f',cpfovs_deg(c));
    [f,gbrs,gaus,sds,indices,info,filters] = applymultiscalegaborfilters(...
        outvec,cpfovs_im(c),-1,1,norient,nphase,.01,2,0);
    if exportimages
        filtind = ceil(prctile(1:size(filters,2),[30 50 60]));
        m = full(max(abs(filters(:))));
        for n = filtind
            outim = intensity2rgb(reshape(filters(:,n+2),...
                [outdim outdim]),gray(1e3),'symmetrical');
            imwrite(outim,fullfile(outdir,sprintf('%s_ex%04d.png',...
                cpname,n)),'PNG');
        end
    end
    % work out location of each filter
    featureLen = size(f,2);
    rawf = squeeze(sqrt(sum(reshape(f,[size(f,1) 2 featureLen/2]).^2,2)));  % images x channels
    % find indices
    nxy = numel(indices{1});
    indarr = false([norient nxy nxy]);
    indarr(:,:,1:floor(nxy/2)) = true;
    leftind{c} = indarr(:);
    assert(numel(leftind{c}) == size(rawf,2));

    % scale by dimensionality
    nfeat(c) = size(rawf,2);
    filterresp{c} = rawf;
    rawf = rawf / sqrt(nfeat(c));
    % build euclidean RDM
    rdv = pdist(rawf);
    outrdms(c) = struct('name',['gwp_' cpname],'RDM',squareform(rdv));
    meanresp(:,c) = mean(rawf,2);
end
leftind = cat(1,leftind{:});
