% fevaled in aamod_pilab_importrsapredictors
%
% rdms = facedist_rsapredictors(subname)
function rdms = facedist_rsapredictors(subname,varargin)

oldpath = path;

% function to add psychtoolbox functions to the path - you may not need this,
% but on the CBU cluster Matlab figure rendering breaks when this is on the path
% so I tend to leave it off...
try
    start_psychtoolbox;
catch err
    assert(strcmp(err.identifier,'MATLAB:undefinedfunction'),'unknown error');
end

rootdir = fileparts(fileparts(mfilename('fullpath')));
subname = subname{1}(1:6);
behaviourdir = fullfile(rootdir,'bids','sourcedata','misc',subname);
assert(exist(behaviourdir,'dir')~=0,'could not find behaviourdir: %s',behaviourdir);

% process stimulusspace - basic predictions across full matrix
stimpath = fullfile(behaviourdir,'stimuli_mri.mat');
ss = loadbetter(stimpath);

rdmfull = struct('name','full','RDM',...
    ss.rdmbyattribute('shapecoef','euclidean'));

% to be replaced with RDMs that preserve the original euclidean
% distances of the stimulus space.
rdmang = struct('name','no_eccentricity_discrimination','RDM',...
    ss.rdmbyattribute('shapecoef','angdist'));
rdmecc = struct('name','no_direction_discrimination','RDM',...
    ss.rdmbyattribute('shapecoef','eccdist'));

% handle absecc - this gets a little hacky but...
resortind = facedist_resortind(24);
[~,reverseind] = sort(resortind);

eucrdm = rdmfull.RDM(resortind,resortind);
angrdm = rdmang.RDM(resortind,resortind);
rdmecc.RDM = rdmecc.RDM(resortind,resortind);
eccfir = rdm2fir(rdmecc);
angfir = rdm2fir(angrdm);
% get rid of the 0 change case since we model this separately below
eccfir(1) = [];
% make a mask for each absolute eccentricity level (this is for subs but we
% can circshift this around to fit)
abseccfir(1) = struct('name','absecc_f1','RDM',...
    [ones(4,4) zeros(4,8); zeros(8,12)]);
% upcast to full 24 x 24 rdm
abseccfir(1).RDM = repmat(abseccfir(1).RDM,[2 2 1]);
abseccfir(1).RDM(diagind(size(abseccfir(1).RDM))) = 0;
abseccfir(2) = struct('name','absecc_f2','RDM',...
    circshift(abseccfir(1).RDM,[4 4]));
abseccfir(3) = struct('name','absecc_f3','RDM',...
    circshift(abseccfir(1).RDM,[8 8]));
alleccfir = [abseccfir eccfir];

% now we can construct the new eccdist model
neweccmodel = RSA(alleccfir,eucrdm);
rdmeccnew.name = 'no_direction_tuning';
rdmeccnew.RDM = asrdmmat(neweccmodel.predictY);

% the old eccdist model
oldeccmodel = RSA(eccfir,eucrdm);
rdmeccold.name = 'no_direction_discrimination';
rdmeccold.RDM = asrdmmat(oldeccmodel.predictY);

% and the old angdist model
angmodel = RSA(angfir,eucrdm);
rdmangold.name = 'no_eccentricity_discrimination';
rdmangold.RDM = asrdmmat(angmodel.predictY);

% (by the way - the old eccentricity model correlates similarly with full
% as the new (tau=.38 for new, tau=.42 for old). This is important.

% combine the rdms
newrdms = [rdmeccnew rdmeccold rdmangold];
% reverse the sorting
for n = 1:numel(newrdms)
    newrdms(n).RDM = newrdms(n).RDM(reverseind,reverseind);
end
% and now we have
rdms = [rdmfull newrdms];

rdms(end+1) = struct('name','vid_corrdist','RDM',...
    ss.rdmbyattribute('videoview','corr'));

clear ss

% indices corresponding to basic predictors
baseinds = 1:4;

% add rating task
pairpath = fullfile(behaviourdir,'data_perceptual_judgment_task.mat');
subdata = loadbetter(pairpath);
assert(length(subdata)==16 || strcmp(subname,'sub-01'),...
    'incorrect number of facepairs sessions')
rdsum = zeros(12,12);
for sess = 1:length(subdata)
    rdsum = rdsum + subdata(sess).res.rdm;
end
rdms(end+1) = struct('name','facepairs','RDM',rdsum);
clear subdata

% binary viewpoint RDM
rdms(end+1) = struct('name','view','RDM',squareform(pdist(...
    [ones(12,1); zeros(12,1)],'meandist')));

% view-specific variants
diagonal = logical(eye(24));
blockf = false(12,12);
blockt = true(12,12);
% parts to mask for right-specific analysis
notright = [blockf blockt; blockt blockt];
notright(diagonal) = false;
viewvariant.name = 'right';
viewvariant.nanmask = notright;
% parts to mask for left-specific
notleft = [blockt blockt; blockt blockf];
notleft(diagonal) = false;
viewvariant(end+1).name = 'left';
viewvariant(end).nanmask = notleft;

% within view
viewvariant(end+1).name = 'within';
viewvariant(end).nanmask = zerodiagonal([blockf blockt; blockt blockf]);

% across view
notacross = [blockt, blockf; blockf blockt];
notacross(diagonal) = false;
viewvariant(end+1).name = 'across';
viewvariant(end).nanmask = notacross;

for r = baseinds
    for v = 1:length(viewvariant)
        rd = rdms(r).RDM;
        rd(viewvariant(v).nanmask) = NaN;
        rdms(end+1) = struct('name',sprintf('view_%s_%s',...
            viewvariant(v).name,rdms(r).name),'RDM',rd);
    end
end

nrelevant = length(rdms);

% sort 
leftind = strfindcell({rdms.name},'view_left');
rightind = strfindcell({rdms.name},'view_right');
acrossind = strfindcell({rdms.name},'view_across');
others = setdiff(1:nrelevant,[acrossind,leftind,rightind]);
% rebuild grouped by view
rdms = [rdms(others) rdms(acrossind) rdms(leftind) rdms(rightind)];

% upcast if necessary
for n = 1:length(rdms)
    if size(rdms(n).RDM,1)==12
        rdms(n).RDM = repmat(rdms(n).RDM,[2 2]);
    end
    assert(numel(unique(rdms(n).RDM(~isnan(rdms(n).RDM))))>1,...
        'invariant RDM: %s',rdms(n).name);
end

path(oldpath);
