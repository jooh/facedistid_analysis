% descriptive spatial statistics for the ROIs - mean and st dev for the peak MNI
% coordinates. 
%
% facedist_rois2mni
function facedist_rois2mni

rootdir = fileparts(fileparts(mfilename('fullpath')));
mridir = fullfile(rootdir,'results_fullsample_realign');
roidir = fullfile(mridir,'rois');
normdir = fullfile(mridir,'aamod_norm_noss_00001');
refvoldir = fullfile(mridir,'facedist','aamod_pilab_mask_importbet_00001');

allrois = dirbetter(fullfile(roidir,'*','spm2roi','*.nii'));

coords = struct;

segs = struct;
refvols = struct;

rois = {'EVC','FFA','OFA','TOS','PPA'};
nroi = numel(rois);

spm_jobman('initcfg');
for thisroi = allrois(:)'
    % find roi
    [~,roiname,~] = fileparts(thisroi.name);
    % skip bilateral rois
    if ~any(strfindcell(roiname,{'r_','l_'})) || ~any(strfindcell(roiname,rois))
        continue;
    end
    % figure out subject
    [~,subname] = fileparts(fileparts(fileparts(thisroi.abspath)));
    % find a seg_sn file
    if ~isfield(segs,subname)
        segpath = dirbetter(fullfile(normdir,subname,'structurals','*_seg_sn.mat'));
        assert(numel(segpath)==1);
        segs.(subname) = segpath.abspath;
    end
    subseg = segs.(subname);
    % find a reference volume
    if ~isfield(refvols,subname)
        refpath = fullfile(refvoldir,subname,'pilab','pilab_mask.nii');
        assert(exist(refpath,'file')~=0);
        refvols.(subname) = refpath;
    end
    subref = refvols.(subname);
    % load up the ROI, find its peak
    xyz = spm_read_vols(spm_vol(thisroi.abspath));
    [m,ind] = max(xyz(:));
    vox = ind2subbetter(size(xyz),ind);
    % and find the mni coordinate of that
    mni = vox2spm_mni(vox,subref,subseg);
    % work out if we need to add a new field to the coords output
    if ~isfield(coords,roiname)
        coords.(roiname) = [];
    end
    coords.(roiname)(end+1,1:3) = mni';
end
% and do some summary stat table I suppose
m = structfun(@mean,coords,'uniformoutput',0);
s = structfun(@std,coords,'uniformoutput',0);

datamat = [];
for roi = 1:nroi
    roistr = rois{roi};
    n = 0;
    for pref = {'l_','r_'}
        n = n + 1;
        prefstr = pref{1};
        datamat(n,:,roi) = m.([prefstr roistr]);
        n = n + 1;
        datamat(n,:,roi) = s.([prefstr roistr]);
    end
end
% now we can do
rowlab = {'left','mean'; ...
    '','st dev'; ...
    'right','mean'; ...
    '','st dev'};
outfile = fullfile(roidir,'roistats.csv');
table2csv(datamat,outfile,'precision',1,...
    'collabels',{'','','x','y','z'},'rowlabels',rowlab,'zlabels',...
    facedist_names(rois{:}));
logstr('saved results to %s\n',outfile);
