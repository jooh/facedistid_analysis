% Main batch analysis of idanddistinct study
% path to root of current dir
rootdir = fileparts(fileparts(mfilename('fullpath')));
% in there we expect to put processed data
[~,wd] = mkdirifneeded(fullfile(rootdir,'results_mri'));
% and we also expect to find the preprocessed data
exportpath = fullfile(rootdir,'aamod_pilab_sharedata_00001',...
    'pilab_export','pilab_exportspec.mat');
assert(exist(exportpath,'file')~=0,'could not find data in location %s',...
    exportpath);

if ~any(strcmp(version,'R2013a'))
    % any non R2013a versions may behave differently
    warning(['This software was developed for Matlab version R2013a ' ...
        '(8.1.0.604). Your version is %s. Reproducibiliy is not ' ...
        'guaranteed.'],version);
end

% necessary?
global defaults
spm('Defaults','FMRI')

clear aap
aap = aarecipe('aap_defaults_idanddistinct.xml',...
    'aap_tasklist_idanddistinct_anonymous.xml');

aap.acq_details.root = wd;

% I think the session numbers are now quite different. No. Everything
% should be the same order (so 2 main and 2 localiser runs alternating).
nloc = 8;
nexp = 8;
ntot = nloc+nexp;
aap = aas_addsubject(aap,'be',[1:ntot]);
aap = aas_addsubject(aap,'cd',[1:ntot]);
aap = aas_addsubject(aap,'cr',[1:ntot]);
aap = aas_addsubject(aap,'mc',[1:ntot]);
aap = aas_addsubject(aap,'mg',[1:ntot]);
aap = aas_addsubject(aap,'rc',[1:ntot]);
aap = aas_addsubject(aap,'sl',[1:ntot]);
aap = aas_addsubject(aap,'td',[1:ntot]);
aap = aas_addsubject(aap,'tg',[1:ntot]);
aap = aas_addsubject(aap,'lt',[1:ntot]);

% TODO: hm. Not sure the epi and dicom headers made it through the import
% process.
%
% so now we need to export the preprocessed stuff. The best way to do this
% would be to have yet another AA module in the preprocessing stream that
% requires all the files as input. Then we can just copy the module
% directory to obtain everything.
% So the export module will need:
% epiBETmask - actually, pilab_mask
% normalisation_seg_sn - DONE
% epi - DEBUG
% epi_dicom_header - DEBUG
% structural_group
%
% the export module should also save a mat file with all these stream names
% and paths so we can easily load this mat up and run the import procedure
% below. Probably we will write a little aas_importdata function for this.

% so now we should be able to do something like
aap = aas_pilab_importdata(aap,exportpath);
% and this function will work out absolute file paths relative to the
% location of the exportspec mat. Nice and neat.

% sessions are in the same order as in original
aap = aas_addsession(aap,'main-1');
aap = aas_addsession(aap,'main-2');
aap = aas_addsession(aap,'localiser-1');
aap = aas_addsession(aap,'localiser-2');
aap = aas_addsession(aap,'main-3');
aap = aas_addsession(aap,'main-4');
aap = aas_addsession(aap,'localiser-3');
aap = aas_addsession(aap,'localiser-4');
aap = aas_addsession(aap,'main-5');
aap = aas_addsession(aap,'main-6');
aap = aas_addsession(aap,'localiser-5');
aap = aas_addsession(aap,'localiser-6');
aap = aas_addsession(aap,'main-7');
aap = aas_addsession(aap,'main-8');
aap = aas_addsession(aap,'localiser-7');
aap = aas_addsession(aap,'localiser-8');

aa_doprocessing(aap);
