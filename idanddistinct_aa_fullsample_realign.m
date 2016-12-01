% Main batch analysis of idanddistinct study
% path to current dir
wd = fullfile(fileparts(fileparts(mfilename('fullpath'))),...
    'results_mri_preprocess');
mkdirifneeded(wd);
%addpath(genpath(cwd));
% (we assume you already have aa4 and aa_custom on path)

% necessary?
global defaults
spm('Defaults','FMRI')

clear aap
aap = aarecipe('aap_defaults_idanddistinct.xml',...
    'aap_tasklist_idanddistinct_preprocess.xml');

aap.acq_details.root = wd;

aap = aas_addsubject(aap,'be',[105:108 205:208 305:308 405:408],[]);
% skipping bad second scan (now deleted compeltely to be safe)
aap = aas_addsubject(aap,'cd',[105:108 305:308 406:409 505:508],[]);
aap = aas_addsubject(aap,'cr',[105:108 205:208 305:308 405:408],[]);
% skipping first session (half was with wrong stims)
% (check that we don't end up importing the T1 / fieldmaps from first
% session now)
% (to be sure, just delete first series completely)
aap = aas_addsubject(aap,'mc',[205:208 305:308 405:408 505:508]);%,[104:110]);
aap = aas_addsubject(aap,'mg',[105:108 205:208 305:308 405:408],[]);
aap = aas_addsubject(aap,'rc',[105:108 205:208 305:308 405:408],[]);
aap = aas_addsubject(aap,'sl',[105:108 205:208 305:308 405:408],[]);
aap = aas_addsubject(aap,'td',[105:108 205:208 305:308 405:408],[]);
% skipping 4th (TD stims) - now deleted to be safe
aap = aas_addsubject(aap,'tg',[105:108 205:208 305:308 505:508],[109:110]);
aap = aas_addsubject(aap,'lt',[106:109 205:208 305:308 405:408],[105]);

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
