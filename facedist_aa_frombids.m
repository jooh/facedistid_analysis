%% AA master script for face study.
% Requirements (these should be in startup.m for parallel compute to work)
% * MATLAB R2013A
% * AA 5
% * pilab_aa 
% * study-specific code (facedistid_analysis repo - hey you're in it!)
%
function facedist_aa_frombids(dolocal,reinit)

% could be persistent but then the variable is only saved when the function
% returns successfully. Which is rarely the case.
global facedist_aap

if (exist('reinit','var') && reinit) || isempty(facedist_aap)
    fprintf('initialising AA session\n');
    aa_ver5

    %% DEFINE SPECIFIC PARAMETERS
    %  Default recipe without model
    facedist_aap=aarecipe('aap_parameters_defaults_CBSU.xml','facedist_aa_frombids_tasklist.xml');
    facedist_aap = aas_configforSPM12(facedist_aap);
    facedist_aap.options.NIFTI4D = 1;										% typical value 0
    % Don't spam me
    % facedist_aap.options.email='johan.carlin@mrc-cbu.cam.ac.uk';

    %% STUDY
    % Directory for analysed data - we assume this code lives in bidsroot/code
    rootdir = fileparts(fileparts(mfilename('fullpath')));
    facedist_aap.acq_details.root = fullfile(rootdir,'derivatives');
    facedist_aap.directory_conventions.analysisid = 'aa'; 

    % Add data
    facedist_aap.directory_conventions.rawdatadir = rootdir;
    facedist_aap.acq_details.numdummies = 0;
    facedist_aap.acq_details.input.combinemultiple = 1;
    facedist_aap.options.autoidentifystructural_choosefirst = 1;
    facedist_aap.options.garbagecollection = 0;
    facedist_aap.options.restorepath = 1; 
    facedist_aap = aas_processBIDS(facedist_aap,[],{'main','localizer'});
else
    fprintf('using cached AA session\n');
end

% Modify standard recipe module selection here if you'd like
% Unless your Torque cluster is configured VERY similarly to the CBU's you may
% struggle to run this code in anything other than localsingle
facedist_aap.options.wheretoprocess ='qsub'; 
if exist('dolocal','var') && dolocal
    facedist_aap.options.wheretoprocess = 'localsingle';
end

%% DO ANALYSIS
aa_doprocessing(facedist_aap);
