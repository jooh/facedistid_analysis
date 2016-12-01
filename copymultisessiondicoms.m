% This function adds data for a given subject name to the raw directory
% varargs are targets (cell of sequence names), outroot (outputs get
% written to outroot/subname).
% Uses a parfor loop to make things less excruciatingly slow.
% varargs:
% targets: cell array of target series
% outroot: directory to put data into
%
% function copymultisessiondicoms(subname,cbucode,session,[varargin])
function copymultisessiondicoms(subname,cbucode,session,varargin)

getArgs(varargin,{'targets',{'MEMPRAGE_4E_1mm_RMS','3D_EPI_2mm',...
    '3D_EPI_2mm_localiser','FieldMapping_32chn','3D_EPI_2mm_B940',...
    '3D_EPI_2mm_localiser_B380'},...
    'outroot','/Users/jc01/research/idanddistinct/raw/dicoms'});

mridata = '/mridata/cbu';

assert(all(session<10),...
    'this method assumes < 10 sessions and < 100 series per session');

subdir = fullfile(outroot,subname);
if ~exist(subdir)
    mkdir(subdir);
end

if iscell(cbucode)
    nsess = length(cbucode);
    assert(nsess == length(session),'n cbucodes must match n session')
else
    nsess = 1;
    cbucode = {cbucode};
    assert(length(session)==1,'n cbucodes must match n session')
end

for c = 1:nsess
    fprintf('---- session %d of %d...\n',c,nsess);
    cb = cbucode{c};
    sess = session(c);
    subonmri = dir(fullfile(mridata,[cb '*']));
    assert(length(subonmri)==1,'failed to find a single mridata dir');
    subonmri = fullfile(mridata,subonmri.name);
    subonmri_sub = dir(fullfile(subonmri,'20*'));
    if length(subonmri_sub) < 1
        error('failed to find any mridata subdir');
    elseif length(subonmri_sub) > 1
        % pick last
        [s,inds] = sort({subonmri_sub.name});
        subonmri_sub = subonmri_sub(inds(end));
    end
    subonmri_sub = fullfile(subonmri,subonmri_sub.name);
    candidates = dir(fullfile(subonmri_sub,'Series*'));
    candidates = {candidates.name};
    for t = 1:length(targets)
        targ = targets{t};
        hits = strfindcell(candidates,targets{t});
        tlen = length(targ);
        for h = 1:length(hits)
            hh = hits(h);
            % confirm that there is an exact match
            [x,cand] = fileparts(candidates{hh});
            candend = cand(end-tlen+1:end);
            if strcmp(candend,targ)
                % find a new series number
                oldnum = str2num(candidates{hh}(8:10));
                assert(oldnum < 100,'this method assumes <100 series');
                newnum = (sess*100) + oldnum;
                % setup an output dir with correct series num
                outdir = fullfile(subdir,sprintf('Series_%03d_%s',...
                    newnum,targ));
                assert(~exist(outdir,'dir'),'duplicate directory detected')
                mkdir(outdir);
                indir = fullfile(subonmri_sub,candidates{hh});
                fprintf('input dir: %s\n',indir);
                fprintf('output dir: %s\n',outdir);
                fprintf('copying series %d as %d\n',oldnum,newnum);
                % update series num in each dicom
                dicoms = dir(fullfile(indir,'*.dcm'));
                parfor d = 1:length(dicoms)
                    metadata = dicominfo(fullfile(subonmri_sub,...
                        candidates{hh},dicoms(d).name));
                    X = dicomread(metadata);
                    assert(metadata.SeriesNumber==oldnum,...
                        'SeriesNumber does not match dir name');
                    metadata.SeriesNumber = newnum;
                    outfile = fullfile(outdir,dicoms(d).name);
                    % write to outdir with correct series num (createmode
                    % copy and WritePrivate are needed to force write of
                    % non-standard header fields (ie Private_X etc)
                    dicomwrite(X,outfile,metadata,'CreateMode','copy',...
                        'WritePrivate',true);
                end
            end
        end
    end
end
fprintf('Done.\n')
