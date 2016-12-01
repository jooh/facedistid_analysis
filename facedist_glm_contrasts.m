% contrasts plugin input for roidata_glm.
%
% contrasts = facedist_glm_contrasts(connames)
function contrasts = facedist_glm_contrasts(names)

% so we need a programmatic way to define contrasts
par = facedist_plotpar;
eccs = par.absecclabels;
dirs = par.dirlabels;
views = {'left','right'};
% so now upcast each of these...
lab.views = views(repinds(1,2,12));
lab.eccs = repmat(eccs(repinds(1,3,4)),[1 2]);
lab.dirs = repmat(dirs,[1 6]);

% so first of all, basic view contrast.
contrasts(1).name = 'contrast_left_o_right';
contrasts(1).convec = zeros(1,24);
contrasts(1).convec(strcmp(lab.views,'left')) = 1;
contrasts(1).convec(strcmp(lab.views,'right')) = -1;
contrasts(1).tail = 'both';
% used later in facedist_dirplot / facedist_abseccplot hack
contrasts(1).cons = [];

% distances for each dir - actually no, just linear here since it's not
% distances.
dirdist = zeromean([1:4]);

viewprefs = struct('skiptarg',{'','right','left'},...
    'prefix',{'all_','view_left_','view_right_'});

for thisv = viewprefs(:)'
    vind = ~strcmp(lab.views,thisv.skiptarg);
    subind = vind & strcmp(lab.eccs,'sub');
    typind = vind & strcmp(lab.eccs,'typical');
    carind = vind & strcmp(lab.eccs,'caricature');
    % means - ecc
    contrasts(end+1) = struct('name',[thisv.prefix 'mean_ecc_sub'],...
        'convec',zeros(1,24),'tail','right','cons',[]);
    contrasts(end).convec(subind) = 1;
    contrasts(end+1) = struct('name',[thisv.prefix 'mean_ecc_typical'],...
        'convec',zeros(1,24),'tail','right','cons',[]);
    contrasts(end).convec(typind) = 1;
    contrasts(end+1) = struct('name',[thisv.prefix 'mean_ecc_caricature'],...
        'convec',zeros(1,24),'tail','right','cons',[]);
    contrasts(end).convec(carind) = 1;
    % main eccentricity slope contrast
    contrasts(end+1).name = [thisv.prefix 'slope_absecc'];
    contrasts(end).tail = 'right';
    contrasts(end).convec = zeros(1,24);
    contrasts(end).convec(subind) = -1;
    contrasts(end).convec(carind) = 1;
    % plug in names for each mean so plot function can find them later
    contrasts(end).cons = {contrasts(end-3:end-1).name};
    % means - dir
    for d = 1:4
        dirind{d} = vind & strcmp(lab.dirs,dirs{d});
        contrasts(end+1) = struct(...
            'name',[thisv.prefix 'mean_dir_' dirs{d}],...
            'convec',zeros(1,24),'tail','right','cons',[]);
        contrasts(end).convec(dirind{d}) = 1;
    end
    % main dir slope
    contrasts(end+1).name = [thisv.prefix 'slope_no_eccentricity_tuning'];
    contrasts(end).convec = zeros(1,24);
    contrasts(end).tail = 'right';
    for d = 1:4
        contrasts(end).convec(dirind{d}) = dirdist(d);
    end
    % insert names for plotting
    contrasts(end).cons = {contrasts(end-4:end-1).name};
    alleccdirind = numel(contrasts);
    % also eccentricity-specific means and slope fits
    for thise = eccs(:)'
        estr = ['_' thise{1} '-' thise{1}];
        % means for this dir + ecc combo
        for d = 1:4
            contrasts(end+1) = contrasts(alleccdirind-5+d);
            contrasts(end).name = [thisv.prefix 'mean_dir_' dirs{d} estr];
            contrasts(end).convec(~strcmp(lab.eccs,thise{1})) = 0;
        end
        % dir slope for this ecc
        contrasts(end+1) = contrasts(alleccdirind);
        contrasts(end).name = [thisv.prefix 'slope' estr];
        % zero out non-matching eccentricies
        contrasts(end).convec(~strcmp(lab.eccs,thise{1})) = 0;
        % add cons for plotting
        contrasts(end).cons = {contrasts(end-4:end-1).name};
    end
end

switch numel(names)
    case 12
        % need to collapse. First drop the view-specific contrasts.
        contrasts(strfindcell({contrasts.name},'left')) = [];
        contrasts(strfindcell({contrasts.name},'right')) = [];
        % then shorten the vectors for the remainder
        for c = 1:numel(contrasts)
            contrasts(c).convec(13:24) = [];
        end
    case 24
        % all good
    otherwise
        error('names array has unsupported length: %d',numel(names));
end

% normalise contrasts so we can compare mean parameter estimates
for c = 1:numel(contrasts)
    posind = contrasts(c).convec>0;
    contrasts(c).convec(posind) = contrasts(c).convec(posind) ./ sum(...
        contrasts(c).convec(posind));
    negind = contrasts(c).convec<0;
    contrasts(c).convec(negind) = contrasts(c).convec(negind) ./ sum(...
        contrasts(c).convec(negind)) .* -1;
end
