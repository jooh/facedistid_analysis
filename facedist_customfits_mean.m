% set up slope fits for the RFX analysis of MeanRSA estimates. Used as
% plugin customfits input for roidata_rfx.
%
% customfit = facedist_customfits_mean(contrasts)
function customfit = facedist_customfits_mean(varargin)

% extract the mean face space distance for each condition
names = varargin{1};
% extract distances
for n = 1:numel(names)
    nstr = names{n};
    if any(strfindcell(nstr,{'slope_','_edi','dirslope','_angscale',...
            'mean_direction_'}))
        distances(n) = NaN;
        continue;
    end
    ind = strfind(nstr,'_m');
    assert(~isempty(ind));
    ind = ind(end);
    distances(n) = str2double(nstr(ind+2:end));
end

% view coding
withinind = strfindcell(names,'view_within');
if isempty(withinind)
    % collapsed view case
    viewprefixes = struct('prefix','','ind',1:numel(names));
else
    viewprefixes = struct('prefix',{'view_within_','view_across_'},...
        'ind',{withinind,strfindcell(names,'view_across')});
end

% wait a minute. 
% I re-introduced 11 below so that we could get all the dir(0) cases in for
% off block diagonals. But it's had the unintended consequence of breaking
% the across view slope code (where dir(0) was excluded). So I think the
% answer is to go back to excluding it for the across view first 3 ecccons.

% nb for the across view case we could also try to fit a slope to the
% within-view case (11)
xstrings = mat2strcell([11,17,25,28],'_direction_%.0f');
ecccons = {'sub-sub','typical-typical','caricature-caricature',...
    'sub-typical','typical-caricature','sub-caricature'};
customfit = emptystruct('name','funhand','cons','x','tail');
for v = 1:numel(viewprefixes)
    pre = viewprefixes(v);
    viewnames = names(pre.ind);
    viewdistances = distances(pre.ind);

    customfit(end+1).name = [pre.prefix 'slope_absecc'];
    customfit(end).funhand = 'fitslope';
    coninds = cellfun(@(x)strfindcell(viewnames,[x '_m']),ecccons(1:3),...
        'uniformoutput',true);
    customfit(end).cons = viewnames(coninds);
    customfit(end).x = viewdistances(coninds)';
    customfit(end).tail = 'right';

    % make direction slope fit for each eccentricity level
    for g = ecccons
        thisg = g{1};
        customfit(end+1).name = [pre.prefix 'slope_' thisg];
        customfit(end).funhand = 'fitslope';
        % NB cell2mat on non-uniform to support empty outputs (no dir(0)
        % case for block diagonals)
        conind = cell2mat(...
            cellfun(@(x)strfindcell(viewnames,[thisg x]),xstrings,...
            'uniformoutput',false)');
        if any(strfindcell(ecccons(1:3),thisg)) && numel(conind)==4
            % drop the across view block diagonal
            conind(1) = [];
        end
        customfit(end).cons = viewnames(conind);
        customfit(end).x = viewdistances(conind)';
        customfit(end).tail = 'right';
    end
    abseccind = (numel(customfit) -3) -[2:-1:0];

    % concatenate the conditions - this is the full euclidean model
    customfit(end+1).name = [pre.prefix 'slope_full'];
    customfit(end).funhand = 'fitslope';
    customfit(end).cons = vertcat(customfit(abseccind).cons);
    customfit(end).x = vertcat(customfit(abseccind).x);
    customfit(end).tail = 'right';

    % average the directions at each eccentricity step - this is the no
    % direction tuning model
    customfit(end+1) = customfit(end);
    customfit(end).name = [pre.prefix 'slope_no_direction_tuning'];
    customfit(end).x = [repmat(mean(customfit(end).x(1:3)),[3 1]); ...
        repmat(mean(customfit(end).x(4:6)),[3 1]); ...
        repmat(mean(customfit(end).x(7:9)),[3 1])];

    % average the eccentricity levels for each direction - this is the no
    % eccentricity tuning model (and is equivalent to averaging the
    % eccentricities and fitting a slope to the resulting delta dirs)
    customfit(end+1) = customfit(end-1);
    customfit(end).name = [pre.prefix 'slope_no_eccentricity_tuning'];
    mx = [mean(customfit(end).x([1 4 7])); ...
        mean(customfit(end).x([2 5 8])); ...
        mean(customfit(end).x([3 6 9]))];
    customfit(end).x = repmat(mx,[3 1]);
    
    % interaction term
    customfit(end+1).name = [pre.prefix 'slope_ofslopes'];
    customfit(end).funhand = 'fitslope';
    customfit(end).cons = {customfit(abseccind).name};
    customfit(end).tail = 'right';
    % not sure if there is a way to produce meaningful units here.
    customfit(end).x = [1:3]';

end

% second pass with uniform X scale
for c = 1:numel(customfit)
    if any(strfindcell(customfit(c).name,{'slope_absecc','slope_full','slope_no_direction_tuning'}))
        % these don't make sense here
        continue
    end
    customfit(end+1) = customfit(c);
    customfit(end).name = [customfit(end).name '_angscale'];
    if any(strfind(customfit(c).name,'slope_ofslopes'))
        % here we need to do make sure we fit the revised slopes instead
        customfit(end).cons = cellfun(@(thisc)[thisc '_angscale'],...
            customfit(c).cons,'uniformoutput',0);
    else
        % for everything else, we make uniform
        ux = unique(customfit(end).x);
        customfit(end).x = arrayfun(@(thisx)find(ux==thisx),...
            customfit(end).x);
    end
end
