% contrasts plugin input for aamod_pilab_rsa_rfx. Used to generate EDI
% contrasts for across view effects in full view condition and slope
% difference tests for full vs no direction tuning models.
%
% contrasts = facedist_rsa_rfx_contrasts_meanslopes(connames)
function contrasts = facedist_rsa_rfx_contrasts_meanslopes(names)

directions = mat2strcell([11 17,25,28],'_direction_%.0f');
ecccons = {'sub-sub','typical-typical','caricature-caricature'};
contrasts = emptystruct('name','conplus','conminus','tail');
withinind = strfindcell(names,'view_within');
if isempty(withinind)
    % collapsed view case
    viewprefixes = struct('prefix','','ind',1:numel(names));
else
    % multiple view case
    viewprefixes = struct('prefix',{'view_within_','view_across_'},...
        'ind',{withinind,strfindcell(names,'view_across')});

    % EDI contrast for each eccentricity level
    for g = ecccons
        thisg = g{1};
        % indices for delta dir at this caricature level
        conind = cellfun(@(x)strfindcell(names,['view_across_' thisg x]),...
            directions,'uniformoutput',true)';
        contrasts(end+1).name = ['view_across_edi_' thisg];
        contrasts(end).conplus = names(conind(2:end));
        contrasts(end).conminus = names(conind(1))';
        contrasts(end).tail = 'right';
    end
    % collapse into big EDI
    contrasts(end+1) = struct('name','view_across_edi','conplus',...
        {[contrasts.conplus]},'conminus',{[contrasts.conminus]},'tail',...
        'right');
end

for v = 1:numel(viewprefixes)
    pre = viewprefixes(v);
    fullind = strfindcell(names,[pre.prefix 'slope_full']);
    nodirind = strfindcell(names,[pre.prefix ...
        'slope_no_direction_tuning']);
    assert(numel(fullind)==1 && numel(nodirind)==1)
    contrasts(end+1).name = [pre.prefix ...
        'slope_full_o_no_direction_tuning'];
    contrasts(end).conplus = names(fullind);
    contrasts(end).conminus = names(nodirind);
    contrasts(end).tail = 'both';

    % also collapse the slopes from each absecc
    % and for both X scales
    for suff = {'','_angscale'}
        suffstr = suff{1};
        contrasts(end+1).name = [pre.prefix 'mean_dirslope' suffstr];
        contrasts(end).con = {};
        for g = ecccons
            thisg = g{1};
            ind = strcmp(names,[pre.prefix 'slope_' thisg suffstr]);
            assert(sum(ind)==1);
            contrasts(end).conplus(end+1) = names(ind);
        end
        contrasts(end).tail = 'right';
    end

    for d = 1:numel(directions)
        if d==1 && ~strcmp(pre.prefix,'view_across_')
            continue
        end
        contrasts(end+1).name = [pre.prefix 'mean' directions{d}];
        n = cellfun(@(thisecc)strfindcell(names,...
            [pre.prefix thisecc directions{d}]),ecccons,...
            'uniformoutput',1);
        contrasts(end).conplus = names(n);
        contrasts(end).tail = 'right';
    end
end

% slope by view interactions and means
if numel(viewprefixes)>1
    withinind = strfindcell(names,'view_within_slope_');
    acrossind = strfindcell(names,'view_across_slope_');
    % I think these should match up but let's confirm
    allok = arrayfun(@(x)isequal(names{withinind(x)}(13:end),...
        names{acrossind(x)}(13:end)),1:numel(withinind));
    assert(all(allok),'mismatched contrasts between views');
    % so assuming that worked we can just do something like...
    for n = 1:numel(withinind)
        contrasts(end+1).name = ['contrast_across-within_' ...
            names{withinind(n)}(13:end)];
        contrasts(end).conplus = names{acrossind(n)};
        contrasts(end).conminus = names{withinind(n)};
        contrasts(end).tail = 'both';
    end
    % and also
    for n = 1:numel(withinind)
        contrasts(end+1).name = ['mean_across_and_within_' ...
            names{withinind(n)}(13:end)];
        contrasts(end).conplus = {names{acrossind(n)},...
            names{withinind(n)}};
        contrasts(end).tail = 'right';
    end
end
