function con = facedist_rsaglm_contrasts(connames)

% either way we want the pairwise tests
con = roidata_allpairwisecontrasts(connames);
withinind = strfindcell(connames,'view_within_');
if ~any(withinind)
    % no need for any of this
    return
end

% now in this case we additionally want to average the view effects
n = numel(con);
for thisind = withinind(:)'
    withinstr = connames{thisind};
    acrossstr = strrep(withinstr,'view_within_','view_across_');
    assert(sum(strcmp(connames,acrossstr))==1)
    con(end+1) = struct('name',strrep(withinstr,'view_within_',''),...
        'conplus',{{withinstr,acrossstr}},'conminus',{{}},'tail','right');
end
% and contrast those too
con = [con roidata_allpairwisecontrasts({con(n+1:end).name})];
