% collapse the left and right viewpoint regressors in the SPM first level model
% (plugin for aamod_firstlevel_model_collapsepredictors).
%
% collapsestruct = facedist_collapsespm(SPM)
function collapsestruct = facedist_collapsespm(SPM)

% make sure the regressor names are consistent
allnames = arrayfun(@(thissess)[thissess(1).U.name],SPM.Sess,'uniformoutput',0);
assert(isequal(allnames{1},allnames{:}));
allnames = allnames{1};

collapsestruct = arrayfun(@(c)struct('name',allnames{c},'targetnames',...
    {allnames([c c+12])}),1:12,'uniformoutput',0);
collapsestruct = horzcat(collapsestruct{:});
