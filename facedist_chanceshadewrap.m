% wrapper for default chanceshade for facedistid.
%
% [h,t] = facedist_chanceshadewrap(x,y,textx)
function [h,t] = facedist_chanceshadewrap(x,y,textx)

xs = size(x);
ys = size(y);
assert(isequal(xs,ys));
x = reshape(x,[xs(1) prod(xs(2:end))]);
y = reshape(y,[xs(1) prod(xs(2:end))]);

[h,t] = chanceshade(x,y,'inftype','parametric',...
    'correction','uncorrected',...
    'textx',textx,...
    'xvalues','vary');
