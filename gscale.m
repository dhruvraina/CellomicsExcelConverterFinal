%Description: Scales Graph Bins
%Author: Dhruv Raina
%Last Edit: 051215
%Takes in X values or Y Values, puts them into bins, and returns bin width
%ToDo:
%1. find the difference between the 98% and the 95% and compare it with the
%top 1percentile. If there is too big a spread, trim edge values.

function varargout = gscale(inVals, binno, matType)

switch matType
    case('cell')
        inValconcat = cell2mat(cellfun(@(x) max(x), inVals, 'UniformOutput', 0));
        inValconcat(inValconcat==max(inValconcat(:))) = 0;
        inValsubmax = max(inValconcat(:));
    case('double')
        inValconcat = inVals;
        inValconcat(inValconcat==max(inValconcat(:))) = 0;
        inValsubmax = max(inValconcat(:));
end
binwidth = floor((inValsubmax/binno)*10)/10;
binmax = inValsubmax;
varargout{1} = binwidth;
varargout{2} = binmax;
