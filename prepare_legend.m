%Takes a vector of numbers or other type and converts it to a list of strings
function [legendStrs] = prepare_legend(numericalLabels, optionalLabel)
	if nargin() == 1
		optionalLabel = [];
	end
	legendStrs = cellstr(num2str(numericalLabels, [optionalLabel '%-d']));
end
