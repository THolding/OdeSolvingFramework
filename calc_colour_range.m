function [colours] = calc_colour_range (values, startColour, endColour)
	%values = [2 3 4]; startColour = [0 0 1]; endColour = [1 1 0];
	minVal = min(values);
	maxVal = max(values);
	valRange = maxVal-minVal;

	colRange = endColour - startColour;

	colours = zeros(length(values), 3);
	
	%Calculate and add each oclour
	for ii = 1:length(values)
		%Calculate colour modifier
		modifier = (values(ii)-minVal)/valRange;
		curColour = startColour + (colRange * modifier);
		colours(ii, :) = curColour;
    end
end
