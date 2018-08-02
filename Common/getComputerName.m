function compName = getComputerName()
% Returns the name of the computer (hostname)
% output string is converted to lower case
% By Noam Omer Dec17

[re, compName] = system('hostname');   
if re ~= 0,
	if ispc
		compName = getenv('COMPUTERNAME');
	else
		compName = getenv('HOSTNAME');
	end;
end;

compName = strtrim(lower(compName));

