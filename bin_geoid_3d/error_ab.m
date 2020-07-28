function error_ab(str)
%ERROR_AB Displays the error message and exits or stops.
%
% ERROR_AB('text')
% Displays the error message 'text'; then, in case of -nodesktop session, 
% Matlab exits to avoid infinite display of error messages, 
% otherwise the script is only stopped.

% Ales Bezdek, bezdek@asu.cas.cz, 1/2015

if ~usejava('Desktop')
   disp(str);
   exit; 
else
   error(str);
end


