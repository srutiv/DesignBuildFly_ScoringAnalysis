function [alphas, CLs, CDs] = parsePolar(path)
%path = './polar/mh121_flap_200000.txt';
CLCD = 0; CDp = 0;
file = textread(path, '%s', 'delimiter', '\n','whitespace', ' ');


% Find where the data begins
i = 12;

% Store all results in arrays
[alphas, CLs, CDs, CDps] = getResults(file,i);
% CLCDs = CLs./CDs;
% 
% idx = find(CDps == min(CDps(:)), 1, 'last');
% fprintf('CDpmin: %f   CD: %f   Angle: %f\n', CDps(idx), CDs(idx), alphas(idx));
% idx = find(CLCDs == max(CLCDs(:)), 1, 'last');
% fprintf('CL/CD max: %f   Angle: %f\n', CLCDs(idx), alphas(idx));
% idx = find(CLs == max(CLs(:)), 1, 'last'); % index of max CL
% fprintf('CLmax: %f    CD: %f   Angle (stall angle): %f\n', CLs(idx), CDs(idx), alphas(idx));
% end
end





% Returns the alpha, CL, and CD results from a specified polar file and
% index i
function [alphas, CLs, CDs, CDps] = getResults(file,i)
alphas = []; CLs = []; CDs = []; CDps = [];
if(i ~= length(file) && ~isempty(file(i+1)))
    j = i+1;
    while (j<=length(file) && ~isempty(file(j)))
        results = textscan(char(file(j)), '%f');
        alphas(j-i) = results{1}(1);
        CLs(j-i) = results{1}(2);
        CDs(j-i) = results{1}(3);
        CDps(j-i) = results{1}(4);
        j = j+1;
    end
end
end

