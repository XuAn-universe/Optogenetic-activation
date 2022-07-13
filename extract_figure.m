function extract_figure(haxis, ~)
%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, June 2019
% xan@cshl.edu
% Version: 1.0
%*---------------------------------------------------------------------*
hf = figure('Color', [1 1 1]);
try
    copyobj(haxis, hf);
    set(gca, 'Position', [0.15 0.13 0.7 0.8]);
catch
    set(haxis, 'ButtonDownFcn', []);
    close(hf);
end