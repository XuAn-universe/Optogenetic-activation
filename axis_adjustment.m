%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Sep 2021
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%%
xl1 = xlim;
yl1 = ylim;
zl1 = zlim;

%%
xl2 = xlim;
yl2 = ylim;
zl2 = zlim;

%%
xl(1) = min(xl1(1), xl2(1));
xl(2) = max(xl1(2), xl2(2));
yl(1) = min(yl1(1), yl2(1));
yl(2) = max(yl1(2), yl2(2));
zl(1) = min(zl1(1), zl2(1));
zl(2) = max(zl1(2), zl2(2));

%%
xlim(xl);
ylim(yl);
zlim(zl);

%%
view = get(gca, 'View');

%%
set(gca, 'View', view);

%%
ylim([0 140]);

%%
pos = get(gcf, 'Position');

%%
set(gcf, 'Position', pos);