function gca_findpeak(x, y)
%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, June 2019
% xan@cshl.edu
% Version: 1.0
%*---------------------------------------------------------------------*
xrange = get(gca, 'XLim');
xrange_peak = x(x >= xrange(1) & x <= xrange(2));
[ypeak, Index] = max(y(x >= xrange(1) & x <= xrange(2)));
xpeak = xrange_peak(Index);
hold on;
plot(xpeak, ypeak, 'ob');
text(xpeak, ypeak, ['\leftarrow' num2str(ypeak, '%.2f')], 'Color', 'b');