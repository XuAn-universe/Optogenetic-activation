function segTH = fit_gaussian2dist(data, nbins, xlabel_text)
[ncount, edges] = histcounts(data, nbins);
edges = movmean_xu(edges, 2);
% fit with a Gaussian function using fminsearch 
Starting_p = [max(ncount), edges(find(ncount == max(ncount), 1)), edges(find(ncount == max(ncount), 1))]; % A, Mu, SD
options = optimset('Display', 'off', 'TolX', 1e-10, 'TolFun', 1e-10, 'MaxFunEvals', 100000, 'MaxIter', 100000);
p_fit = fminsearch(@fit_singleGaussian, Starting_p, options, edges, ncount);
curve_fit = p_fit(1) * exp( - (edges - p_fit(2)).^2 / p_fit(3).^2 );
Mu_fit = p_fit(2);
SD_fit = p_fit(3);
segTH(1) = Mu_fit + SD_fit*2;
segTH(2) = Mu_fit + SD_fit*3;
figure;
histogram(data, nbins, 'Normalization', 'count');
hold on;
plot(edges, curve_fit, '-r');
yl = ylim;
line([segTH; segTH], [yl' yl'], 'Color', [1 0 0], 'LineStyle', '--');
xlabel(xlabel_text);
ylabel('Count');
title([num2str(sum(data > segTH(1))) ' Trials > 2SD; ' num2str(sum(data > segTH(2))) ' Trials > 3SD; ']);

function result = movmean_xu(data, nwindow)
result = zeros(1, numel(data)-nwindow+1);
for i = 1:numel(data)-nwindow+1
    datainwindow = data(i:i+nwindow-1);
    result(i) = mean(datainwindow);
end