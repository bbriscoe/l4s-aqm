# blameshift_unresp_plot_stats.m
# Plot marking probability of two unresponsive bursty flows to illustrate blame-shifting
#  2 plots: sojourn-based and EST-based marking
#
# Called manually after blameshift_unresp.m run in non-qt_mode

# ToDo: Work out why sojourn plot is smaller within the window

set(0, "defaultlinelinewidth", 1.5);
## ToDo: The following are in ~/.octaverc but don't seem to affect the title or xlabel, ylabel
##set(0, "defaulttextfontsize", 18);
##set(0, "defaulttextfontname", "Times New Roman")
##set(0, "defaultaxesfontsize", 18);

linestyles = {":"; "--"; "-"};
linewidths = [2.5; 1.5; 1.5];
colours = [0.6, 0.6, 0.6;     # light grey
           1  , 0  , 1  ;     # magenta
           0.2, 0.2, 1  ;     # blue
           0.6, 0  , 0.6;     # purple
           0  , 0  , 0  ];    # black
slS = cast(lambdaSum, "int16");

# Plot p_s
figure(1)
for n = 1 : lambdaSum-1
  i_style = mod(n-slS,3)+1;
  ltag = strcat("p_{sa}: λ_a=^{", num2str(n), "}/_{", num2str(lambdas), "}");
  plot(beta(1,:), p_stats(n,:,1,1,1), 
       "DisplayName", ltag, ...
       "linestyle", linestyles{i_style}, ...
       "linewidth", linewidths(i_style), ...
       "color", colours(floor((n-1)/3)+1,:) );
  hold on;
endfor
hold off;
axis([0 double(betaSum)/betas 0 1.04]);
h_tit = title(["Sojourn marking of two unresponsive flows, a & b\n" ...
       "Capacity fractions, λ_a & λ_b;   Utilization, Σλ = " ...
       num2str(double(lambdaSum)/lambdas*100) "%;   " ...
       "Burst sizes β_a & β_b, where Σβ = " num2str(double(betaSum)/betas*100) ...
       "% of marking threshold"]);
h_xl = xlabel(["Normalized burst size of flow a,  β_a"]);
h_yl = ylabel("Sojourn-based marking probability of flow a, p_{sa}\n\n");
h_leg = legend();
set(h_tit, 'fontname', 'Times New Roman', 'fontsize', 20)
set(h_xl, 'fontname', 'Times New Roman', 'fontsize', 20)
set(h_yl, 'fontname', 'Times New Roman', 'fontsize', 20)
set(h_leg, 'location', 'northeastoutside', ...
    'fontname', 'Times New Roman', 'fontsize', 14)
set(gca, "outerposition", [0 0 1 1])
savefile = strrep(savefile, "_p_stats", "_p_s");
print([savefile ".pdf"]);

# Plot p_e
figure(2)
for n = 1 : lambdaSum-1
  i_style = mod(n-slS,3)+1;
  ltag = strcat("p_{ea}: λ_a=^{", num2str(n), "}/_{", num2str(lambdas), "}");
  plot(beta(1,:), p_stats(n,:,1,2,1), 
       "DisplayName", ltag, ...
       "linestyle", linestyles{i_style}, ...
       "linewidth", linewidths(i_style), ...
       "color", colours(floor((n-1)/3)+1,:) );
  hold on;
endfor
hold off;
axis([0 double(betaSum)/betas 0 1.04]);
h_tit = title(["Expected Service Time (EST) marking of two unresponsive flows, a & b\n" ...
       "Capacity fractions, λ_a & λ_b;   Utilization, Σλ = " ...
       num2str(double(lambdaSum)/lambdas*100) "%;   " ...
       "Burst sizes β_a & β_b, where Σβ = " num2str(double(betaSum)/betas*100) ...
       "% of marking threshold"]);
h_xl = xlabel(["Normalized burst size of flow a,  β_a"]);
h_yl = ylabel("EST-based marking probability of flow a,  p_{ea}\n\n");
h_leg = legend();
set(h_tit, 'fontname', 'Times New Roman', 'fontsize', 20);
set(h_xl, 'fontname', 'Times New Roman', 'fontsize', 20);
set(h_yl, 'fontname', 'Times New Roman', 'fontsize', 20);
set(h_leg, 'location', 'northeastoutside', ...
    'fontname', 'Times New Roman', 'fontsize', 14)
set(gca, "outerposition", [0 0 1 1])
savefile = strrep(savefile, "_p_s", "_p_e");
print([savefile ".pdf"]);
