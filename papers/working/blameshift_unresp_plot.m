# blameshift_unresp_plot.m
# Plot time series of two unresponsive bursty flows to illustrate blame-shifting
#  Optionally display ECN marking of each flow
#
# Called manually after blameshift_unresp.m

# Input parameter
est_marking = true(1);    # plot est-marked traffic darker if true

h_a = area(qt_out(:,1), qt_out(:,2:6));
set(gca, "ygrid", "on", "gridalpha", 1, "xgrid", "off", "linewidth", 2, "ytick", [0 1]);
h_tit = title(["Two unresponsive flows, a & b"; ...
       "Capacity fractions, λ_a = " num2str(lambda(i,1)) "/" num2str(lambdas) ...
       ", λ_b = " num2str(lambda(i,2)) "/" num2str(lambdas) ...
       " (Σλ = " num2str(double(lambdaSum)/lambdas*100) "%);"; ...
       "Burst queue delays β_a = " num2str(beta(1,j)*100) ...
       "%, β_b = "  num2str(beta(2,j)*100) "% " ...
       "(Σβ = " num2str(double(betaSum)/betas*100) "%)"]);
h_xl = xlabel("time, t  (normalized to marking threshold = 1)\n");
h_yl = ylabel("queue delay, q  (normalized to marking threshold = 1)\n\n");
set(h_a, "edgecolor", "none");
set(h_tit, 'fontname', 'Times New Roman', 'fontsize', 20)
set(h_xl, 'fontname', 'Times New Roman', 'fontsize', 20)
set(h_yl, 'fontname', 'Times New Roman', 'fontsize', 20)
set(h_a(1), "facecolor", [0.7 0.7 0.7]);
set(h_a(3), "facecolor", [1.0 0.8 0.8]);
set(h_a(5), "facecolor", [0.7 0.7 0.7]);
if (est_marking)
  set(h_a(2), "facecolor", [0.4 0.4 0.4]);
  set(h_a(4), "facecolor", [0.7 0.5 0.5]);
  h_leg = legend("q_a unmarked", "q_a EST-marked ", "q_b unmarked", "q_b EST-marked ", ...
                 "location", "northeastoutside");
else
  set(h_a(2), "facecolor", [0.7 0.7 0.7]);
  set(h_a(4), "facecolor", [1.0 0.8 0.8]);
  h_leg = legend("q_a  ", "", "q_b  ", "location", "northeastoutside");
endif
set(h_leg, 'fontname', 'Times New Roman', 'fontsize', 14)
savefile = [savepre, "_qt_out_", savemid, savesuf];
print([savefile ".pdf"]);
