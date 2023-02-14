# blameshift_unresp_plot_qt.m
# Plot time series of two unresponsive bursty flows to illustrate blame-shifting
#  Optionally display ECN marking of each flow
#
# Copyright (c) 2022-23 Bob Briscoe <research@bobbriscoe.net>
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

# Usage:
# 1. Edit the input parameters directly in this file and save.
# 2. Run this script from octave after running blameshift_unresp.m with
#     qt_mode on. See blameshift_unresp.m for details.

# Input parameters
est_marking = true(1);    # plot est-marked traffic darker if true
p_tab = true(1);         # write a table of values of p if true

figure(9);
h_a = area(qt_out(:,1), qt_out(:,2:6));
set(gca, "ygrid", "on", "gridalpha", 1, "xgrid", "off", "linewidth", 2, ...
         "ytick", [0 1]);
h_tit = title(["Two unresponsive flows, a \\& b;          " ...
       "_{Phase shift, φ = " num2str((phi(k)+smidgen)/phis*100) "%}"; ...
       "_{Capacity fractions, λ_a = " num2str(lambda(i,1)) ...
       "/" num2str(lambdas) ...
       ", λ_b = " num2str(lambda(i,2)) "/" num2str(lambdas) ...
       " (Σλ = " num2str(double(lambdaSum)/lambdas*100) "%);      " ...
       "Burst queue delays β_a = " num2str(beta(1,j)*100) ...
       "%, β_b = "  num2str(beta(2,j)*100) "% " ...
       "(Σβ = " num2str(double(betaSum)/betas*100) "%)}"]);
h_xl = xlabel("normalized time, t\n");
h_yl = ylabel("normalized queue delay, q\n\n");
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
  h_leg = legend("q_a unmarked", "q_a EST-marked ", ...
                 "q_b unmarked", "q_b EST-marked ", ...
                 "location", "northeastoutside");
else
  set(h_a(2), "facecolor", [0.7 0.7 0.7]);
  set(h_a(4), "facecolor", [1.0 0.8 0.8]);
  h_leg = legend("q_a  ", "", "q_b  ", "location", "northeastoutside");
endif
set(h_leg, 'fontname', 'Times New Roman', 'fontsize', 14)
# Hacky table of values of p
if (p_tab)
  p_text = ["\ta\t\tb";
            "p_s   ", num2str(p(1,:,1,1));
            "p_e   ", num2str(p(1,:,2,1))];
  text(1, 1.06*double(betaSum)/betas, p_text, "fontsize", 14);
endif
print([savefile ".pdf"], "-dpdfcairo");
