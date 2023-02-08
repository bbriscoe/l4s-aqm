# blameshift_unresp_plot_stats.m
# Plot marking probability of two unresponsive bursty flows to illustrate blame-
#  shifting.
# Plots up to four marking metrics using sojourn-based and/or EST-based marking
#  with and/or without max-min ranges, thus outputting up to 16 plots.

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
# 1. Edit the options line directly in this file and save.
# 2. Run this script from octave after running blameshift_unresp.m with
#     qt_mode off. See blameshift_unresp.m for details.

# Option processing
# plot options
# option{1,1} flow:      1 = p_a;     2 = Δp = p_a - p_b
# option{1,2} approach:  1 = sojourn; 2 = EST
# option{1,3} metric:    1 = p;       2 = λp
# option{1,4} statistic: 1 = mean;    2 = mean & max-min
# Default: option = {1:2, 1:2, 1:2, 1};
option = {1, 1, 1, 1};

set(0, "defaultlinelinewidth", 1.5);
## ToDo: The following are in ~/.octaverc but don't seem to affect the title or
##  xlabel, ylabel
##set(0, "defaulttextfontsize", 18);
##set(0, "defaulttextfontname", "Times New Roman")
##set(0, "defaultaxesfontsize", 18);

linestyles = {":"; "--"; "-"};
linewidths_s = [2.5 ; 1.5 ; 1.5 ]; # For on screen
linewidths_p = [2.0; 1.0; 1.0]; # For print -dpdfcairo
colours        = [0.6 , 0.6 , 0.6 ;     # light grey
                  1   , 0.2 , 1   ;     # magenta
                  0.2 , 0.2 , 1   ;     # blue
                  0.6 , 0   , 0.6 ;     # purple
                  0   , 0   , 0   ];    # black
# Create map of brighter but related colours for shaded fill
scale = 8;
colours(:,:,2) = (colours(:,:,1)+scale)/(1+scale);
slS = cast(lambdaSum, "int16");

for statistic = option{1,4}
  for metric = option{1,3}
    for flow = option{1,1}
      for approach = option{1,2}
        ## 1 var_tag = "p_{sa}";
        ## 2 var_tag = "p_{ea}";
        ## 3 var_tag = "Δp_{s}";
        ## 4 var_tag = "Δp_{e}";
        ## 5 var_tag = "λ_ap_{sa}";
        ## 6 var_tag = "λp_{ea}";
        ## 7 var_tag = "Δ(λp_{s})";
        ## 8 var_tag = "Δ(λp_{e})";
        fig = 4*(metric-1)+2*(flow-1)+(approach-1) + 1
        if (approach == 1)
          appr_tag = "s";
          tit_tag = "Sojourn marking with";
          if (flow == 1)
            ylab_tag = "Sojourn marking";
          else
            ylab_tag = "Diff. betw. sojourn marking";
          endif
        else
          appr_tag = "e";
          tit_tag = "Expected Service Time (EST) marking with";
          if (flow == 1)
            ylab_tag = "EST marking";
          else
            ylab_tag = "Diff. betw. EST marking";
          endif
        endif
        if (metric == 1)
          lam_tag = "";
        else
          lam_tag = "λ";
        endif
        if (flow == 1)
          ylabm_tag = " probability of";
          if (metric == 2)
            ylabm_tag = " rate of";
          endif
          ylabf_tag = " flow a, ";
          delt_tag = del_tag = "";
          flow_tag = "a}";
          lylim = 0;
        else
          ylabm_tag = " probabilities";
          ylabf_tag = [" " lam_tag "p_a - " lam_tag "p_b, "];
          delt_tag = del_tag = "Δ";
          flow_tag = "}";
          if (metric == 2)
            ylabm_tag = " rates";
            delt_tag = [del_tag "("];
            flow_tag = [flow_tag ")"];
          endif
          lylim = -1.04;
        endif
        file_tag = ["_" del_tag lam_tag "p_" appr_tag];
        var_tag = [delt_tag lam_tag "p_{" appr_tag flow_tag];
        stag = "";

        figure(fig, "paperunits", "centimeters", "paperposition", ...
                    [0, 0, 27.25, 13.75])
        # Plot loop
        for n = 1 : lambdaSum-1
          i_style = mod(n-slS,3)+1;
          ltag = strcat(var_tag, ": λ_a=^{", num2str(n), "}/_{", ...
          num2str(lambdas), "}");
          if (statistic == 2)
            # Fill an area from max to min behind the subsequent plot of mean
            stag = "_{min-max} ";
            fill([beta(1,:), fliplr(beta(1,:))], ...
                 [p_stats(n,:,flow,approach,metric,2), ...
                  fliplr(p_stats(n,:,flow,approach,metric,3))], ...
                 colours(floor((n-1)/3)+1,:,2), ...
                 "DisplayName", [stag ltag], ...
                 "linestyle", "none",...
                 "facealpha", 0.5);
            stag = "    _{mean} ";
            hold on;
          endif
          h_plot(n) = plot(beta(1,:), p_stats(n,:,flow,approach,metric,1), ...
               "DisplayName", [stag ltag], ...
               "linestyle", linestyles{i_style}, ...
               "linewidth", linewidths_p(i_style), ...
               "color", colours(floor((n-1)/3)+1,:,1) );
               # ToDo: bug: when (LambdaSum-1) is not a multiple of 3
               #  linestyles don't sync with line colours.
          hold on;
        endfor
        hold off;
        axis([0 double(betaSum)/betas lylim 1.04]);
        h_tit = title([tit_tag " two unresponsive flows, a \\& b\n" ...
               "_{Capacity fractions, λ_a \\& λ_b : utilization, Σλ = " ...
               num2str(double(lambdaSum)/lambdas*100) "%;   " ...
               "Burst sizes β_a \\& β_b : Σβ = " ...
               num2str(double(betaSum)/betas*100) ...
               "% of marking threshold}"]);
        h_xl = xlabel(["Normalized burst size of flow a,  β_a"]);
        h_yl = ylabel([ylab_tag ylabm_tag ylabf_tag var_tag "\n\n"]);
        h_leg = legend();
        set(h_tit, 'fontname', 'Times New Roman', 'fontsize', 20)
        set(h_xl, 'fontname', 'Times New Roman', 'fontsize', 20)
        set(h_yl, 'fontname', 'Times New Roman', 'fontsize', 20)
        set(h_leg, 'location', 'northeastoutside', ...
            'fontname', 'Times New Roman', 'fontsize', 14)
        set(gca, "outerposition", [0 0 1 1])
        printfile = strrep(savefile, "_p_stats", file_tag);
        print([printfile ".pdf"], '-dpdfcairo');
        for n = 1 : lambdaSum-1
          set(h_plot(n), "linewidth", linewidths_s(mod(n-slS,3)+1));
        endfor
      endfor
    endfor
  endfor
endfor
