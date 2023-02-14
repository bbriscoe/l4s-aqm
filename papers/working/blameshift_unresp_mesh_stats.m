# blameshift_unresp_mesh_stats.m
# Plot marking probability of two unresponsive bursty flows to illustrate blame-
#  shifting.
# Plots up to four marking metrics using sojourn-based and/or EST-based marking,
#  thus outputting up to 8 plots.

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
# mesh options
# option{1,1} flow:      1 = p_a;     2 = Δp = p_a - p_b
# option{1,2} approach:  1 = sojourn; 2 = EST
# option{1,3} metric:    1 = p;       2 = λp
# Default: option = {1, 1:2, 1};
option = {1, 1:2, 1};

for metric = option{1,3}
  for flow = option{1,1}
    for approach = option{1,2}
      fig = 4*(metric-1)+2*(flow-1)+(approach-1) + 10
      if (approach == 1)
        appr_tag = "s";
        tit_tag = "Sojourn marking with";
        if (flow == 1)
          zlab_tag = "Soj. marking";
        else
          zlab_tag = "Diff. betw. soj. marking";
        endif
      else
        appr_tag = "e";
        tit_tag = "Expected Service Time (EST) marking with";
        if (flow == 1)
          zlab_tag = "EST marking";
        else
          zlab_tag = "Diff. betw. EST marking";
        endif
      endif
      if (metric == 1)
        lam_tag = "";
      else
        lam_tag = "λ";
      endif
      if (flow == 1)
          zlabm_tag = " probability";
        if (metric == 2)
          zlabm_tag = " rate";
        endif
        zlabf_tag = " of flow a, ";
        delt_tag = del_tag = "";
        flow_tag = "a}";
        lzlim = 0;
      else
        zlabm_tag = " probabilities";
        zlabf_tag = [" of " lam_tag "p_a - " lam_tag "p_b, "];
        delt_tag = del_tag = "Δ";
        flow_tag = "}";
        if (metric == 2)
          zlabm_tag = " rates";
          delt_tag = [del_tag "("];
          flow_tag = [flow_tag ")"];
        endif
        lzlim = 1.2;
      endif
      file_tag = ["_3d_" del_tag lam_tag "p_" appr_tag];
      var_tag = [delt_tag lam_tag "p_{" appr_tag flow_tag];

      figure(fig)
      h_mesh = mesh(beta(1,:)', double(lambda(1:(lambdaSum-1),1))/lambdas, ...
                    p_stats(:,:,flow,approach,metric,1));
      axis([-Inf Inf -Inf Inf lzlim 1.2]);
      h_xl = xlabel ("\nNormalized burst size of flow a, β_a");
      h_yl = ylabel ("\nCapacity share of flow a, λ_a");
      h_zl = zlabel ([zlab_tag zlabm_tag zlabf_tag var_tag "\n"]);
      h_tit = title ([tit_tag " two unresponsive flows, a \\& b"; ...
                     "_{Capacity fractions, λ_a \\& λ_b : utilization, Σλ = " ...
                     num2str(double(lambdaSum)/lambdas*100) "%;   " ...
                     "Burst sizes β_a \\& β_b : Σβ = " ...
                     num2str(double(betaSum)/betas*100) ...
                     "% of marking threshold}"]);
      set(h_tit, 'fontname', 'Times New Roman', 'fontsize', 20)
      set(h_xl, 'fontname', 'Times New Roman', 'fontsize', 20)
      set(h_yl, 'fontname', 'Times New Roman', 'fontsize', 20)
      set(h_zl, 'fontname', 'Times New Roman', 'fontsize', 20)
      printfile = strrep(savefile, "_p_stats", file_tag);
      print([printfile ".pdf"], '-dpdfcairo');
    endfor
  endfor
endfor
