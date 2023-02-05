# blameshift_unresp_plot_stats.m
# Plot marking probability of two unresponsive bursty flows to illustrate blame-shifting
#  2 plots: sojourn-based and EST-based marking
#
# Called manually after blameshift_unresp.m run in non-qt_mode

set(0, "defaultlinelinewidth", 1.5);
## ToDo: The following are in ~/.octaverc but don't seem to affect the title or xlabel, ylabel
##set(0, "defaulttextfontsize", 18);
##set(0, "defaulttextfontname", "Times New Roman")
##set(0, "defaultaxesfontsize", 18);

linestyles = {":"; "--"; "-"};
linewidths = [2.5; 1.5; 1.5];
colours        = [0.6 , 0.6 , 0.6 ;     # light grey
                  1   , 0   , 1   ;     # magenta
                  0.2 , 0.2 , 1   ;     # blue
                  0.6 , 0   , 0.6 ;     # purple
                  0   , 0   , 0   ];    # black
colours(:,:,2) = [0.95, 0.95, 0.95;     # light grey
                  1   , 0.95, 1   ;     # magenta
                  0.9 , 0.9 , 1   ;     # blue
                  0.9 , 0.8 , 0.9 ;     # purple
                  0.9 , 0.9 , 0.9 ];    # black
slS = cast(lambdaSum, "int16");

# Option processing
# plot options
# option{1,1} flow:      1 = p_a;     2 = Δp = p_a - p_b
# option{1,2} approach:  1 = sojourn; 2 = EST
# option{1,3} metric:    1 = p;       2 = λp
# option{1,4} statistic: 1 = mean;    2 = mean & std. dev
# Default: option = {1:2, 1:2, 1:2, 1};
option = {1, 2, 1, 2};

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

        # Plot loop
        figure(fig)
        for n = 1 : lambdaSum-1
          i_style = mod(n-slS,3)+1;
          ltag = strcat(var_tag, ": λ_a=^{", num2str(n), "}/_{", num2str(lambdas), "}");
          if (statistic == 2)
            # Fill an area from max to min behind the mean plotted next
            fill([beta(1,:), fliplr(beta(1,:))], ...
                 [p_stats(n,:,flow,approach,metric,2), ...
                  fliplr(p_stats(n,:,flow,approach,metric,3))], ...
                 colours(floor((n-1)/3)+1,:,2), "linestyle", "none");
            hold on;
          endif
          plot(beta(1,:), p_stats(n,:,flow,approach,metric,1), ...
               "DisplayName", ltag, ...
               "linestyle", linestyles{i_style}, ...
               "linewidth", linewidths(i_style), ...
               "color", colours(floor((n-1)/3)+1,:,1) );
          hold on;
        endfor
        hold off;
        axis([0 double(betaSum)/betas lylim 1.04]);
        h_tit = title([tit_tag " two unresponsive flows, a & b\n" ...
               "Capacity fractions, λ_a & λ_b;   Utilization, Σλ = " ...
               num2str(double(lambdaSum)/lambdas*100) "%;   " ...
               "Burst sizes β_a & β_b, where Σβ = " num2str(double(betaSum)/betas*100) ...
               "% of marking threshold"]);
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
        print([printfile ".pdf"]);
      endfor
    endfor
  endfor
endfor
