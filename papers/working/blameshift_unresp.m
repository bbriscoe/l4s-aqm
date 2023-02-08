#! /bin/octave -qf
# blameshift_unresp.m
# Scan the blame-shifting problem space

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
#  1. Edit primary and secondary parameters directly in this file and save.
#  2. Run either from within octave as
#       >> blameshift_unresp
#     or from the command line as
#       $ octave blameshift_unresp.m
#  3. A binary (.bin) file of the results is saved at the end, in subdirectory
#      octave_data if it exists, otherwise in the current directory.
#  4. To plot output at any time, this saved binary can be loaded into octave:
#      >> load <file.bin>
#      But this is unnecessary if plotting is done immediately after running
#      this script, while the output is still in memory.
#     Then a second script must be run from octave, depending respectively on
#      whether qt_mode was set to true(1) or false(1) in this script:
#       >> blameshift_unresp_plot_qt
#       >> blameshift_unresp_plot_stats
#      See the source of each of these scripts for their usage.

# Testing:
#  octave 5.2.0 under Ubuntu 20.04 LTS
#  octave 7.3.0 under Win 10
# Not tested with MatLab

## ToDo:
## * Tidy code into functions before publication

clear
debug_mode = true(1);  # Debug mode

# Primary parameters
lambdaSum = 1;		# Utilization
betaSum = 5/4;	# Max combined normalized burst delay (wrt marking threshold)
if (lambdaSum > 1)
  error("utilization parameter 'lambdaSum' cannot exceed 100%");
endif

# Secondary parameters
# lambda, beta & phi are held as integers 
#  with implicit denominators lambdas, betas & phis, resp.
#  (except after setup, when beta is cast to double and downscaled), 
lambdas = 16;		# no. of divisions of capacity share, lambda, if lambdaSum=1
betas = 32;		# no. of divisions of normalized burst delay, beta, if betaSum=1
phis = 4;		# no. of divisions of phase shift, phi, in 360deg
smidgen = 0.123456789;  # To avoid unrealistic degree of exact phase lock
#
# set qt_mode to true(1) to produce one time series of the queue
# set qt_mode to false(1) to scan parameter space and produce marking statistics
qt_mode = true(1);
if (qt_mode)
  i_lambda = 2; # index of lambda to plot if in qt_mode
  i_beta = 6;  # index of beta   to plot if in qt_mode
  i_phi = 1;    # index of phi    to plot if in qt_mode
endif

savesuf = ["_", strftime("%Y%m%d-%H%M%S", localtime (time ()))];
savedir = "";
if (exist("octave_data", "dir"))
  savedir = "octave_data/";
endif
if (qt_mode)
  savefile = [savedir mfilename(), "_qt_out_", ...
              "Σλ", num2str(lambdaSum), "_Σβ", num2str(betaSum), ...
              "_", num2str(i_lambda), "λ", num2str(lambdas), ...
              "_", num2str(i_beta), "β", num2str(betas),...
              "_", num2str(i_phi), "φ", num2str(phis), ...
              savesuf];
else
  savefile = [savedir mfilename(), "_p_stats_", ...
              "Σλ", num2str(lambdaSum), "_Σβ", num2str(betaSum), ...
              "_β", num2str(betas), ...
              savesuf];
endif
clear savedir savesuf;

# Upscale primary parameters to integers
lambdaSum *= lambdas;
betaSum *= betas;
lambdaSum = cast(lambdaSum, "uint16");
betaSum = cast(betaSum, "uint16");
if (qt_mode)
  if (i_lambda < 1 || i_lambda > lambdaSum-1)
    error("capacity share index parameter 'i_lambda' outside valid range");
  endif
  if (i_beta < 1 || i_beta > betaSum-1)
    error("burst size index parameter 'i_beta' outside valid range");
  endif
  if (i_phi < 1 || i_phi > phis)
    error("burst size index parameter 'i_phi' outside valid range");
  endif
endif

# Fill vectors for scanning the problem space
# Two flows, a&b are represented by slices 1&2
# The second flow's capacity share and burst size use up the remainder of
#  lambdaSum and betaSum
lambda = 1 : lambdaSum - 1;
lambda = [lambda; lambdaSum - lambda]';
beta = 1 : betaSum - 1;
beta = [beta; betaSum - beta];
phi = 0 : phis - 1;
#
#                     [ 1, ... beta(1,:) ... betasSum - 1 ]
#  _              _   _         beta(2,:)                _ ]
# |                | |                                    |
# | 1              | |                                    ||
# | .              | |                                    ||
# | .              | |                                    ||
# | .              | |                                    ||
# | lambda(:,1)    | |         ti(:,:,1)                  ||
# | .              | |                                    ||
# | .              | |                                    ||
# | .              | |                                    ||
# |_LambdaSum - 1 _| |_                                  _||
#  |_lambda(:,2)   _| |_            ti(:,:,2)             _|

# Interval between bursts, ti = beta/lambda
# Compute matrix of intervals (in unscaled time units)
ti(:,:,1) = (double(beta(1,:)) ./ double(lambda(:,1)) * lambdas/betas);
ti(:,:,2) = (double(beta(2,:)) ./ double(lambda(:,2)) * lambdas/betas);
if (min(min(min(ti))) == 0)
  error("At least one interval between bursts 'ti' is zero");
endif

# No need to run time beyond when the pattern of the two cyclic bursts repeats,
#  which is the lowest common multiple of the intervals of flows a & b
#  t_max = lcm_frac(ti(:,:,1), ti(:,:,2)), or for any one slice of the matrix
#        = lcm_frac(beta_1/lambda_1, beta_2/lambda_2)
#  where lcm_frac is a notional lcm function that would work on fractions
#  and beta & lambda are the fractional values of the input parameters.
# In practice the parameters of an lcm function must be integers. 
# So, multiply both input parameters by betas*lambdas * lambda_1*lambda_2
#  and divide the result by the same.
#  t_max = lcm(betas*lambdas*beta_1*lambda_2, betas*lambdas*beta_2*lambda_1) /
#                             betas*lambdas * lambda_1*lambda_2

ti_tmp(:,:,1) = beta(1,:) .* lambda(:,2);
ti_tmp(:,:,2) = beta(2,:) .* lambda(:,1);
t_max = lcm(uint32(ti_tmp(:,:,1)), uint32(ti_tmp(:,:,2)));
if (lambdas >= betas)
  t_max *= lambdas/betas;  # multiply before saturation check
endif
clear ti_tmp;
if (max(max(t_max)) == intmax("uint32"))
  error("lcm() saturated at least one uint32 element of 't_max'");
endif
# Cast t_max to FP and downscale
t_max = double(t_max) ./ double(lambda(:,1) .* lambda(:,2));
if (lambdas < betas)
  t_max *= lambdas/betas;  # multiply after saturation check
endif
# Cast beta to FP and downscale
beta = double(beta) ./ betas;

if (!qt_mode)
  i_lambda = 1 : lambdaSum - 1;
  i_beta = 1 : betaSum - 1;
  i_phi = 1 : phis;
  p_stats = zeros(lambdaSum-1,betaSum-1,2,2,2,3);
endif

debug_mode && (i_bug = 0);
for (i = i_lambda)
  for (j = i_beta)
    if (!qt_mode)
      printf(".");
    endif
    # i_freq indexes the flow with more frequent bursts (or smaller burst size
    #  in case of a tie)
    if (ti(i,j,1) != ti(i,j,2))
      [~, i_freq] = min(ti(i,j,:));
    else
      [~, i_freq] = min(beta(:,j));
    endif
    i_rare = !(i_freq-1) + 1;
    
    # As long as utilization < 100% and bursts are evenly spaced, there cannot 
    #  be >1 sequence of packets from each flow in the queue at once.
    # So a vector of 2 q sizes and a 2-state integer 'i_head' indexing the flow
    #  at the head are sufficient to fully describe the queue at any one time

    # Marking prob, p, for each flow & each approach
    if (!qt_mode)
      # For each phi, the 2x2 matrix of p is accumulated as time is scanned,
      #  then normalized once t_max is reached
      #          [ a , b ]
      #  _    _   _     _
      # |  1   | |       |_
      # |  .   | |       | |
      # |  .   | |       | |
      # | phi  | |       | |
      # |  .   | |       | |
      # |  .   | | p_s   | |
      # |_phis_| |_     _| |
      #            |_p_e  _|
      #
      p = zeros(phis,2,2);
    else
      p = zeros(1,2,2);
    endif
    
    # k : index of phi
    for (k = i_phi)
      # Variable definitions
      #  t :          current time
      #  i_head :     which flow is at the head of the queue
      #  t_burst[2] : start time of next burst from each flow
      #  i_bNxt :     which flow's burst is next
      
      # Time is scanned in two passes:
      # 1) to find where the q will be shortest (using only matrix ops)
      # 2) to derive the queue process & marking event-by-event in a while loop
      
      # #1 time scan
      # Set the origin (t=0 & q=0) where the combined queue will be shortest, 
      #  assuming utilization <= 100%, and assuming a model queue that starts 
      #  full enough to never empty, draining it at the avg arrival rate so that
      #  the pattern of burst arrivals maintains a standing q.
      # Only rare bursts and the freq bursts just before them ('pre-rare') 
      #  need to be scanned, because the q at the start of any other freq burst 
      #  will be greater.
      # 0) Pick any phase difference to be the first between a pre-rare freq 
      #    burst and the subsequent rare burst
      # 0) 'Seed' the initial phase shift using fraction phi of freq interval
      t_delta0 = ti(i,j,i_freq) * (phi(k)+smidgen)/phis;
      # Vector of rare burst times
      t_rare = t_delta0 : ti(i,j,i_rare) : t_max(i,j);
      # 1) Create vector t_rare of each phase shift from each pre-rare to rare
      #    burst
      t_rare = rem(t_rare, ti(i,j,i_freq));
      # 2) Create vector q_rare of the no. (say n) of freq bursts from each pre-
      #    rare burst to the next pre-rare burst. 
      q_rare = ceil((ti(i,j,i_rare) - t_rare) / ti(i,j,i_freq));
      # 3) Then the change in the queue from each pre-rare burst to the next can 
      #    be calculated as: 
      #    1*rare burst + n*(freq burst - drain per freq interval)
      #    This will give a sequence of -ve and +ve changes that sum to zero
      drain_freq = ti(i,j,i_freq) * double(lambdaSum) / lambdas;
      q_rare = beta(i_rare,j) + q_rare .* (beta(i_freq,j) - drain_freq);
      # 4) Convert these to a cumulative sum. 
      #    This will give the queue at the start of each pre-rare burst relative
      #    to the q at the start. 
      #    The min of these is candidate #1 for the location of the origin.
      q_rare = cumsum(q_rare);
      # 5) Then the queue at the start of each rare burst can be derived by 
      #    adding the (pre-rare) freq burst and subtracting the drain rate over
      #    the phase shift still stored in in t_rare.
      #    The min of these is candidate #2 for the location of the origin.
      q_rare(2,:)= q_rare + beta(i_freq,j) - t_rare * double(lambdaSum)/lambdas;
      # Find the mins of both vectors, and the indices of their locations
      [qr_min, i_qr_min] = min(q_rare,[],2);
      #
      # The outputs from the first pass are:
      #  t_burst{2} : all the times when each flow will burst, starting when q=0
      #               (a cell array containing 2 vectors of different lengths)
      #  i_tb[2] :    index of the next burst in each vector
      #  i_bNxt :     index of the next flow to burst
      #  i_head :     index of the flow currently at the head of the q
      if (qr_min(1) < qr_min(2))
        # t=0, q=0 at start of freq burst
        t_burst{i_freq} = 0 : ti(i,j,i_freq) : t_max(i,j);
        t_burst{i_rare} = t_rare(i_qr_min(1)) : ti(i,j,i_rare) : t_max(i,j) ...
                        + t_rare(i_qr_min(1)); # just for limit test
        if (t_burst{i_rare}(1) > 0)
          i_head = i_freq;
        else
          # Tie-break if initial bursts coincide
          [~, i_head] = min(beta(:,j));
        endif
      else
        # t=0, q=0 at start of rare burst
        t_burst{i_rare} = 0 : ti(i,j,i_rare) : t_max(i,j);
        t_burst{i_freq} = ti(i,j,i_freq) - t_rare(i_qr_min(2)) ... # >= 0
                             : ti(i,j,i_freq) : t_max(i,j) ...
             + ti(i,j,i_freq) - t_rare(i_qr_min(2)); # just for limit test
        i_head = i_rare;
      endif
      i_bNxt = i_head;
      i_tb = ones(2,1);
      !debug_mode && clear('t_delta0', 't_rare', 'q_rare', 'drain_freq', ...
                            'qr_min', 'i_qr_min');
      
      # #2 time scan
      # Initialize and define variables
      t = 0;              # current time
      q = zeros(2,1);     # queue delay contributed by each flow
      ott = 0;            # whether combined queue is over the threshold (q>1)
      qt_mode && (i_event = 0);      # Event index
      #
      # The qt_out matrix is extended a row at a time, as time is scanned:
      #  * qt_out(:,1) : time of next event; 
      #  * qt_out(:,2:3) : resulting per flow queue contribution at that time; 
      #  * qt_out(:,4) : which flow is at the head of the q; 0:1 means q(1:2);
      #  * qt_out(:,5) : whether queue is above threshold (q>=1) after t.
      while (t < t_max(i,j))  # redundant 'cos break always preempts
        # The contributions from each flow to the queue are piecewise linear
        #  between 'events', where an 'event' is a discontinuity in one of the
        #  contributions or when the queue crosses the marking threshold
        # 
        t_emptyNxt = t + q(i_head); # Time the head q will next empty
        q_xs = sum(q) - 1;          # Excess q above threshold (can be -ve)

        # Define t_burstNxt for readability & brevity, but take care! it's not a
        #  macro; it uses the values of i_bNxt and i_tb now, not when it is used
        t_burstNxt = t_burst{i_bNxt}(i_tb(i_bNxt)); # time of next burst
        # Check whether combined q falls below threshold 
        #  before head q empties or burst arrives
        if ( (ott) && ...
             ( ((t_emptyNxt <= t_burstNxt) && (q_xs <= q(i_head))) || ...
               ((t_emptyNxt > t_burstNxt) && (q_xs <= t_burstNxt - t)) ...
             ) ...
           )
          # Combined queue has fallen below threshold before any other event
          t += q_xs;
          # Add to p_e
          p(k, i_head, 2) += q_xs;
          q(i_head) -= q_xs;
          qt_mode && (qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott]);
          ott = 0;
          # The earlier condition (q_xs <= t_burstNxt - t) could have been
          #  changed to < to suppress the following qt_out in the case when it
          #  seems redundant when q drains exactly to the threshold just
          #  before a burst. But, for robustness, if (q==1) ott=0, even if q
          #  is about to increase again.
          qt_mode && (qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott]);
          # Check for q(i_head) emptying just as tail crosses threshold
          if (t >= t_emptyNxt)
            # Point i_head to other flow's queue if it's non-negative
            i_other = !(i_head-1) + 1;
            if (q(i_other) > 0)
              qt_mode && (qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott]);
              i_head = i_other;
              # No need to update i_other - not read elsewhere
            endif
          endif
        else
          if ((q(i_head) > 0) && (t_emptyNxt <= t_burstNxt))
            # q(i_head) has emptied (defer any simultaneous burst to next event)
            t = t_emptyNxt;
            if (ott)
              # Altho head flow has emptied, tail is over threshold
              #  so add time since previous event to p_e
              p(k, i_head, 2) += q(i_head);
            endif
            q(i_head) = 0;
            # Point i_head to other flow's queue if it's non-negative
            i_other = !(i_head-1) + 1;
            if (q(i_other) > 0)
              qt_mode && (qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott]);
              i_head = i_other;
              # No need to update i_other - not read elsewhere
            else
              # Both q's empty
              i_head = i_bNxt;
            endif
            qt_mode && (qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott]);
          else
            # Burst has arrived
            if (q(i_head) > 0)
              # Drain head queue since last event
              q(i_head) -= t_burstNxt - t;  # resulting q(i_head) will be >0
              if (ott)
                # Add to p_e
                p(k, i_head, 2) += t_burstNxt - t;
              endif
              q_xs = sum(q) - 1;
            endif
            t = t_burstNxt;
            qt_mode && (qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott]);
            # The precision of t is 1 eps 'cos it is assigned from a range.
            # However, an octave bug loses another eps when within a script
            # The alternative of testing for the end of the vector led to very
            #  complex code
            if (t >= t_max(i,j) - 2*eps(t))
              break
            endif
            # Add burst to tail
            #  but first check whether combined queue rises above threshold
            delta_q = beta(i_bNxt,j);
            if ( (ott == 0) && (q_xs + beta(i_bNxt,j) > 0) )
              # Care! if (q_xs + beta == tiny) (e.g tiny = eps()), after q is
              #  incremented below, the tiny value disappears and q_xs =
              #  (sum(q) - 1) == 0.
              #  So subsequently, we test the result (ott), not q_xs again.
              q(i_bNxt) -= q_xs;
              delta_q += q_xs;
              qt_mode && (qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott]);
              ott = 1;
            endif
            q(i_bNxt) += delta_q;
            if (ott)
              # Add to p_s
              p(k, i_bNxt, 1) += delta_q;
            endif
            qt_mode && (qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott]);
            # Prepare for next burst
            # This flow's burst is complete, so set t_burstNxt for the same flow
            # Earlier break prevents overflow here
            i_tb(i_bNxt)++;
            # Calc arrival time and flow id of next burst, handling tie if nec.
            [~, i_bNxt] = min([t_burst{1}(i_tb(1)), t_burst{2}(i_tb(2))]);
            if (t_burst{1}(i_tb(1)) == t_burst{2}(i_tb(2)))
              [~, i_bNxt] = min(beta(:,j));
            endif
          endif    
        endif
      endwhile
      if (debug_mode && (sum(q) > 0))
        qt_bug(++i_bug,:) = [t, q(1), q(2), i_head-1, t_max(i,j), i, j, k];
      endif
    endfor
    
    # Add another dimension: p[phis,2,2,2], such that
    # 1) phase shift,    1 : phis
    # 2) flow,           a, b
    # 3) approach,       soj, est
    # 4) marking metric, p (marking prob), λp (normalized marking rate)
    p(:,:,:,2) = p ./t_max(i,j);
    p(:,:,:,1) = p(:,:,:,2) ./ double(lambda(i,:))*lambdas;
    
    if (!qt_mode)
      # (a - b) replaces b's column
      p(:,2,:,:) = p(:,1,:,:)-p(:,2,:,:);
      # p_stats(i,j,flow-id,approach,metric,statistic), where
      # i       : index of lambda
      # j       : index of beta
      # flow-id : 1,2 = flow a,(a-b)
      # approach: 1,2 = soj,est
      # metric  : 1,2 = p, λp
      # statistic : 1,2,3 = mean,max,min
      p_stats(i,j,:,:,:,1) = mean(p(1:phis,:,:,:));
      p_stats(i,j,:,:,:,2) =  max(p(1:phis,:,:,:));
      p_stats(i,j,:,:,:,3) =  min(p(1:phis,:,:,:));
    endif

  endfor
  !qt_mode && (printf("%d%%\n", double(i)/double(lambdaSum)*100));
endfor

if (qt_mode)
  # post-process qt_out matrix
  # To simplify plotting, split up flows into EST-marked and unmarked columns
  #  and split flow 'a' into head and tail columns, resulting in:
  # 1: time of next event
  # 2: q_a if at head; unmarked
  # 3: q_a if at head; marked
  # 4: q_b unmarked
  # 5: q_b marked
  # 6: q_a if at tail
  # 7: i_head-1
  # 8: ott
 qt_out = [qt_out(:,1), ... 
           qt_out(:,2) .*  !qt_out(:,4) .* !qt_out(:,5), ...
           qt_out(:,2) .*  !qt_out(:,4) .*  qt_out(:,5), ...
           qt_out(:,3) .* !(qt_out(:,4) .*  qt_out(:,5)), ...
           qt_out(:,3) .*  (qt_out(:,4) .*  qt_out(:,5)), ...
           qt_out(:,2) .*   qt_out(:,4), ...
           qt_out(:,4:5)];
  save("-binary", [savefile ".bin"], "savefile", "qt_out", "p");
else
  save("-binary", [savefile ".bin"], "savefile", ...
       "lambdaSum", "lambdas", "betaSum", "betas", "beta", "p_stats");
endif

## ============================================================================
## # Manifest of variables
## # =====================
## 
## # Input Parameters
## lambda           [lambdaSum-1, 2         ]
## beta             [2,           betaSum-1 ]
## phi              [phis        ]
## qt_mode          logical
## debug_mode       logical
## 
## # Limits
## lambdaSum
## betaSum
## lambdas
## betas
## phis
##
## # Constant
## smidgen
## 
## # Indices (mostly prefixed by 'i_')
## i_lambda, i : index of lambda (range and single value)
## i_beta,   j : index of beta   (range and single value)
## i_phi,    k : index of phi    (range and single value)
## i_event
## i_freq
## i_rare
## i_head
## i_other
## i_bNxt
## i_qr_min         [2]
## i_tb             [2]
## i_bug
##
## # Time intervals
## t
## t_delta0
## t_rare           [1, nnn]
## ti               [lambdaSum-1, betaSum-1,      2 ]
## t_max            [lambdaSum-1, betaSum-1 ]
## t_burst          {2  }
## t_burstNxt
## t_emptyNxt
##
## # Queue delays
## drain_freq
## q_rare           [2, nnn]
## qr_min           [2]
## q_xs
## delta_q
##
## # Outputs
## q                [2  ]
## ott
## qt_out           [mmm,         8]   ([mmm, 5] until post-processed)
## qt_bug           [rrr,         8]
## p                [phis+1,      2,              2 ]
## p_stats          [lambdaSum-1, betaSum-1, 2, 2, 3]
##
## # Strings
## savedir
## savesuf
## savefile
## ============================================================================


