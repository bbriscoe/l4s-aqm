#! /bin/octave -qf
# Scan the blame-shifting problem space
# Copyright (c) 2022-23 Bob Briscoe 
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

## ToDo:
## * Check whether code completes properly at each t_max limit
##   - scan #1 & #2
##   - q = 0 prior to last event (investigate bug-causes of non-zero cases)
##   - qt_mode plot cases where q!=0 at t_max
## * Output table of values of p on qt_mode plots
## * Understand why 
##   - p_e dips between the two main discontinuities where one burst = 1
## * Tidy code into functions before publication

clear

# Primary parameters
lambdaSum = 3/4;		# Utilization
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


# set qt_mode to true(1) to produce one time series of the queue
# set qt_mode to false(1) to scan parameter space and produce marking statistics
qt_mode = false(1);
if (qt_mode)
  i_lambda = 11; # index of lambda to plot if in qt_mode
  i_beta = 39;  # index of beta   to plot if in qt_mode
  i_phi = 4;    # index of phi    to plot if in qt_mode
endif

savepre = [mfilename()];
savem1 = ["Σλ", num2str(lambdaSum), ...
              "_Σβ", num2str(betaSum)];
savem2 = ["_β", num2str(betas)];
savesuf = ["_", strftime("%Y%m%d-%H%M%S", localtime (time ()))];

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

# p_stats(i,j,flow-id,approach,statistic), where
# i       : index of lambda
# j       : index of beta
# flow-id : 1,2 = flow a,b
# approach : 1,2 = soj,est
# statistic : 1,2 = mean,std-dev
p_stats = zeros(lambdaSum-1,betaSum-1,2,2,2);

if (!qt_mode)
  i_lambda = 1 : lambdaSum - 1;
  i_beta = 1 : betaSum - 1;
  i_phi = 1 : phis;
endif

for i = i_lambda
  if (!qt_mode)
    printf("%d%%\n", double(i)/double(lambdaSum)*100);
  endif
  for j = i_beta
    if (!qt_mode)
      printf(".");
      ##printf("%d,", j);
    endif
    # i_freq indexes the flow with more frequent bursts (or smaller burst size in
    #  case of a tie)
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
    
    p = zeros(phis,2,2);    # Marking prob for each flow & each approach 
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
    
    # k : index of phi
    for k = i_phi
      clear qt_out;
      ## Debug; ToDo: Remove once fixed
      i_bug = 0;
      
      # Variable definitions
      #  t :            current time
      #  i_head :       which flow is at the head of the queue
      #  t_burst[2] : start time of next burst from each flow
      #  i_next_burst : which burst is next
      
      # Time is scanned in two passes:
      # 1) to find where the q will be shortest
      # 2) to derive the actual queue process
      # Define the origin (t=0 & q=0) at the point when the combined queue is 
      #  shortest, assuming no standing queue, given utilization <= 100%
      
      # #1 time scan
      # q cannot be shorter than that after the max phase shift betw different
      #  bursts of different flows
      # Reason: In a cycle of duration t_max before the arrival pattern repeats,
      #  wherever is chosen as t=0 the same total amount of load arrives.
      #  So q will be shortest after the longest spread of arrivals.
      #  That is from a rare burst to the freq burst that is the shortest phase
      #  shift before the next rare burst.
      # Only rare bursts need to be scanned to find the min phase shift back to 
      #  the prev freq burst, 'cos the max phase shift has to be between bursts 
      #  of different types.
      t = 0;
      c = 0;                                   # rare burst counter
      # 'Seed' the initial phase shift using fraction phi of the freq interval
      t_delta0 = ti(i,j,i_freq) * (phi(k)+smidgen)/phis; # start at a freq burst
      # Vector of rare burst times
      t_delta = t_delta0 : ti(i,j,i_rare) : t_max(i,j);
      # Replace with vector of phase shifts from previous freq burst
      t_delta = rem(t_delta, ti(i,j,i_freq));
      t_delta_min = min(t_delta);
      clear t_delta;
      # Whether q=0 at the start or end of the min phase shift from freq to rare
      #  depends on whether the freq burst is large enough to keep the queue 
      #  busy over t_delta_min.
      # In each case, i_next_burst and i_head are pointed to the flow that 
      #  bursts at the origin and the time until the next burst for each flow, 
      #  t_burst[2], is set.
      if (t_delta_min < beta(i_freq,j))
        # q=0 at freq burst before min phase shift
        i_next_burst = i_head = i_freq;
        t_burst(i_freq) = 0;
        t_burst(i_rare) = t_delta_min;
      else
        # q=0 at rare burst after min phase shift
        i_next_burst = i_head = i_rare;
        t_burst(i_rare) = 0;
        t_burst(i_freq) = ti(i,j,i_freq) - t_delta_min;
      endif
      # Tie-break when bursts coincide: take smaller first
      if (t_delta_min <= 0)   # <= to be robust 
        [~, i_next_burst] = min(beta(:,j));
        i_head = i_next_burst;
        t_burst(i_rare) = t_burst(i_freq) = 0;
      endif
      t_next_burst = 0;
      
      # #2 time scan
      t = 0;
      q = zeros(2,1);     # queue delay contributed by each flow
      ott = 0;            # whether combined queue is over the threshold (q>=1)
      (qt_mode) && i_event = 0;      # Event index
      #
      # The qt_out matrix is filled as time is scanned:
      #  * qt_out(:,1) : time of next event; 
      #  * qt_out(:,2:3) : resulting per flow queue contribution at that time; 
      #  * qt_out(:,4) : which flow is at the head of the q; 0:1 means q(1:2);
      #  * qt_out(:,5) : whether queue is above threshold (q>=1) after t.
      while (t <= t_max(i,j))
      
        # The contributions from each flow to the queue are piecewise linear
        #  between 'events', where an 'event' is a discontinuity in one of the
        #  contributions or when the queue crosses the marking threshold
        # 
        t_next_empty = t + q(i_head); # Time the head q will next empty
        q_xs = sum(q) - 1;            # Excess q above threshold (could be -ve)
        # Check whether combined q falls below threshold 
        #  before head q empties or burst arrives
        if ( (ott) && ...
             ( ((t_next_empty <= t_next_burst) && (q_xs <= q(i_head))) || ...
               ((t_next_empty > t_next_burst) && (q_xs <= t_next_burst - t)) ...
             ) ...
           )
          # Combined queue has fallen below threshold before any other event
          t += q_xs;
          if (t > t_max(i,j))
            break
          endif
          # Add to p_e
          p(k, i_head, 2) += q_xs;
          q(i_head) -= q_xs;
          (qt_mode) && qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott];
          ott = 0;
          # The earlier condition (q_xs <= t_next_burst - t) could have been
          #  changed to < to suppress the following qt_out in the case when it
          #  seems redundant when q drains exactly to the threshold just
          #  before a burst. But, for robustness, if (q==1) ott=0, even if q
          #  is about to increase again.
          (qt_mode) && qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott];
          # Special case of q(i_head) emptying just as tail crosses threshold
          if (t >= t_next_empty)
            # Point i_head to other flow's queue if it's non-negative
            i_other = !(i_head-1) + 1;
            if (q(i_other) > 0)
              (qt_mode) && qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott];
              i_head = i_other;
              # No need to update i_other - not read elsewhere
            endif
          endif
        else
          if ((q(i_head) > 0) && (t_next_empty <= t_next_burst))
            # q(i_head) has emptied (defer any simultaneous burst to next event)
            t = t_next_empty;
            if (t > t_max(i,j))
              break
            endif
            if (ott)
              # Altho head flow has emptied, tail is over threshold
              #  so add time since previous event to p_e
              p(k, i_head, 2) += q(i_head);
            endif
            q(i_head) = 0;
            # Point i_head to other flow's queue if it's non-negative
            i_other = !(i_head-1) + 1;
            if (q(i_other) > 0)
              (qt_mode) && qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott];
              i_head = i_other;
              # No need to update i_other - not read elsewhere
            else
              # Both q's empty
              i_head = i_next_burst;
            endif
            (qt_mode) && qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott];
          else
            # Burst has arrived
            if (q(i_head) > 0)
              # Drain head queue since last event
              q(i_head) -= t_next_burst - t;  # resulting q(i_head) will be >0
              if (ott)
                # Add to p_e
                p(k, i_head, 2) += t_next_burst - t;
              endif
              q_xs = sum(q) - 1;
            endif
            t = t_next_burst;
            if (t > t_max(i,j))
              break
            endif
            (qt_mode) && qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott];
            # Add burst to tail
            #  but first check whether combined queue rises above threshold
            delta_q = beta(i_next_burst,j);
            if ( (ott == 0) && (q_xs + beta(i_next_burst,j) > 0) )
              # Take care! if (q_xs + beta == tiny), e.g tiny = 1.1102e-16
              #  after q incremented as below, q_xs = (sum(q) - 1) == 0;
              q(i_next_burst) -= q_xs;
              delta_q += q_xs;
              (qt_mode) && qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott];
              ott = 1;
            endif
            q(i_next_burst) += delta_q;
            if (ott)
              # Add to p_s
              p(k, i_next_burst, 1) += delta_q;
            endif
            (qt_mode) && qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott];
            # Prepare for next burst
            # Set t_burst for next burst from this flow
            t_burst(i_next_burst) += ti(i,j,i_next_burst);
            # Calc arrival time and flow id of next burst, handling tie if nec.
            [t_next_burst, i_next_burst] = min(t_burst);
            if (t_burst(1) == t_burst(2))
              [~, i_next_burst] = min(beta(:,j));
            endif
          endif    
        endif
      endwhile
      # Debug: End of each time series; ToDo: Remove once fixed
      if (sum(q) > 0)
        qt_bug(++i_bug,:) = [t, q(1), q(2), i_head-1, t_max(i,j), i, j, k];
      endif
    endfor
    
    p .*= lambdas/t_max(i,j);
    # In solely the following line, the abbreviated assignment operator (./=) 
    #  can't be used, because it errors within a script, 
    #  even though it is fine when run manually.
    p = p ./ cast(lambda(i,:),"double");
    
    if (!qt_mode)
      # p_stats(i,j,flow-id,approach,statistic), where
      # i       : index of lambda
      # j       : index of beta
      # flow-id : 1,2 = flow a,b
      # approach : 1,2 = soj,est
      # statistic : 1,2 = mean,std-dev
      p_stats(i,j,:,:,1) = mean(p(1:phis,:,:));
      p_stats(i,j,:,:,2) = std(p(1:phis,:,:),1);
    endif

  endfor
endfor

data_dir = "";
if (exist("octave_data", "dir"))
  data_dir = "octave_data/";
endif
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
  savefile = [savepre, "_qt_out_", savem1, ...
              "_", num2str(i_lambda), "λ", num2str(lambdas), ...
              "_", num2str(i_beta), "β", num2str(betas),...
              "_", num2str(i_phi), "φ", num2str(phis), ...
              savesuf];
  save("-binary", [data_dir savefile ".bin"], "savefile", "qt_out");
else
  savefile = [savepre, "_p_stats_", savem1, savem2, savesuf];
  save("-binary", [data_dir savefile ".bin"], "savefile", ...
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
## i_next_burst
##
## # Time intervals
## t
## t_delta0
## t_delta          [1, nnn]
## t_delta_min
## ti               [lambdaSum-1, betaSum-1,      2 ]
## t_max            [lambdaSum-1, betaSum-1 ]
## t_burst          [2  ]
## t_next_burst
## t_next_empty
##
## # Queue delays
## q_xs
## delta_q
##
## # Outputs
## q                [2  ]
## ott
## qt_out           [mmm,         8]
## p                [phis+1,      2,              2 ]
## p_stats          [lambdaSum-1, betaSum-1, 2, 2, 2]
##
## # Strings
## savefile
## savepre
## savem1
## savem2
## savesuf
## ============================================================================


