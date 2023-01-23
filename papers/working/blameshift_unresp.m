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
## * Investigate why the queue has brief empty periods in the following:
##    blameshift_unresp_qt_out_lambdaSum100_betaSum125_betas400_i1_j50_k1_20230109-154003
## * Check whether code completes properly at each t_max limit (scan #1 & #2)
## * Understand why 
##   - p_e dips between the two discontinuities
##   - qt_mode: plot at and either side of discontinuities 
## * Tidy code into functions before publication

clear

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
lambdas = 16;		# no. of steps of capacity share, lambda
betas = 400;		# no. of steps of normalized burst delay, beta
phis = 8;		# no. of steps of phase shift, phi, in 360deg

# set qt_mode to true(1) to produce one time series of the queue
# set qt_mode to false(1) to scan parameter space and produce marking statistics
qt_mode = true(1);
if (qt_mode)
  i_lambda = 1; # index of lambda to plot if in qt_mode
  i_beta = 50;  # index of beta   to plot if in qt_mode
  i_phi = 2;    # index of phi    to plot if in qt_mode
  if (i_lambda < 1 || i_lambda > lambdas)
    error("capacity share index parameter 'i_lambda' outside valid range");
  endif
  if (i_beta < 1 || i_beta > betas)
    error("burst size index parameter 'i_beta' outside valid range");
  endif
  if (i_phi < 1 || i_phi > phis)
    error("burst size index parameter 'i_beta' outside valid range");
  endif
endif

savepre = [mfilename()];
savemid = ["lambdaSum", num2str(lambdaSum*100), ...
              "_betaSum", num2str(betaSum*100), "_betas", num2str(betas)];
savesuf = ["_", strftime("%Y%m%d-%H%M%S", localtime (time ()))];

# Upscale primary parameters to integers
lambdaSum *= lambdas;
betaSum *= betas;
lambdaSum = cast(lambdaSum, "uint16");
betaSum = cast(betaSum, "uint16");

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
    printf(".");
    ##printf("%d,", j);
    # i_freq indexes the flow with more frequent bursts or smaller burst size in
    #  case of a tie
    # i_rare indexes the other flow
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
      # No q will be shorter than that after the max phase shift betw different
      #  bursts of different flows
      # Reason: Wherever is chosen as t=0, the same total amount of load arrives
      #  in a cycle of duration t_max before the pattern repeats.
      #  So q will be shortest after the longest time without any arrivals of
      #  the other flow, irrespective of which burst is last, because the queue
      #  depends on all the accumulated arrivals in the busy period before it 
      #  empties.
      #
      # Only freq bursts need to be scanned to find the max phase shift:
      #  * either back to the last rare burst 
      #  * or forward to the next rare burst
      #  'cos the max phase shift has to be between bursts of different types
      #
      t = 0;  ## ToDo: Not needed?
      c = 0;                                   # rare burst counter
      t_delta_max = zeros(1,2);                # max phase shift so far
      # 'Seed' the initial phase shift using fraction phi of the freq interval
      t_delta0 = ti(i,j,i_freq) * phi(k)/phis;
      while (t <= t_max(i,j))
        t = t_delta0 + c * ti(i,j,i_rare);     # time of each rare burst
        t_delta(1) = rem(t, ti(i,j,i_freq));   # shift back to last freq burst
        t_delta(2) = ti(i,j,i_freq) - t_delta(1); # shift fwd to next freq burst
        t_delta_max = max(t_delta, t_delta_max);
        ## ToDo: might be able to set i_next_burst here, sthg like:
        ## [t_next_burst, i_next_burst] = min(t_burst);
        c++;
      endwhile
      clear c;

      # Set defaults for case where t=0, q=0 at the start of a rare burst.
      # In each case, i_head is pointed to the flow that bursts at the origin
      #  and t_burst[2], the time until the next burst for each flow, is set.
      if ( t_delta_max(2) > t_delta_max(1) )
        # q=0 at start of first freq burst
        i_next_burst = i_head = i_freq;
        t_burst(i_freq) = 0;
        t_burst(i_rare) = ti(i,j,i_rare) - t_delta_max(2);    # >= 0
      else
        # q=0 at start of first rare burst (or t_delta_max == 0, i.e. sync'd)
        i_next_burst = i_head = i_rare;
        t_burst(i_rare) = 0;
        t_burst(i_freq) = ti(i,j,i_freq) - t_delta_max(1);    # >= 0
      endif
      # Special case when bursts coincide
      if (max(t_delta_max) == ti(i,j,i_freq))
        [~, i_next_burst] = min(beta(:,j));
        t_burst(i_rare) = t_burst(i_freq) = 0;
      endif
      t_next_burst = 0;
      ## ToDo: I think this is all predictable from which of the above cases we're in
      ## t_next_burst = 0   # in all case
      ## i_next_burst = i_head    # in all cases
      ##
      # Calc arrival time and flow id of next burst
##      [t_next_burst, i_next_burst] = min(t_burst);

      # If there's a tie, point i_next_burst to the index of the smaller burst
      #  (otherwise min() points i_next_burst to the lower index)
##      if (t_burst(1) == t_burst(2))
##        [~, i_next_burst] = min(beta(:,j));
##      endif
      
      # #2 time scan
      t = 0;
      i_event = 0;        # Event index
      q = zeros(2,1);     # queue delay contributed by each flow
      ott = 0;            # whether combined queue is over the threshold (q>=1)
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
        if ( (q_xs > 0) && ...
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
            if (qt_mode)
              qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott];
            endif
            ott = 0;
            # The earlier condition (q_xs <= t_next_burst - t) could have been
            #  changed to < to suppress the following qt_out in the case when it
            #  seems redundant when q drains exactly to the threshold just
            #  before a burst. But, for robustness, if (q==1) ott=0, even if q
            #  is about to increase again.
            if (qt_mode)
              qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott];
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
              if (qt_mode)
                qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott];
              endif
              i_head = i_other;
              # No need to update i_other - not used elsewhere
            endif
            if (qt_mode)
              qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott];
            endif
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
            if (qt_mode)
              qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott];
            endif
            # Add burst to tail
            #  but first check whether combined queue rises above threshold
            delta_q = beta(i_next_burst,j);
            if ( (ott == 0) && (q_xs + beta(i_next_burst,j) > 0) )
              q(i_next_burst) -= q_xs;
              delta_q += q_xs;
              if (qt_mode)
                qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott];
              endif
              ott = 1;
            endif
            q(i_next_burst) += delta_q;
            if (ott)
              # Add to p_s
              p(k, i_next_burst, 1) += delta_q;
            endif
            if (qt_mode)
              qt_out(++i_event,:) = [t, q(1), q(2), i_head-1, ott];
            endif
            
            # Prepare for next burst
            # Set t_burst for next burst from this flow
            t_burst(i_next_burst) += ti(i,j,i_next_burst);
            if (q(i_head) == 0)
              # Queue of head flow emptied, just as burst in other flow arrived
              i_head = i_next_burst;
            endif
            # Calc arrival time and flow id of next burst, handling tie if nec.
            [t_next_burst, i_next_burst] = min(t_burst);
            if (t_burst(1) == t_burst(2))
              [~, i_next_burst] = min(beta(:,j));
            endif
          endif    
        endif
      
      endwhile
    endfor
    
    p .*= lambdas/t_max(i,j);
    # In solely the following line, the abbreviated assignment operator (./=) 
    #  can't be used, because it errors within a script, 
    #  even though it is fine when run manually.
    p = p ./ cast(lambda(i,:),"double");
    # Append in the last 2 rows the mean & std-dev of the other rows
    
    # p_stats(i,j,flow-id,approach,statistic), where
    # i       : index of lambda
    # j       : index of beta
    # flow-id : 1,2 = flow a,b
    # approach : 1,2 = soj,est
    # statistic : 1,2 = mean,std-dev
    p_stats(i,j,:,:,1) = mean(p(1:phis,:,:));
    p_stats(i,j,:,:,2) = std(p(1:phis,:,:),1);

  endfor
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
  # Convert to stacked plot
  ## Commented out: no need because built in to area() function
  ##qt_out(:,2:4) = cumsum(qt_out(:,2:4), 2);
  savefile = [savepre, "_qt_out_", savemid, ... 
              "_i", num2str(i_lambda), ...
              "_j", num2str(i_beta), ...
              "_k", num2str(i_phi), ...
              savesuf];
  save("-binary", [savefile ".bin"], "savefile", "qt_out");
else
  savefile = [savepre, "_p_stats_", savemid, savesuf];
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
## 
## # Limits
## lambdaSum
## betaSum
## lambdas
## betas
## phis
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
## # Counters
## c
## 
## # Time intervals
## t
## t_delta0
## t_delta
## t_delta_max
## ti               [lambdaSum-1, betaSum-1,      2 ]
## t_max            [lambdaSum-1, betaSum-1 ]
## t_burst          [2  ]
## t_next_burst
## t_next_empty
##
## # Queue delays
## q_xs
## delta_q
## q_xs
##
## # Outputs
## q                [2  ]
## ott
## qt_out           [nnn,         8]
## p                [phis+1,      2,              2 ]
## p_stats          [lambdaSum-1, betaSum-1, 2, 2, 2]
##
## # Strings
## savefile
## savepre
## savemid
## savesuf
## ============================================================================


