# Scan Blame Shifting problem space

## ToDo:
## * Debug (i,j,k) = (17,2,5) underutilization
## * Test a few other (i,j,k)
## * averaging over all phase angles
## * choosing which points in the parameter space to output
## * for loops
## * pretty plotting

clear 

## ============================================================================
## # Manifest of variables
## # =====================
## # Scalars
## # _______
## 
## # Limits
## lambdaSum
## betaSum
## 
## # Steps
## lambdas
## betas
## phis
## 
## # Pointers
## i_event
## i_freq
## i_rare
## i_head
## i_other
## i_burst_next
## 
## # Times
## t
## t_burst_next
## 
## # Matrices
## # ________
## 
## # Input Parameters
## lambda[lambdaSum-1, 2]
## beta[2, betaSum-1]
## phi[phis]
## 
## # Time intervals
## t_max[lambdaSum-1, betaSum-1]
## ti[lambdaSum-1, betaSum-1, 2]
## 
## q[2]
## t_burst[2]
## 
## output[???] {
## t,
## q[2],
## head
## }
## ============================================================================

# Primary parameters
lambdaSum = 1;		# Utilization
betaSum = 5/4;		# Max combined normalized burst delay (wrt marking threshold)
if (lambdaSum > 1)
  error("utilization parameter 'lambdaSum' cannot exceed 100%");
endif

# Secondary parameters
# lambda, beta & phi are held as integers 
#  with implicit denominators lambdas, betas & phis, resp.
#  (except after setup, when beta is cast to double and downscaled), 
lambdas = 32;		# no. of steps of capacity share, lambda
betas = 16;		# no. of steps of normalized burst delay, beta
phis = 8;		# no. of steps of phase shift, phi, in 360deg

# Upscale primary parameters to integers
lambdaSum *= lambdas;
betaSum *= betas;
lambdaSum = cast(lambdaSum, "uint16");
betaSum = cast(betaSum, "uint16");

# Fill vectors for scanning the problem space
# Two flows, a&b are represented by slices 1&2
# The second flow's capacity share and burst size use up the remainder of
#  lambdaSum and betaSum
lambda = linspace(1, lambdaSum - 1, lambdaSum - 1);
lambda = [lambda; lambdaSum - lambda]';
beta = linspace(1, betaSum - 1, betaSum - 1);
beta = [beta; betaSum - beta];
phi = linspace(0, phis - 1, phis);
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
ti(:,:,1) = double(lambdas/betas * beta(1,:)) ./ double(lambda(:,1));
ti(:,:,2) = double(lambdas/betas * beta(2,:)) ./ double(lambda(:,2));

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
## lambda = [31/32;1/32] beta = [2/16,18/16]
## t_max = lcm(16/32*1/32*2/16, 16/32*31/32*18/16) / (16*32*31/32*1/32)
##       = lcm(1*2,31*18) * 2 / (31*1)

ti_tmp(:,:,1) = beta(1,:) .* lambda(:,2);
ti_tmp(:,:,2) = beta(2,:) .* lambda(:,1);
t_max = lcm(uint16(ti_tmp(:,:,1)), uint16(ti_tmp(:,:,2)));
if (lambdas >= betas)
  t_max *= lambdas/betas;  # multiply before saturation check
endif
clear ti_tmp;
if (max(max(t_max)) == intmax("uint16"))
  error("lcm() saturated at least one uint16 element of 't_max'");
endif
# Cast t_max to FP and downscale
t_max = double(t_max) ./ double(lambda(:,1) .* lambda(:,2));
if (lambdas < betas)
  t_max *= lambdas/betas;  # multiply after saturation check
endif
# Cast beta to FP and downscale
beta = double(beta) ./ betas;

## For fast testing, temporarily pick test values of i,j
#  ToDo: (to be replaced by nested for loops)
i = 17; # index of lambda
j = 2;  # index of beta

# i_freq indexes the more frequent slice 
#  (i.e. the flow with lower interval between bursts)
#  strictly i_freq is the slice that is not less frequent, i.e. might be same
# i_ indexes the other flow
[~, i_freq] = min(ti(i,j,:));
i_rare = !(i_freq-1) + 1;
## ToDo: catch case when both the same interval

# As long as utilization < 100% and bursts are evenly spaced, 
#  there cannot be >1 sequence of packets from each flow in the queue at once
# So a vector of 2 queue sizes and a 2-state integer 'i_head' indexing the flow
#  at the head are sufficient to fully describe the queue at any one time

## For fast testing, temporarily pick a test value of k
# ToDo: replace with for k = 1 : phis
k = 5;  # index of phi

# Variable definitions
#  t :            current time
#  i_head :       which flow is at the head of the queue
#  t_burst[2] : start time of next burst from each flow
#  i_burst_next : which burst is next

# Define t=0 & q=0 at the point when the combined queue is smallest
#  assuming no standing queue, given utilization <= 100%
# Each rare burst is scanned to find the max phase shift back to the last freq 
#  burst that is still draining (if any). Then, at the start of this freq burst,
#  t=0, q=0, t_burst(i_freq)=0.
#  Then t_burst(i_rare) should be set to the max phase shift, t_delta_max.
# If no rare bursts arrive at a non-empty queue, at the start of any rare burst 
#  t=0, q=0, t_burst(i_rare)=0.
# Depending on which case, i_head is pointed to the flow that bursts at t=0
t = 0;
c = 0;                                         # rare burst counter
t_delta0 = ti(i,j,i_freq) * (1 - phi(k)/phis); # initial phase shift
t_delta_max = 0;                               # max overlapping phase shift
t_burst(i_freq) = t_delta0;
i_head = i_rare;
while (t <= t_max(i,j))
  t = t_delta0 + c * ti(i,j,i_rare);           # time of each rare burst
  t_delta = rem(t, ti(i,j,i_freq));            # phase shift at each rare burst
  if ( (t_delta > t_delta_max) && (t_delta <= beta(i_freq,j)) )
    t_delta_max = t_delta;
  endif
  c++;
endwhile
t_burst(i_rare) = t_delta_max;
if (t_delta_max > 0)
  # queue empty at start of first i_freq burst (or of both)
  t_burst(i_freq) = 0;
  i_head = i_freq;
endif

# Now run through time again to calculate the queue
t = 0;
i_event = 0;    # Event index
q = zeros(2,1); # queue delay contributed by each flow
while (t <= t_max(i,j))

  # Calc arrival time and flow id of next burst
  [t_burst_next, i_burst_next] = min(t_burst);
  # If there's a tie, point i_burst_next to the index of the smaller burst
  #  (otherwise min() points i_burst_next to the lower index)
  if (t_burst(1) == t_burst(2))
    [~, i_burst_next] = min(beta(:,j));
  endif
  
  # The contributions from each flow to the queue are piecewise linear between
  # 'events', where an 'event' is a discontinuity in one of the contributions
  # Output matrix:
  #  * output(:,1) : time of next event; 
  #  * output(:,2:3) : resulting per flow queue contribution at that time; and 
  #  * output(:,4) : zero (placeholder)
  #  * output(:,5) : which flow is at the head of the queue; 0:1 means q(1:2).
  if ((q(i_head) > 0) && (t + q(i_head) <= t_burst_next))
    # q(i_head) has emptied (defer any simultaneous burst to next event)
    t += q(i_head);
##    if (t > t_max(i,j))
##      break
##    endif
    q(i_head) = 0;
    # Point i_head to other flow's queue if it's non-negative
    i_other = !(i_head-1) + 1;
    if (q(i_other) > 0)
      output(++i_event,:) = [t, q(1), q(2), 0, i_head-1];
      i_head = i_other;
    endif
    output(++i_event,:) = [t, q(1), q(2), 0, i_head-1];
  else
    # Burst arrives
    if (q(i_head) > 0)
      # Drain head queue since last event
      q(i_head) -= t_burst_next - t;  # resulting q(i_head) will be positive
    endif
    t = t_burst_next;
##    if (t > t_max(i,j))
##      break
##    endif
    output(++i_event,:) = [t, q(1), q(2), 0, i_head-1];
    # Add burst to tail
    q(i_burst_next) += beta(i_burst_next,j);    
    output(++i_event,:) = [t, q(1), q(2), 0, i_head-1];
    # Set t_burst for next burst from this flow
    t_burst(i_burst_next) += ti(i,j,i_burst_next);
    if (q(i_head) == 0)
      # Queue of head flow emptied, just as burst in other flow arrived
      i_head = i_burst_next;
    endif
  endif

endwhile
# ToDo: endfor
