### Research Question
How to apportion blame for a queue with arrivals in different bursts?

### Approach
Comparison of different ECN marking algorithms applied to two unresponsive flows 
with different burst amplitudes and different regular intervals between bursts. 
See the overview presentation of the results (below) for explanation and justification of the approach.

### Results
* Overview: "[How to Apportion Blame for a Queue with Arrivals in Bursts?](https://github.com/bbriscoe/l4s-aqm/blob/master/presents/sigqdyn-overview.pdf)"
* Detailed results "Signalling Bursty Queue Delay; Expt 1.1" using one of a choice of four metrics:
  * marking probability of flow a, [_p<sub>a</sub>_](https://github.com/bbriscoe/l4s-aqm/blob/master/presents/sigqdyn-p-expt1_1.pdf)
  * difference between marking probabilities of flow a and b, [_Δp_](https://github.com/bbriscoe/l4s-aqm/blob/master/presents/sigqdyn-Δp-expt1_1.pdf)
  * marking rate of flow a, [_λp<sub>a</sub>_](https://github.com/bbriscoe/l4s-aqm/blob/master/presents/sigqdyn-λp-expt1_1.pdf)
  * difference between marking rates of flow a and b, [_Δλp_](https://github.com/bbriscoe/l4s-aqm/blob/master/presents/sigqdyn-Δλp-expt1_1.pdf)
* The technical report on this research: "Rapid Signalling of Queue Dynamics"<br>
\[ [Latest work-in-progress version](https://github.com/bbriscoe/l4s-aqm/blob/master/papers/sigqdyn_tr.pdf) | [Latest published version](https://arxiv.org/abs/1904.07044) \]<br>
(neither are necessarily as up to date as the Overview presentation above).

### Usage
To repeat the experiments:
 1. Edit primary and secondary parameters directly in the file [`blameshift_unresp.m`](https://github.com/bbriscoe/l4s-aqm/blob/master/papers/working/blameshift_unresp.m) and save.
 2. Run either from within [octave](https://octave.org/) as<br>
      `>> blameshift_unresp`<br>
    or from the command line as<br>
      `$ octave blameshift_unresp.m`
 3. A binary (.bin) file of the results is saved at the end, in subdirectory `octave_data` if it exists, otherwise in the current directory.
 4. To plot output at any time, this saved binary can be loaded into octave:<br>
     `>> load <file.bin>`<br>
     But this is unnecessary if plotting is done immediately after running this script, while the output is still in memory.
 5.  Then a second script must be run from octave, depending on how `qt_mode` was set in the original script:<br>
      `>> blameshift_unresp_plot_qt     `# time series plot (requires `qt_mode = true(1)`)<br>
      `>> blameshift_unresp_plot_stats  `# 2-D plot of parameter space (requires `qt_mode = false(1)`)<br>
      `>> blameshift_unresp_mesh_stats  `# 3-D plot of parameter space (requires `qt_mode = false(1)`)<br>
     See the source of each of these scripts for their usage.

### Testing
* octave 5.2.0 under Ubuntu 20.04 LTS
* octave 7.3.0 under Win 10

Not tested from MatLab
