function [SR]=get_raster_single(trigs, spikes,ppms,time_before_trig,time_after_trig)
h=(reshape(spikes,1,[])-reshape(trigs,[],1))/ppms;
[trial,time]=find(h>=-time_before_trig & h<=time_after_trig);
time2=h(h>=-time_before_trig & h<=time_after_trig);
SR=[time2,trial];
