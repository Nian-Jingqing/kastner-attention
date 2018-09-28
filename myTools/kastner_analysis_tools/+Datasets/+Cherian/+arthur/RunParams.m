classdef RunParams < LFADS.RunParams
   properties
       % Add additional custom parameters here. The default you assign to
       % them will be used when computing the hash value. Any params whose value
       % differs from the default will be included in the hash value, to allow new
       % parameters to be added without invalidating old hashes. So choose
       % the default once and don't change it. If you decide to use another
       % value later by default, override it in the constructor instead.
       doPBT = false
       PBTscript = ''
       % Experiment block selector
       i_block = 1;
       % 'center-out' or 'random-walk'
       dataset_type = 'center-out';
       % 'align-chop' or 'overlap-chop'
       trialize_method = 'align-chop';
       % 'spikes' or 'emg'
       run_type = 'spikes';
       % *center-out options
       ap_threshold = 0.25;
       ap_max_speed_min = 10;
       pre_ap_ms = 200;
       post_ap_ms = 700;
       % *random-walk options*
       trial_time_ms = 500;
       trial_olap_ms = 100;
   end
   
   methods
       function par = RunParams()
          % adjust whatever defaults you like. Modificiations made here
          % function identically to modifications made by the user, meaning
          % that the value will be incorporated in the hash only if the
          % value differs from the property definition in the class.
          
       end 
   end
end