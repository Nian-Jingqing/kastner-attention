function setPaths( user )

switch user
  case 'fzhu23'
    % add paths for Feng
    addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
    pulvinarBase = '/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention';
    
  case 'cpandar'
    % add paths for Feng
    addpath('/snel/home/cpandar/c/feng_lfads-run-manager/src');
    pulvinarBase = '/snel/home/cpandar/c/Kastner_Attention';
end

addpath( fullfile( pulvinarBase, '/Matlab_Offline_Files_SDK/') );
addpath( fullfile( pulvinarBase, '/myTools/kastner_analysis_tools') )
addpath( fullfile( pulvinarBase, '/myTools') )
addpath( fullfile( pulvinarBase, '/myTools/rawDataProcessing') )
addpath( fullfile( pulvinarBase, '/codes') )
addpath( fullfile( pulvinarBase, '/codes/broadbandProcessing') )
addpath( fullfile( pulvinarBase, '/codes/postAnalysisCodes') )
% jpca paths
addpath( fullfile( pulvinarBase, '/myTools/jPCA_tools') )
addpath( fullfile( pulvinarBase, '/myTools/jPCA_tools/fromMarksLibraries') )
addpath( fullfile( pulvinarBase, '/myTools/jPCA_tools/CircStat2010d') )
