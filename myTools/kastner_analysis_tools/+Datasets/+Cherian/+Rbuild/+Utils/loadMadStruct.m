function madstruct=loadMadStruct(filename,loadmode); % load a MADF filename, 
% file can be loaded
% without loading all the data or the spike waveforms 
%
% usage:
% >> madstruct=loadmad('myfile')      % load the whole thing
% >> madstruct=loadmad('myfile',1)    % load header only
% >> madstruct=loadmad('myfile',2)    % load all except spike waveforms
%
%
% filename: duh
% floaddata: T/F flag for loading entire file (1) or just header (0)
% 

% 2004-09-16 v0.999
% can now load whole thing, header only, or all but spikewaveforms.
%
% known issues:if you have custom fields in there, only loadmode==0 will load them. 
% workaround: load them yourself, then assign them to madstruct.

% 2004-09-08 v0.99
% refining. 

% 2004-09-05 v0.9
% working on this, trying to make it do something clever
% loads, and converts the top level info field in the file back into a
% field inside the appropriate struct array for that set of signals.

% Since signal data isn't modified
% here, it shouldn't be copied in memory (matlab uses 'lazy' pass by value)
%
version=.999;

switch nargin
    case 1
        loadmode=0; % load everything
    case 2
        % umm, don't do anything, this is fine.
    otherwise
        disp('You need 1 or 2 args or something');
        disp('madstruct=loadmad(filename,loadmode)');
        madstruct=[];
        return
end % switch

if isempty(intersect(loadmode,[ 0 1 2]))
    disp('I don''t know what you want me to do');
    disp('usage:');
    disp('madstruct=loadmad(''myfile'')      % load the whole thing');
    disp('madstruct=loadmad(''myfile'',1)    % load header only');
    disp('madstruct=loadmad(''myfile'',2)    % load all except spike waveforms');
end

if ~(exist(filename,'file')) % hope it's a file
    disp(sprintf('%s ain''t no file I ever heard of',filename));
    madstruct=0;
    return
end % if file exists

load(filename,'savedmadversion'); % just load that to make sure it's there
switch 1 % do the first true one of these:
    case ~exist('savedmadversion')
        disp('there may be a problem. loadmad() can only load mad files saved with savemad()');
    case savedmadversion~=version
        disp('different version found. this could be interesting.');
    case savedmadversion<1
        disp('INFO: MadStruct Loader v1.0');
    otherwise
        disp('This is a perfectly fine version');
end % switch

sanctionedfields={'header', 'analog', 'levels', 'events', 'spikes'};
splitfields={'analog', 'levels', 'events', 'spikes'};

% --- now, load stuff:
switch loadmode
    
    % load the whole thing
    case 0
        disp('INFO: Loading File...');
        madstruct=load(filename);     % load whole thing, plus xtras? but I don't want that savedmadversion thing in there.
        
        madstruct=rmfield(madstruct,{'savedmadversion'}); % I refuse to believe that this will copy the structure in memory
        if isfield(madstruct,'README') % won't be there in older versions
            madstruct=rmfield(madstruct,{'README'}); % I refuse to believe that this will copy the structure in memory
        end
        disp('INFO: File load complete.')

    % header only, even though other stuff might be in there
    case 1
        madstruct=load(filename,'header', 'analoginfo', 'levelsinfo', 'eventsinfo', 'spikesinfo');
        
    % header & spikes, no waveforms. so everything except 'spikes'
    case 2
        madstruct=load(filename,'header', 'analoginfo', 'levelsinfo', 'eventsinfo', 'spikesinfo', ...
            'analog', 'levels', 'events', 'spikestimestamps');
end % switch

% because matlab is being stupid about the events field - Im temporarlily
% using the kluge below to make things work - AC 8-3-05

% apparently some files do not have this problem - so we shall build in a
% routine to check for this at some point - until then , we simply comment 
% out the kluge code until we have a problem ...so the program doesnt
% crash when it cant find it.

% madstruct.events=madstruct.events_;
% madstruct = rmfield(madstruct, 'events_');


% umm, convert info fields back into the main fields
% this happens for all loadmodes


for fi=1:length(splitfields)
    %madstruct.(splitinfofields{fi}).info=%struct();
    for si=1:length(madstruct.([splitfields{fi} 'info'])) % go through each element in the signal struct array
        madstruct.(splitfields{fi})(si).info=madstruct.([splitfields{fi} 'info'])(si);
    end
    madstruct=rmfield(madstruct,[splitfields{fi} 'info']);
end % fi

%if we load the data, spiketimestamps need to be put back
if (loadmode==0 | loadmode==2) % oh yeah, 
    for si=1:length(madstruct.spikestimestamps)
        madstruct.spikes(si).timestamps=madstruct.spikestimestamps{si};
    end % for si
    madstruct=rmfield(madstruct,'spikestimestamps'); % again, I hope this isn't too stupid
end % if floaddata

allfields=fieldnames(madstruct);
unsanctionedfields=setdiff(allfields,sanctionedfields);
if ~isempty(unsanctionedfields)
    disp(sprintf('with %d unknown fields',length(unsanctionedfields)));
    disp(unsanctionedfields');
end

end % loadmad