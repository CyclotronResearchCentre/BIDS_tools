% 
% Bits of code to test the crc_BIDS_select.m functionalities.
% 
% Rely on the example BIDS data set from:
%   https://github.com/INCF/BIDS-examples
% 
% FILTER options
% - ffilt     : filtering details, expressed as a structure
%     .rootDir:     root directory path of BIDS dataset [def. current dir]
%     .SubjType:    type of subject to consider [def. 'all']
%     .SubjInd: 	index of subjects to consider or 'all' [def. 'all']
%     .SessInd:     index/name of session to consider [def. '']
%     .TaskLab:     label of task to consider [def. '']
%     .ImgMod:      name of imaging modality [def. 'bold']
%     .DatMod:      name of structured data field to return [def. 'events']
%     .ProcLev:     level of the data processing ('raw', 'derivative' or
%                   'results') [def. 'raw']
%     .FnPrefx:     required prefix to filename [def. '']
%     .RegExp:      regular expression for file selection [def. '']
%     .ResetBIDS:   force the reload the BIDS structure [def. false]
% 
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

% Where are the example data?
if ispc
    rDir = '';
else
    rDir = '/Users/cphillips/Dropbox/Work/3_data/BIDS-examples-master/';
end

% Load everything and check 
[BIDSall,lDir] = crc_BIDS_checkDS(rDir);


% Just pick one that works, e.g. 3rd directory
% 2 bold sessions for 3 tasks.
idat = 3;
B = BIDSall{idat};

% Get BOLD from all subject, all runs, for a single task
ffilt = struct( ...
    'ImgMod', {{'bold'}}, ...
    'TaskLab', {{'mixedeventrelatedprobe'}} );
[fn_out,nr_out] = crc_BIDS_select(ffilt,B)

% Get BOLD from all subject, all runs, for 2 tasks
ffilt = struct( ...
    'ImgMod', {{'bold'}}, ...
    'TaskLab', {{'mixedeventrelatedprobe','probabilisticclassification'}} );
[fn_out,nr_out] = crc_BIDS_select(ffilt,B)


