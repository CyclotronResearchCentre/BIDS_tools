% 
% Bits of code to test the crc_BIDS_select.m functionalities.
% 
% Rely on the example BIDS data set from:
%   https://github.com/INCF/BIDS-examples
% 
% FILTER options
% - ffilt     : filtering details, expressed as a structure
% #subjects
%     .rootDir:     root directory path of BIDS dataset [def. current dir]
%     .SubjType:    type of subject to consider [def. 'all']
%     .SubjInd: 	index of subjects to consider or 'all' [def. 'all']
% #images, sessions & runs
%     .SessInd:     index/name of session to consider [def. ''=all]
%     .RunInd:      index of runs to consider [def. []=all]
%     .TaskLab:     label of task to consider [def. ''=all]
%     .ImgMod:      name of imaging modality [def. 'func']
%     .ImgType:     image type for structural [def. ''=all]
% #others
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
    % Put here the data from the example repo on GItHub, cf. help.
    rDir = '/Users/cphillips/Dropbox/Work/3_data/BIDS-examples-master/';
end

% Load everything and check 
[BIDSall,lDir] = crc_BIDS_checkDS(rDir);


% Just pick one that works, e.g. 3rd directory
% 2 func sessions for 3 tasks.
idat = 3;
B = BIDSall{idat};

% Get functional from all subject, all runs, for a single task
ffilt = struct( ...
    'ImgMod', {{'func'}}, ...
    'TaskLab', {{'mixedeventrelatedprobe'}} );
[fn_out,nr_out] = crc_BIDS_select(ffilt,B)

% Get functional from all subject, all runs, for 2 tasks
ffilt = struct( ...
    'ImgMod', {{'func'}}, ...
    'TaskLab', {{'mixedeventrelatedprobe','probabilisticclassification'}} );
[fn_out,nr_out] = crc_BIDS_select(ffilt,B)

% Get functional from all subject, for 1 task and 2nd run
ffilt = struct( ...
    'ImgMod', {{'func'}}, ...
    'TaskLab', {{'deterministicclassification'}}, ...
    'RunInd', 2);
[fn_out,nr_out] = crc_BIDS_select(ffilt,B)

% Get functional from all subject, for 2 tasks and 1st run
ffilt = struct( ...
    'ImgMod', {{'func'}}, ...
    'TaskLab', {{'mixedeventrelatedprobe','probabilisticclassification'}}, ...
    'RunInd', 1);
[fn_out,nr_out] = crc_BIDS_select(ffilt,B)

% Get functional from subject #4 & #6, for all tasks and 1st run
ffilt = struct( ...
    'ImgMod', {{'func'}}, ...
    'SubjInd', [4 6], ...
    'TaskLab', '', ...
    'RunInd', 1);
[fn_out,nr_out] = crc_BIDS_select(ffilt,B)

% Get anatomical from all subject, all types
ffilt = struct( ...
    'ImgMod', {{'anat'}} );
[fn_out,nr_out] = crc_BIDS_select(ffilt,B)

% Get anatomical from all subject, only T1w
ffilt = struct( ...
    'ImgMod', {{'anat'}} , ...
    'ImgType', {{'T1w'}} );
[fn_out,nr_out] = crc_BIDS_select(ffilt,B)

% Get anatomical from all subject, both T1w & inplaneT2
ffilt = struct( ...
    'ImgMod', {{'anat'}} , ...
    'ImgType', {{'T1w','inplaneT2'}} );
[fn_out,nr_out] = crc_BIDS_select(ffilt,B)

% Get event from subject #5, for 1 task and 2nd run
ffilt = struct( ...
    'ImgMod', {{'func'}}, ...
    'SubjInd', [5 7], ...
    'TaskLab', {{'deterministicclassification'}}, ...
    'RunInd', 2, ...
    'DatType', 'val', ...
    'DatField', {{'events'}});
[v_out,nr_out] = crc_BIDS_select(ffilt,B)