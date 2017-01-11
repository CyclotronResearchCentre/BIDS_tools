function [fn_out,nr_out] = crc_BIDS_select(ffilt)
% Function to select images/data from a BIDS dataset based on some filters.
%
% The filters are defined by the fields of the 'ffilt' structure:
% FORMAT
%   [fn_out,nr_out] = crc_BIDS_select(ffilt)
% or
%   st_out = crc_BIDS_select(ffilt)
%
% INPUT
% ffilt : filtering details, expressed as a structure
%     .rootDir:     root directory path of BIDS dataset (path)
%     .SubjType:    type of subject to consider [def. '']
%     .SubjInd: 	index of subjects to consider or 'all' [def. 'all']
%     .SessInd:     index/name of session to consider [def. '']
%     .TaskLab:     label of task to consider [def. '']
%     .ImgMod:      name of imaging modality [def. 'bold']
%     .DatMod:      name of structured data field to return [def. 'events']
%     .ProcLec:     level of the data processing ('raw', 'derivative' or 
%                   'results') [def. 'raw']
%     .FnPrefx:     required prefix to filename [def. '']
%     .RegExp:      regular expression for file selection [def. '']
% 
% All key-names should follow BIDS nomenclature.
%
% OUTPUT
% fn_out : char array with all the requested (full path) filenames
% n_out  : number of such filenames returned
% or
% st_out : (array of) structure(s) with requested data
%
%__________________________________________________________________________
%
% BIDS (Brain Imaging Data Structure): http://bids.neuroimaging.io/
%   The brain imaging data structure, a format for organizing and
%   describing outputs of neuroimaging experiments.
%   K. J. Gorgolewski et al, Scientific Data, 2016.
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

%% 1. Initializing and some checks
% Ensures the BIDS complete structure is extracted only once
persistent BIDS 

% Defaults for the filter
ffilt_def = struct( ...
    rootDir, [], ...    % root directory path of BIDS dataset
    SubjType, [], ...   % type of subject to consider
    SubjInd, 'all', ... % index of subjects to consider
    SessInd, '', ...    % index of session to consider
    TaskLab, '', ...    % label of task to consider
    ImgMod, {{'bold'}}, ... % imaging modality
    DatMod, '', ...     % name of structured data field to return
    ProcLec, 'raw', ... % origin of the data (raw, derivative or results)
    FnPrefx, '', ...    % required prefix to filename
    RegExp, '' ...      % regular expression for file selection
    );

% Filling in defaults
ffilt = crc_check_flag(ffilt_def,ffilt);

% Loading in the whole BIDS stucture for the 1st call
if isempty(BIDS)
    BIDS = spm_BIDS(ffilt.rootDir);
end

%% #. Extracting the requested filenames/data
fn_out = '';

%% #. Preparing the output
% Number of filenames returned
if nargout==2
    nr_out = size(fn_out,1);
end

end