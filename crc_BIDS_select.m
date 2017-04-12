function [fn_out,nr_out] = crc_BIDS_select(ffilt,BIDS_spm)
% Function to select images/data from a BIDS dataset based on some filters.
%
% The filters are defined by the fields of the 'ffilt' structure:
% FORMAT
%   [fn_out,nr_out] = crc_BIDS_select(ffilt,BIDS_spm)
% or
%   st_out = crc_BIDS_select(ffilt,BIDS_spm)
%
% INPUT
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
% - BIDS_spm  : a BIDS-structure as extracted with spm_BIDS (see function 
%               in recent SPM12 distribution) or path name to it (this 
%               overlaods the 'rootDir' in the ffilt input!!!)
%
% All key-names should follow BIDS nomenclature.
%
% OUTPUT
% fn_out : char array with all the requested (full path) filenames
% n_out  : number of such filenames returned
% or
% st_out : (array of) structure(s) with requested data
%
% Notes regarding the filtering options
% * SubjType:
%   - if omited, then all types are considered
%   - if only one type, pass it as a char
%   - for multiple types, use a celle array of chars.
%     For example ffilt.SubjType = {{'ctl','pat'}} for 'ctl' and 'pat'.
% * SubjInd:
%   - if omited or set to 'all' (or another char), then no filtering
%   - if (array of) number(s), then keeping only those matching
% * ResetBIDS:
%   after the first call to crc_BIDS_select, the BIDS structure is saved
%   as a persistent variable. In order to reset this structure, i.e. for
%   the reload, then set this flag to 'true'. This would mainly be used if
%   data have changed on disk (unlikely but...) or when accessing another
%   BIDS data set during the same Matlab session (very possible).
%
% EXAMPLES:
% Use the cloned example data sets from
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
if nargin<1, ffilt = struct; end
ffilt_def = struct( ...
    'rootDir', pwd, ...   % root directory path of BIDS dataset
    'SubjType', {{''}}, ...   % type of subject to consider
    'SubjInd', 'all', ... % index of subjects to consider
    'SessInd', '', ...    % index of session to consider
    'TaskLab', '', ...    % label of task to consider
    'ImgMod', {{'bold'}}, ... % imaging modality
    'DatMod', '', ...     % name of structured data field to return
    'ProcLev', 'raw', ... % origin of the data (raw, derivative or results)
    'FnPrefx', '', ...    % required prefix to filename
    'RegExp', '', ...     % regular expression for file selection
    'ResetBIDS', false ...% force the reload the BIDS structure or not
    );
% Filling in defaults
ffilt = crc_check_flag(ffilt_def,ffilt);

% Load the BIDS structure
if ffilt.ResetBIDS, BIDS = []; end % Reset the BIDS-structure
if nargin<2
    % Loading in the whole BIDS stucture for the 1st call
    if isempty(BIDS)
        BIDS = spm_BIDS(ffilt.rootDir);
    end
else
    if isstruct(BIDS_spm)
        BIDS = BIDS_spm;
    elseif ischar(BIDS_spm)
        ffilt.rootDir = BIDS_spm;
        if isempty(BIDS)
            BIDS = spm_BIDS(ffilt.rootDir);
        end
    end
end

% Note:
% spm_BIDS does check for the existence and validity of the BIDS-directory
% provided, no need to redo this check here then?

% Check number of subject's dirs match participants info
if nrSubj~=numel(BIDS.participants.participant_id)
    error('Inconsistent number of subjects in BIDS directory.');
end

%% Extracting the requested filenames/data
fn_out = ''; %#ok<*NASGU>
nrSubj = numel(BIDS.subjects);

% Extracting the subjects' type and index from the original their label
[pTyp,pInd] = extract_participant_label_ind(BIDS.participants.participant_id);

% Find the list of subjects based on SubjType and SubjInd filter
l_Subj = 1:nrSubj;
if ~isempty(ffilt.SubjType)
    if isempty(ffilt.SubjType) || strcmp(ffilt.SubjType,'all') || ...
            isempty(ffilt.SubjType{1}) || strcmp(ffilt.SubjType{''},'all')
        % No filtering needed
    else
        % reduce list to those with matching (each) SubjType
        fSubjType = ffilt.SubjType;
        if ischar(fSubjType), fSubjType = cellstr(fSubjType); end
        to_keep = [];
        for ii=1:numel(fSubjType)
            to_keep = [to_keep ; find(strcmp(fSubjType{ii},pTyp))];
        end
        l_Subj = l_Subj(to_keep);
    end
    
end

if ischar(ffilt.SubjInd) % Not reducing list of chars just letting know
    if ~strcmp(ffilt.SubjInd,'all') % that 'all' should be used.
        warning('Wrong subject index, using ''all''.');
    end
elseif ~any(isnan(pInd)) % no NaN's -> apply index filter
    tmp = intersect(ffilt.SubjInd,pInd(l_Subj));
    to_keep = [];
    for ii = tmp'
        to_keep = [to_keep ; find(pInd(l_Subj)==ii)]; %#ok<*AGROW>
    end
    l_Subj = l_Subj(to_keep);
end

% Selecting files!
fn_out = '';
for ii = l_Subj
    for kk = 1:numel(ffilt.ImgMod)
        switch ffilt.ImgMod{kk}
            case 'bold'
                % Use BIDS structure & filter on 'TaskLab'
                ffilt_func = struct(...
                    'TaskLab',ffilt.TaskLab );
                tmp = get_bold_data(BIDS.subjects(ii).func,ffilt_func);
% NOTE:
% this needs to be extended in order to account for other filtering option
% such as 'AcqLab', 'RunLab' and 'RecLab', the acquisition, reconstruction 
% and run flags.

%                 tmp = fullfile(BIDS.subjects(ii).path, 'func', ...
%                     BIDS.subjects(ii).func.filename);
            otherwise
                warning('Unknown image type.')
                tmp = [];
        end
        if ~isempty(tmp)
            fn_out  = char(fn_out,tmp);
        end
    end
end



% Note:
% selecting files or data (array/structure) should be placed in separate
% subfunctions to clarify the code
% -> do that when dealing with data stuff

%% #. Preparing the output
% remove 1st line that is empty
if size(fn_out,1)>1, fn_out(1,:) = []; end

% Number of filenames returned
if nargout==2
    nr_out = size(fn_out,1);
end

end

% ________________________________________________________________________
%
%% SUBFUNCTIONS
% ________________________________________________________________________

function [pTyp,pInd] = extract_participant_label_ind(pLabel)
% Function to extract/decompose the participants' label (pLabel, char)
% into a 'type' (pTyp, char) and 'index' (pInd, integer).
%
% Principle:
% We assume that the index is placed at the end of the subject's label,
% i.e.the last continuous set of numbers (0 to 9) will converted as the
% subject's index (-> pInd), the 1st part will be its type (-> pTyp).
% The notion of last numbers is useful to handle the case where a number is
% used in the subject's label itself, e.g. the case of a disease stage.
% Afterwards, the list of indexes is checked. If a few are exactly the
% same (for the same pTyp!), then these are not representing a numeral
% index, e.g. the labels end with a number like a few subjects in 'Stage3'/
% 'Stage2'. So the indexes are removed (pInd = []) and the char labels
% (pTyp) are returned.
%
% INPUT
% pLabel : single subject label (char), char array or cell array of char
%
% OUPUT
% pTyp : cell array with the char part of the pLabel, i.e. subject's type
% pInd : array with the index part of the pLabel
%
% For example:
%   sub-ctrl012 -> label = 'ctrl012' -> pTyp = 'ctrl' and pInd = 12
%   it would return arrays in case of an array of labels.

% Turn into cell array of char
if ischar(pLabel), pLabel = cellstr(pLabel); end
nrLabel = numel(pLabel);

pTyp = cell(nrLabel,1);
pInd = zeros(nrLabel,1);

% All pLabel should start with 'sub-'
lCheck = cell2mat(strfind(pLabel,'sub-'));
if nrLabel~=numel(lCheck) || ~all(lCheck==1)
    error('Subject/participants label not fitting BIDS specs.')
end

% Extract pInd and pTyp values
for ii=1:nrLabel
    % split the label one at a time.
    tmp = pLabel{ii};
    lab_ii = tmp(5:end);
    l_numbers = regexp(lab_ii,'[0123456789]');
    if isempty(l_numbers) || l_numbers(end)~=numel(lab_ii)
        % no numbers or no numbers in last position -> pInd =  0;
        pInd(ii) = NaN;
        pTyp{ii} = lab_ii;
    else
        % Find last set of continuous numbers
        dl_numbers = diff(l_numbers);
        jj = 1; % we know we end with at least 1 number
        while ~isempty(dl_numbers) && dl_numbers(end)==1
            jj=jj+1; dl_numbers(end)=[];
        end
        % -> subject index (pInd) and the rest goes into pTyp
        pInd(ii) = str2double(lab_ii(end-jj+1:end));
        if jj<numel(lab_ii)
            pTyp{ii} = lab_ii(1:end-jj);
        else
            pTyp{ii} = '';
        end
    end
end
% All pInd should now be either a number or NaN.

% Check pInd's values, numbers or NaN's, according to pTyp
u_pTyp = unique(pTyp);
nu_pTyp = numel(u_pTyp);
ok_ind = zeros(nu_pTyp,1);
for ii=1:numel(u_pTyp)
    % l_ii = find(strcmp(u_pTyp{ii},pTyp));
    % ok_ind = check_indexes(pInd(l_ii));
    ok_ind = check_indexes(pInd(strcmp(u_pTyp{ii},pTyp)));
end

if ~ok_ind
    warning('Inconsistent numbering of participants'' label.');
end

end

%%

function ok_ind = check_indexes(pInd)
% Function to check that all the pInd are different or NaN.

if all(isnan(pInd))
    % if all NaN no check needed.
    ok_ind = true;
else
    if numel(pInd)~=numel(unique(pInd)) || any(isnan(pInd))
        % if not all unique or some NaN -> not OK
        ok_ind = false;
    else
        ok_ind = true;
    end
end

end

%%

function fn = get_bold_data(func_str,ffilt_func)
% Function to extract the list of functional data based on ffilt_func
% 
% INPUT
% - func_str    : structure array with the bold data
% - ffilt_func  : filtering structure for the functional data
%   .TaskLab    : task label filter
%
% NOTE:
% this needs to be extended in order to account for other filtering option
% such as 'AcqLab', 'RunLab' and 'RecLab', the acquisition, reconstruction 
% and run flags.
% -> add other fields in ffilt_func

lTasks = char(func_str.task);
lAcq = char(func_str.acq);
lRec = char(func_str.rec);
lRun = char(func_str.run);

end

% Add function to split runs according to numbers!