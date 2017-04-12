function [BIDSall,lDir] = crc_BIDS_checkDS(rootDir)
% Checking a series of BIDS-ds.
% Use those provided as an example from the BIDS account:
% https://github.com/INCF/BIDS-examples
%_______________________________________________________________________
% Copyright (C) 2017 Cyclotron Research Centre

% Written by C. Phillips.
% Cyclotron Research Centre, University of Liege, Belgium

if nargin==0
    rootDir = pwd;
end

% Get list of subdirs
Dcontent = dir(rootDir);
lDir = char(Dcontent([Dcontent.isdir]).name);
% Matlab always returns current (.) and upper (..) directories à la Unix
% -> remove them
if size(lDir,1)>2
    lDir = lDir(3:end,:);
else
    lDir = [];
end
nDir = size(lDir,1);


BIDSall = cell(nDir,1);
for ii = 1:nDir
    try
        BIDSall{ii} = spm_BIDS(fullfile(rootDir,deblank(lDir(ii,:))));
    catch
        BIDSall{ii} = 'problem';
    end
end

end
