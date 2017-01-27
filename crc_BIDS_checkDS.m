function [BIDSall,lDir] = crc_BIDS_checkDS(rootDir)
% Checking a series of BIDS-ds.

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
