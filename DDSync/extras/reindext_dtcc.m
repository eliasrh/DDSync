% reindex_dtcc.m
% -------------------------------------------------------------------------
% Reads 'catalog.txt', reindexes the event IDs 1,2,3,… in order of appearance,
% writes both 'catalog_new.txt' and 'dt_new.cc' with those new indices.

clearvars; clc;

%% 1) Read catalog.txt and build ID→newIndex map
fidCat = fopen('catalog_long.txt','r');
if fidCat<0
    error('Couldn''t open catalog.txt');
end

oldIDs = [];
catalogLines = {};
while ~feof(fidCat)
    tline = fgetl(fidCat);
    catalogLines{end+1} = tline; %#ok<AGROW>
    if isempty(tline), continue; end
    parts = strsplit(strtrim(tline));
    evtID = str2double(parts{end});     % last column
    oldIDs(end+1) = evtID;              %#ok<AGROW>
end
fclose(fidCat);

% Keep only the first occurrence of each ID
[uniqueIDs, ia] = unique(oldIDs, 'stable');

% Build a map: original ID → new index (1-based)
idxMap = containers.Map(uniqueIDs, 1:numel(uniqueIDs));

%% 2) Write catalog_new.txt with reindexed IDs
fidCatIn  = fopen('catalog_long.txt','r');
fidCatOut = fopen('catalog.txt','w');
if fidCatIn<0 || fidCatOut<0
    error('Couldn''t open catalog.txt or create catalog_new.txt');
end
while ~feof(fidCatIn)
    tline = fgetl(fidCatIn);
    if ~ischar(tline), break; end
    % Match everything up to last ID, then the ID
    tokens = regexp(tline, '^(.*\S)\s+(\d+)\s*$', 'tokens','once');
    if isempty(tokens)
        fprintf(fidCatOut, '%s\n', tline);
    else
        body = tokens{1}; oldID = str2double(tokens{2});
        newID = idxMap(oldID);
        fprintf(fidCatOut, '%s %d\n', body, newID);
    end
end
fclose(fidCatIn);
fclose(fidCatOut);

%% 3) Process dt.cc and write dt_new.cc
fidIn  = fopen('dt_long.cc','r');
fidOut = fopen('dt.cc','w');
if fidIn<0 || fidOut<0
    error('Couldn''t open dt.cc or dt_new.cc');
end

while true
    tline = fgetl(fidIn);
    if ~ischar(tline), break; end    % EOF
    
    % Detect lines with a "#" event‐pair
    tok = regexp(tline, '^(?<pre>\s*)#\s*(\d+)\s+(\d+)(?<post>.*)$', 'names');
    if ~isempty(tok)
        % parse original IDs
        nums = sscanf(tline, '# %ld %ld', 2).';

        old1 = nums(1); old2 = nums(2);
        new1 = idxMap(old1);
        new2 = idxMap(old2);
        newLine = sprintf('%s#   %d   %d%s', tok.pre, new1, new2, tok.post);
        fprintf(fidOut, '%s\n', newLine);
    else
        fprintf(fidOut, '%s\n', tline);
    end
end

fclose(fidIn);
fclose(fidOut);

fprintf('Done. Wrote catalog_new.txt and dt_new.cc with %d unique events.\n', numel(uniqueIDs));
