% Move 7th column to the end, preserving "plain" numeric formatting (no e-notation)
infile  = 'catalog_long.txt';
outfile = 'catalog_long_reordered.txt';

fin  = fopen(infile,'r');
assert(fin>0, 'Could not open %s', infile);
fout = fopen(outfile,'w');
assert(fout>0, 'Could not open %s', outfile);

% Read line-by-line so we can preserve text-style formatting in output
while true
    line = fgetl(fin);
    if ~ischar(line), break; end
    if isempty(strtrim(line))
        fprintf(fout, '\n');
        continue
    end

    % Parse numeric fields (space-delimited)
    vals = sscanf(line, '%f').';
    if numel(vals) < 7
        % If something odd happens, pass line through unchanged
        fprintf(fout, '%s\n', line);
        continue
    end

    % Move column 7 to the end (keep all columns; count unchanged)
    vals = [vals(1:6), vals(8:end), vals(7)];

    % Print with integer formatting for cols 1-6, "plain" for the rest
    % (adjust decimals if you want fixed width; this avoids scientific notation)
    fprintf(fout, '%4d %02d %02d %02d %02d %02d', round(vals(1:6)));
    for k = 7:numel(vals)
        fprintf(fout, ' %s', num2str(vals(k), '%.15g')); % no e-format unless truly unavoidable
    end
    fprintf(fout, '\n');
end

fclose(fin);
fclose(fout);
