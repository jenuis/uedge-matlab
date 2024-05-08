function version_strs = getversion(varargin)
    %% Parse input arguments
    p = inputParser;
    p.addRequired('FileNames', @iscellstr);
    p.addParameter('VersionKeyword', 'version', @ischar);
    p.addParameter('PrintWarning', true, @islogical);
    p.addParameter('PrintVersion', true, @islogical);
    p.addParameter('PrintPrefix', '', @ischar);
    parse(p, varargin{:});

    filenames = p.Results.FileNames;
    version_keyword = p.Results.VersionKeyword;
    print_warning = p.Results.PrintWarning;
    print_version = p.Results.PrintVersion;
    print_prefix = p.Results.PrintPrefix;
    
    version_strs = {};
    %% Process each file
    for i = 1:length(filenames)
        filename = filenames{i};
        % Find the file on MATLAB's path
        try
            fullpath = which2(filename);
        catch
            fullpath = which(filename);
        end
        if isempty(fullpath)
            if print_warning
                warning('File %s not found on MATLAB''s path.', filename);
            end
            continue;
        end
        % Open the file and find the version line
        version_strs{i} = '';
        fid = fopen(fullpath, 'r');
        while(1)
            line = fgetl(fid);
            line = lower(strip(line));
            if line(1) == '%' && contains(line, version_keyword)
                version_strs{i} = strtrim(line(2:end));
                break
            end
            
        end
        fclose(fid);
        if isempty(version_strs{i})
            if print_warning
                warning('Can not find "%s" in "%s".', version_keyword, filename)
            end
            continue
        end
        % print version line
        if ~print_version
            continue
        end
        disp([print_prefix ' '...
            strrep(version_strs{i}, version_keyword, ['"' filename '" version'])])
    end
end