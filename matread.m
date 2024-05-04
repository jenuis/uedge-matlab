function [attr_value, mf] = matread(fpath_or_dir, attr_name, raise_error)
    %% [attr_value, mf] = matread(fpath_or_dir, attr_name, raise_error)
    
    %% check arguments
    if nargin < 3
        raise_error = false;
    end
    
    if isstruct(fpath_or_dir)
        fpath_or_dir = abspath(fpath_or_dir);
    end
    
    if ~exist(fpath_or_dir, 'file')
        warning_or_error(['Mat-file not exist: "' fpath_or_dir '"!'], raise_error)
        return
    end
    
    attr_value = [];
    %% read attribute value
    mf = matfile(fpath_or_dir);
    if ~strcmpi(who(mf), attr_name)
        warning_or_error(['No attribute "' attr_name '"!'], raise_error);
        return
    end
    
    attr_value = mf.(attr_name);
end

function warning_or_error(msg, raise_error)
    if raise_error
        error(msg)
    else
        warning(msg)
    end
end