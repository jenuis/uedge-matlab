function varargout = which2(mfile)
    %% check arguments
    [~, mfile] = fileparts(mfile);
    %% get $MATLABPATH
    matlabpath_env = getenv('MATLABPATH');
    paths = strsplit(matlabpath_env, ':');
    %% locate m file
    varargout = {};
    for i = 1:length(paths)
        if isempty(paths{i})
            continue
        end
        full_path = [paths{i} filesep mfile '.m'];
        if exist(full_path, 'file') == 2 || exist(full_path, 'class') == 8
            if nargout == 0
                disp(full_path)
            else
                varargout = {full_path};
            end
            return;
        end
    end
    
end