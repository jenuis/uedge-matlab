function path = abspath(dir_struct)
    if length(dir_struct) > 1
        path = {};
        for i=1:length(dir_struct)
            path{i} = abspath(dir_struct(i));
        end
        return
    end
    
    assert(isfield(dir_struct, 'folder'), '"dir_struct" has no field "folder"!')
    assert(isfield(dir_struct, 'name'), '"dir_struct" has no field "name"!')
    path = fullfile(dir_struct.folder, dir_struct.name);
end