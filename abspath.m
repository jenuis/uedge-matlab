function paths = abspath(dir_struct, single_element_output_char)
    %% check arguments
    if nargin < 2
        single_element_output_char = true;
    end
    %% call array fun to gen paths
    paths = arrayfun(@(d) fullfile(d.folder, d.name), dir_struct, 'UniformOutput', false);
    %% check output
    if single_element_output_char && length(paths) == 1
        paths = paths{1};
    end
end