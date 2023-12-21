function num = dimnum(value, warning_off)
    %% check argument
    if nargin < 2
        warning_off = false;
    end
    assert(isnumeric(value), 'Input value is not a numerical type!')
    %% get size of value
    value_size = size(value);
    num = length(value_size);
    num_real = sum(value_size > 1);
    %% number of dimension >= 3
    if num >= 3
        if ~warning_off && num > num_real
            warning(['"Value" can be shrinked to "' num2str(num_real) 'd" !'])
        end
        return
    end
    %% 2D
    if num == num_real
        return
    end
    %% 0D or 1D
    num = num_real;
end