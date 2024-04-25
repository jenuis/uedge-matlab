% Author: Xiang LIU@ASIPP
% E-mail: xliu@ipp.ac.cn
% Created: 2023-12-15
% Version: V 0.1.5
classdef uedgerun < handle
    properties(Access=private)
        script_run % temp, script will be deleted after calling self.run
        script_input_diff % temp, script will be deleted after calling self.run
    end
    
    properties(Constant)
        file_save_prefix = 'savedt';
        file_image_prefix = 'images';
        file_mesh_prefix = 'mesh';
        file_extension = '.hdf5';
    end
    
    properties
        script_input
        script_image
        script_rundt
        script_runinit
        
        input_diff = struct();
        
        file_init
        file_save
    end
    
    methods(Static)
        function file = check_existence(files, raise_error)
            %% check arguments
            if nargin < 2
                raise_error = 0;
            end
            
            if ischar(files)
                files = {files};
            end
            %% check script existence
            file = '';
            flag_exist = 0;
            for i=1:length(files)
                tmp = files{i};
                if exist(tmp, 'file')
                    file = tmp;
                    flag_exist = 1;
                    break
                end
            end
            %% raise error check
            if raise_error && ~flag_exist
                error(['"' char(files) '" Not Found!'])
            end
        end
        
        function script_name = script_name_gen(varargin)
            Args.Prefix = 'uedgerun';
            Args = parseArgs(varargin, Args);
            
            uuid = char(java.util.UUID.randomUUID);
            script_name = [Args.Prefix '_' uuid '.py'];
        end
        
        function script_save(file_script, contents)
            fid = fopen(file_script, 'w');
            fprintf(fid, '%s', contents);
            fclose(fid);
        end
        
        function input_diff_check(input_diff, varargin)
            Args.ValueFieldName = 'value';
            Args.UedgeVarName = 'uedgevar';
            Args = parseArgs(varargin, Args);
            
            assert(isstruct(input_diff), '"input_diff" is not a struct!')
            fnames = fieldnames(input_diff);
            if isempty(fnames)
                return
            end
            for i=1:length(fnames)
                fn = fnames{i};
                field = input_diff.(fn);
                assert(isfield(field, Args.ValueFieldName), ['Value field not found: "' fn '"'])
                assert(isfield(field, Args.UedgeVarName), ['Uedge variable field not found: "' fn '"'])
            end
        end
        
        function str = input_diff_gen_str(input_diff, varargin)
            Args.Delimiter = '-';
            Args.ValueFieldName = 'value';
            Args = parseArgs(varargin, Args);
            
            vfn = Args.ValueFieldName;
            
            fnames = sort(fieldnames(input_diff));
            str = {};
            for i=1:length(fnames)
                fn = fnames{i};
                str{i} = [fn num2str(input_diff.(fn).(vfn))];
            end
            
            if ~isempty(str)
                str = strjoin(str, Args.Delimiter);
            end
        end
        
        function str = parse_file_name(file_name, varargin)
            Args.Delimiter = '_';
            Args = parseArgs(varargin, Args);
            
            [~, file_name] = fileparts(file_name);
            file_name_splits = strsplit(file_name, Args.Delimiter);
            assert(length(file_name_splits) == 3, 'Fail to parse!')
            str = file_name_splits{2};
        end
        
        function file_name = generate_file_name(input_diff, varargin)            
            Args.Prefix = 'savedt';
            Args.Suffix = 'successful';
            Args.Delimiter = '_';
            Args.Extension = '.hdf5';
            Args = parseArgs(varargin, Args);
            
            name_strs = {Args.Prefix};
            input_diff_str = uedgerun.input_diff_gen_str(input_diff);
            if ~isempty(input_diff_str)
                name_strs{end+1} = input_diff_str;
            end
            name_strs{end+1} = Args.Suffix;
            
            file_name = strjoin(name_strs, Args.Delimiter);
            file_name = [file_name Args.Extension];
        end
        
        function line = input_fix_line(line)
            %% fix print for py2.X
            if contains(line, 'print') && ~contains(line, '(') && ~contains(line, ')')
                line_rep = strrep(line, 'print', 'print ');
                line = regexprep(line_rep, 'print\s+(.*?)\s*$', 'print($1)');
                line = regexprep(line, 'print(.*?),\s*$', 'print($1, end="")');
                
                if ~strcmp(line, line_rep)
                    disp(['-' line_rep])
                    disp(['+' line])
                end
            end
            %% using hdf5_restore
            if contains(line, 'restore') && ~contains(line, 'hdf5_restore')
                disp(['- ' line])
                line = ['from uedge.hdf5 import *; ' strrep(line, 'restore', 'hdf5_restore')];
                disp(['+ ' line])
            end
            %% no fix
        end
        
        function input_fix(input_script)
            contents = {};
            fid = fopen(input_script);
            while(1)
                %% end line, exit
                line = fgetl(fid);
                if ~ischar(line)
                    break
                end
                %% fix line
                line = uedgerun.input_fix_line(line);
                %% store to contents
                contents{end+1} = line;
            end
            fclose(fid);
            %% fix oldseec
            if ~any(contains(contents, 'oldseec'))
                line = 'bbb.oldseec = 0';
                contents{end+1} = line;
                disp(['+ ' line])
            end
            %% join to text
            contents = strjoin(contents, '\n');
            %% save
            uedgerun.script_save(input_script, contents);
        end
        
        function contents = input_modify(input_script, varargin)
            %% check arguments
            input_script = uedgerun.check_existence(input_script, 1);
            vararg_len = length(varargin);
            if vararg_len == 1 && iscell(varargin)
                varargin = varargin{1};
                vararg_len = length(varargin);
            end
            assert(mod(vararg_len, 2)==0, 'The input arguments are not paired!')
            %% scan and modify
            contents = {};
            fid = fopen(input_script);
            while(1)
                %% end line, exit
                line = fgetl(fid);
                if ~ischar(line)
                    break
                end
                %% add to contents
                contents{end+1} = line;
                %% check if this line is an input
                if ~contains(line, '=')
                    continue
                end
                %% get input key
                arg_name_org = strsplit(line, '=');
                arg_name_org = arg_name_org{1};
                %% compare and modify input
                for i=1:(vararg_len/2)
                    arg_name = varargin{2*i-1};
                    arg_val = varargin{2*i};
                    if ~strcmpi(strip(arg_name), strip(arg_name_org))
                        continue
                    end
                    line_new = [arg_name ' = ' num2str(arg_val)];
                    contents{end} = line_new;
                    disp(['- ' line])
                    disp(['+ ' line_new])
                end
            end
            fclose(fid);
            %% join to text
            contents = strjoin(contents, '\n');
        end
        
        function clean(tar_dir)
            if nargin == 0
                tar_dir = '.';
            end
            
            answer = input('Sure to clean running cache? Please ensure there is no running task! [yes/no]', 's');
            if ~contains(lower(answer), 'y')
                return
            end
            
            file_patterns = {'uedgerun_*.py', 'uedgeinput_*.py'};
            for i=1:length(file_patterns)
                pattern = file_patterns{i};
                file_list = dir(fullfile(tar_dir, pattern));
                for j=1:length(file_list)
                    f = abspath(file_list(j));
                    disp(['Remove: "' f '"'])
                    delete(f)
                end
            end
        end
    end
    methods
        function self = uedgerun(input_script, file_init, varargin)
            %% uedgerun(input_script, file_init, 'ImageScript', '', 'RundtScript', 'rdcontdt.py')
            
            %% check arguments
            if nargin == 0
                return
            end
            
            Args.ImageScript = {'image_save_new.py', 'image_save.py', '../image_save_new.py', '../image_save.py'};
            Args.RundtScript = {'rdcontdt.py', '../rdcontdt.py'};
            Args.RuninitScript = {'rdinitdt.py', '../rdinitdt.py'};
            
            Args = parseArgs(varargin, Args);
            %% set properties
            self.script_input = self.check_existence(input_script, 1);
            self.file_init = self.check_existence(file_init, 1);
            self.script_image = self.check_existence(Args.ImageScript);
            self.script_rundt = self.check_existence(Args.RundtScript, 1);
            self.script_runinit = self.check_existence(Args.RuninitScript, 1);
        end
        
        function args = input_diff_update(self, scan_name, scan_value, uedgevar_names, uedgevar_amplifications)
            %% recursive call
            if nargin < 2
                fnames = fieldnames(self.input_diff);
                args = {};
                for i=1:length(fnames)
                    scan_name = fnames{i};
                    args = [args self.input_diff_update(scan_name)];
                end
                return
            end
            if nargin < 4
                assert(isfield(self.input_diff, scan_name), 'Input "scan_name" not exist in "input_diff"!')
                scan = self.input_diff.(scan_name);
                if nargin > 2
                    scan.value = scan_value;
                end
                args = self.input_diff_update(scan_name, scan.value, scan.uedgevar.names, scan.uedgevar.amps);
                return
            end
            %% check arguments
            assert(nargin > 4, 'Input argments not enough!')
            if ischar(uedgevar_names)
                uedgevar_names = {uedgevar_names};
            end
            assert(iscellstr(uedgevar_names), '"uedgevar_names" should be cellstr or char!')
            assert(isnumeric(uedgevar_amplifications),'"uedgevar_amplifications" should be numeric!')
            assert(length(uedgevar_names) == length(uedgevar_amplifications), '"uedge_variables_*" length does not match!')
            scan = struct();
            scan.value = scan_value;
            scan.uedgevar.names = uedgevar_names;
            scan.uedgevar.amps = uedgevar_amplifications;
            %% update scan
            if isfield(self.input_diff, scan_name)
                scan_stored = self.input_diff.(scan_name);
                if ~isequal(scan_stored, scan)
                    if scan_stored.value ~= scan.value
                        disp(['"input_diff.' scan_name '.value" changed: ' num2str(scan_stored.value) ' -> ' num2str(scan.value)])
                    end

                    uv_names = scan_stored.uedgevar.names;
                    uv_amps = scan_stored.uedgevar.amps;
                    for i=1:length(uv_names)
                        [flag, index] = haselement(scan.uedgevar.names, uv_names{i});
                        if flag
                            if uv_amps(i) ~= scan.uedgevar.amps(index)
                                disp(['"input_diff.' scan_name '.uedgevar.' uv_names{i} '" amp changed: ' num2str(uv_amps(i)) ' -> ' num2str(scan.uedgevar.amps(index))])
                            end
                            continue
                        end
                        scan.uedgevar.names{end+1} = uv_names{i};
                        scan.uedgevar.amps(end+1) = uv_amps(i);
                    end
                end
            end
            self.input_diff.(scan_name) = scan;
            %% output args
            args = {};
            for i=1:length(scan.uedgevar.names)
                args{2*i-1} = scan.uedgevar.names{i};
                args{2*i} = scan.value * scan.uedgevar.amps(i);
            end
        end
        
        function file_input_new = script_input_diff_gen(self)
            file_input_new = self.script_name_gen('prefix', 'uedgeinput');
            args = self.input_diff_update();
            if isempty(args)
                warning('"input_diff" is empty, input file not being modfied!')
            end
            contents = self.input_modify(self.script_input, args{:});
            self.script_save(file_input_new, contents);
            self.script_input_diff = file_input_new;
        end
        
        function file_run = script_run_gen(self, varargin)
            %% check arguments
            Args.PrintPrefix = '>>>>> ';
            Args.ID = [];
            Args = parseArgs(varargin, Args);
            
            disp_prefix = Args.PrintPrefix;
            if ~isempty(Args.ID)
                disp_prefix = ['[' num2str(Args.ID) ']' disp_prefix];
            end
            %% gen file names
            file_run = self.script_name_gen();
            %% collect variables to be write into the script
            profile_init = self.check_existence(self.file_init, 1);
            profile_save = self.file_save;
            if isempty(profile_save)
                profile_save = self.generate_file_name(self.input_diff);
            end
            rd_in_script = self.script_input_diff_gen();
            image_script = self.script_image;
            rundt_script = self.check_existence(self.script_rundt, 1);
            runinit_script = self.check_existence(self.script_runinit, 1);
            %% gen run script content
            contents = {
                'import os', ...
                '', ...
                ['profile_init = "' profile_init '"'], ...
                ['profile_save = "' profile_save '"'], ...
                ['rd_in_script = "' rd_in_script '"'], ...
                ['image_script = "' image_script '"'], ...
                ['rundt_script = "' rundt_script '"'], ...
                ['runinit_script = "' runinit_script '"'], ...
                '', ...
                'exec(open(rd_in_script).read())', ...
                ['assert os.path.exists(profile_init), "' disp_prefix 'Init file not found: " + profile_init'], ...
                'if profile_init == profile_save:', ...
                ['    print("' disp_prefix 'Rerun file: " + profile_init)'],...
                'else:',...
                ['    print("' disp_prefix 'Loading init file: " + profile_init)'], ...
                'hdf5_restore(profile_init)', ...
                '', ...
                'bbb.exmain()', ...
                ['assert bbb.iterm == 1, "' disp_prefix 'bbb.iterm != 1"'], ...
                '', ...
                ['print("' disp_prefix 'Exec: " + rundt_script)'],...
                'exec(open(runinit_script).read())', ...
                'exec(open(rundt_script).read())', ...
                ['print("' disp_prefix 'Saving profile ...")'],...
                'hdf5_save(profile_save)', ...
                'casedir, _ = os.path.split(profile_save)', ...
                ''};
            if ~isempty(image_script)
                contents{end+1} = ['image_save = os.path.join(casedir, "' self.file_image_prefix self.file_extension '")'];
                contents{end+1} = ['if "' self.file_save_prefix '" in profile_save:'];
                contents{end+1} = ['    image_save = profile_save.replace("' self.file_save_prefix '", "' self.file_image_prefix '")'];
                contents{end+1} = ['print("' disp_prefix 'Exec: " + image_script)'];
                contents{end+1} = 'exec(open(image_script).read())';
            end
            contents = strjoin(contents, '\n');
            %% save script
            self.script_save(file_run, contents)
            self.script_run = file_run;
        end
        
        function status = run(self)
            assert(~isempty(self.script_run), 'Call "script_run_gen" to generate the run script!')
            file_run = self.check_existence(self.script_run, 1);
            status = system(['python3 ' file_run]);
            delete(file_run)
            delete(self.script_input_diff)
        end
    end
end