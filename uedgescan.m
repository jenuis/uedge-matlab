% Author: Xiang LIU@ASIPP
% E-mail: xliu@ipp.ac.cn
% Created: 2023-12-16
% Version: V 0.1.17
classdef uedgescan < handle   
    properties(Constant)
        scan_field_value = 'value';
        scan_field_uedgevar = 'uedgevar';
        scan_field_names = 'names';
        scan_field_amps = 'amps';
        job_file_prefix = 'job';
    end
    
    properties
        work_dir
        scan = struct()
        script_input
        script_image
        file_init
        jobinfo % a collection of all job files, used as cache for massive run, where job_get_files consumes significant time
    end
    
    methods(Static)
        function input_diff_list = scan_traverse(scan, value_field_name)
            %% check arguments
            if nargin < 2
                value_field_name = uedgescan.scan_field_value;
            end
            %% return if scan is alreay OK
            fnames = fieldnames(scan);
            fnames_len = length(fnames);
            values_len = zeros(fnames_len);
            for i=1:fnames_len
                fn = fnames{i};
                values_len(i) = length(scan.(fn).(value_field_name));
            end
            
            multiple_values_inds = values_len > 1;
            fnames_traverse = fnames(multiple_values_inds);
            if isempty(fnames_traverse)
                input_diff_list = {scan};
                return
            end
            %% recursive call           
            input_diff_list = {};
            fname_traverse = fnames_traverse{1};
            values_traverse = scan.(fname_traverse).(value_field_name);
            for i=1:length(values_traverse)
                scan_tmp = scan;
                scan_tmp.(fname_traverse).(value_field_name) = values_traverse(i);
                input_diff_list = [input_diff_list uedgescan.scan_traverse(scan_tmp, value_field_name)];
            end
        end
        
        function input_diff_list_disp(input_diff_list, value_field_name)
            if nargin < 2
                value_field_name = uedgescan.scan_field_value;
            end
            
            for i=1:length(input_diff_list)
                disp(uedgerun.input_diff_gen_str(input_diff_list{i}, 'ValueFieldName', value_field_name))
            end
        end
        
        function job_file_init = file_init_find_nearest(work_dir, scan, input_diff)
            %% job_file_init = uedgescan.file_init_find_nearest(work_dir, scan, input_diff)
            % Match only one scan parameter for input_diff in sequence
            % and find the nearest savedt file.
            
            %% init output
            job_file_init = '';
            fname_value = uedgescan.scan_field_value;
            %% find the nearest savedt file
            % scan all space is significantly faster than loading
            % job file. it is also robust than parsing file name by
            % regular expression. still time consuming.
            scan_names = fieldnames(input_diff);
            for i=1:length(scan_names)
                %% no check dummpy scan parameter
                scan_name = scan_names{i};
                scan_values = scan.(scan_name).(fname_value);
                if length(scan_values) < 2
                    continue
                end
                %% get savedt files
                input_diff_tmp = input_diff;
                input_diff_tmp.(scan_name).(fname_value) = '*';
                savedt_files = dir(fullfile(work_dir, ...
                        uedgerun.generate_file_name(input_diff_tmp)));
                if isempty(savedt_files)
                    continue
                end
                if length(savedt_files) == 1
                    job_file_init = abspath(savedt_files);
                    return
                end
                %% find nearest value
                [~, inds] = sort(abs(scan_values - input_diff.(scan_name).(fname_value)));
                scan_values = scan_values(inds);
                init_files = cell(length(scan_values), 1);
                for j=1:length(scan_values)
                    input_diff_tmp.(scan_name).(fname_value) = scan_values(j);
                    init_files{j} = uedgerun.generate_file_name(input_diff_tmp);
                end
                [~, inds] = intersect(init_files, {savedt_files(:).name});
                job_file_init = fullfile(work_dir, init_files{min(inds)});
                return
            end
        end
        
        function file_init_disp(work_dir)
            savedt_files = abspath(...
                dir( fullfile(work_dir, [uedgerun.file_save_prefix '*.mat']) )...
                );
            for i=1:length(savedt_files)
                file_save = savedt_files{i};
                job = matread(file_save, 'job');
                disp('---------------------------------------------------')
                if isfield(job, 'file_init_used')
                    file_init = job.file_init_used;
                else
                    file_init = job.file_init;
                end
                [~, file_init_name] = fileparts(file_init);
                [~, file_save_name] = fileparts(file_save);
                disp(['Init: ' file_init_name])
                disp(['Save: ' file_save_name])
            end
        end
        
        function job = run_basic(job, varargin)
            %% check arguments
            Args.SkipExist = true;
            Args.ID = [];
            Args.SaveJob = true;
            Args.FailDirName = 'fail';
            Args.Debug = false;
            Args = parseArgs(varargin, Args, {'SaveJob', 'SkipExist', 'Debug'});
            
            assert(length(job) == 1, 'This function only run single job!')
            %% process arguments
            script_image = '';
            if isfield(job, 'script_image')
                script_image = job.script_image;
            end
            
            if isempty(Args.ID) && isfield(job, 'id')
                Args.ID = job.id;
            end
            
            work_dir = '.';
            if isfield(job, 'work_dir')
                work_dir = job.work_dir;
            end
            
            disp_prefix = [uedgerun.print_prefix ' '];
            if ~isempty(Args.ID)
                disp_prefix = ['[' num2str(Args.ID) ']' disp_prefix];
            end
            %% check fail
            fail_dir = fullfile(work_dir, Args.FailDirName);
            if ~exist(fail_dir, 'dir')
                mkdir(fail_dir);
            end
            
            job_file_fail = uedgerun.generate_file_name(job.input_diff, 'suffix', 'fail', 'extension', '.mat');
            job_file_fail = fullfile(fail_dir, job_file_fail);
            if exist(job_file_fail, 'file')
                job.status = NaN;
                job.reason = 'skipped, failed';
                disp([disp_prefix 'Skip failed job: ' job_file_fail])
                return
            end
            %% check or find the init file
            job_file_init = '';
            if contains(lower(job.file_init), 'near')
                self = matread(['scan_' work_dir '.mat'], 'self');
                job_file_init = uedgescan.file_init_find_nearest(work_dir, self.scan, job.input_diff);
                if isempty(job_file_init)
                    job.status = NaN;
                    job.reason = 'pending, nearest init file not ready';
                    disp([disp_prefix 'Waiting for nearest init file, temporarily skipped ...'])
                    return
                end
                job.file_init_used = job_file_init;
            else
                job_file_init = uedgerun.check_existence({job.file_init, fullfile(work_dir, job.file_init)});
            end
            %% run job            
            if isempty(job_file_init) || ~exist(job_file_init, 'file')
                %% omit run
                job.status = NaN; % no init file
                job.reason = 'unavailable';
                disp([disp_prefix 'No init file: ' job_file_init])
            else
                %% set input_diff
                ur = uedgerun(job.script_input, job_file_init, 'ImageScript', script_image);

                fnames = fieldnames(job.input_diff);
                for i=1:length(fnames)
                    fn = fnames{i};
                    field = job.input_diff.(fn);
                    ur.input_diff_update(fn, field.value, field.uedgevar.names, field.uedgevar.amps);
                end
                %% check existence
                ur.file_save = fullfile(work_dir, ur.generate_file_name(ur.input_diff));
                if Args.SkipExist && exist(ur.file_save, 'file')
                    job.status = NaN;
                    job.reason = 'skipped, successful';
                    disp([disp_prefix 'Skip exist: ' ur.file_save])
                    return
                end
                %% run
                ur.script_run_gen('ID', Args.ID);
                tic
                if Args.Debug
                    job.status = 0;
                    job.reason = '';
                    fclose(fopen(ur.file_save, 'w')); % dummy savedt file
                    pause(1);
                else
                    [job.status, job.reason] = ur.run();
                end
                job.elapsed_time = toc;
                disp([disp_prefix 'Total time used: ' num2str(job.elapsed_time/60) ' minutes!'])
            end
            %% analyze job            
            reason = lower(job.reason);
            % killed by "kill -9"
            if contains(reason, 'kill')
                disp([disp_prefix 'Job KILLED at "' ur.input_diff_gen_str(ur.input_diff) '"'])
                return
            end
            % convergence status
            if isempty(reason)
                job_file = strrep(ur.file_save, ur.file_extension, '.mat');
                run_status = 'successful';
            else
                job_file = job_file_fail;
                run_status = 'failed';
            end
            
            if Args.SaveJob
                save(job_file, 'job')
            end
            disp([disp_prefix 'Converged status: ' upper(run_status)])
        end
        
        function status = run_job_file(job_file, varargin)
            %% check if file exist
            status = '';
            fpath = abspath(job_file);
            if ~exist(fpath, 'file')
                warning(["File not exist: " fpath])
                return
            end
            
            disp_prefix = [uedgerun.print_prefix ' '];
            jobs = matread(fpath, 'job');
            if ~isempty(jobs(1).id)
                disp_prefix = ['[' num2str(jobs(1).id) ']' disp_prefix];
            end
            work_dir = jobs(1).work_dir;
            %% start diary
            [~, job_file_name] = fileparts(job_file.name);
            diary_file = fullfile(job_file.folder, [job_file_name '_log.txt']);
            diary_file = uedgerun.get_increment_file_name(diary_file);
            diary(diary_file);
            diary on
            disp([disp_prefix 'Created diary file: "' diary_file '"'])
            %% check if the job is running
            if uedgescan.run_lock(fpath, 'check')
                disp([disp_prefix 'Skip running job: "' job_file.name '"!'])
                status = 'running';
                return
            end
            %% run jobs
            uedgescan.run_lock(fpath, 'lock');
            for i=1:length(jobs)
                disp([disp_prefix '-------------------------------------'])
                job_new = uedgescan.run_basic(jobs(i), varargin{:});
                reason = lower(job_new.reason);
                if isempty(reason)
                    continue
                end
                
                if contains(reason, 'kill')
                    break
                end
            end
            %% delete job file
            disp([disp_prefix 'Completed: "' fpath '"'])
            delete(fpath)
            uedgescan.run_lock(fpath, 'unlock');
            %% stop diary
            diary off
        end
        
        function status = run_lock(lock_file, mode)
            %% check arguments
            if nargin < 2
                mode = 'check';
            else
                mode = lower(mode);
            end
            
            assert(haselement({'lock', 'unlock', 'check'}, mode), '"mode" should be in {"lock", "unlock", "check"}')
            lock_file_new = [lock_file '.lock']; % avoiding accident overwrite
            %% check
            if strcmpi(mode, 'check')
                status = exist(lock_file_new, 'file');
                return 
            end
            %% write lock file
            if strcmpi(mode, 'lock')
                assert(~uedgescan.run_lock(lock_file, 'check'), 'already locked!')
                uedgerun.script_save(lock_file_new, '');
                return
            end
            %% delete lock file
            assert(strcmpi(mode, 'unlock'), 'should be the only left mode!')
            assert(uedgescan.run_lock(lock_file, 'check'), 'lock file should exist!')
            delete(lock_file_new)
        end
        
        function log_print_new(file_id, varargin)
            %% check arguments
            Args.PrintPrefix = '';
            Args = parseArgs(varargin, Args);
            %% get new content range
            current_position = ftell(file_id);
            fseek(file_id, 0, 'eof');
            end_position = ftell(file_id);
            %% print new content
            if end_position > current_position
                fseek(file_id, current_position, 'bof');
                contents = fread(file_id, end_position - current_position, '*char')';
                if ~isempty(Args.PrintPrefix)
                    contents = strsplit(contents, '\n');
                    contents = contents(contains(contents, Args.PrintPrefix));
                    contents = strjoin(contents, '\n');
                end
                fprintf(contents);
            end
        end
        
        function logfile_new = log_trim(logfile, varargin)
            %% check arguments
            Args.ExtraPrintPrefix = {};
            Args = parseArgs(varargin, Args);
            
            if ischar(Args.ExtraPrintPrefix)
                Args.ExtraPrintPrefix = {Args.ExtraPrintPrefix};
            end
            
            print_prefix_list = {uedgerun.print_prefix};
            print_prefix_list = [print_prefix_list Args.ExtraPrintPrefix];
            
            assert(exist(logfile, 'file'), 'File not exist!')
            %% check logfile_new
            [path, name, ext] = fileparts(logfile);
            logfile_new = fullfile(path, [name '_trim' ext]);
            if exist(logfile_new, 'file')
                warning('File has already trimmed! Current file will be overwritten!')
            end
            %% trim log 
            fh = fopen(logfile, 'r');
            diary(logfile_new);
            diary on
            while(~feof(fh))
                l = fgetl(fh);
                for i=1:length(print_prefix_list)
                    if contains(l, print_prefix_list{i})
                        disp(l)
                    end
                end
            end
            diary off
        end   
    end
    
    methods(Access=private)
        function job_file_init = jobinfo_get_init_profile(self, jobid)
            job_file_init = '';
            
            jobinfo_id = [self.jobinfo(:).id];
            inds = jobinfo_id == jobid;
            if sum(inds) == 0
                return
            end
            
            job_file_init = {self.jobinfo(inds).file_init};
            job_file_init = job_file_init{1};
        end
    end
    
    methods
        function self = uedgescan(work_dir, script_input, file_init, varargin)
            %% check arguments
            self.work_dir = work_dir;
            self.scan_load();
            
            if nargin == 1    
                return
            end
            
            Args.ImageScript = '';
            Args = parseArgs(varargin, Args);
            %% set properties
            self.work_dir = work_dir;
            if ~exist(self.work_dir, 'dir')
                mkdir(self.work_dir);
            end
            self.script_input = uedgerun.check_existence(script_input, 1);
            self.file_init = uedgerun.check_existence(file_init, 1);
            self.script_image = uedgerun.check_existence(Args.ImageScript);
            %% copy files for uedgestat and uedgedata to work
            copyfile(self.script_input, self.work_dir)
            copyfile('mesh.hdf5', self.work_dir)
        end
        
        function is_loaded = scan_load(self)
            scan_file = ['scan_' self.work_dir '.mat'];
            is_loaded = false;
            if ~exist(scan_file, 'file')
                disp(['"' scan_file '" not exist!'])
                return
            end
            file_self = matread(scan_file, 'self');
            fnames = fieldnames(file_self);
            for i=1:length(fnames)
                fn = fnames{i};
                if ~isequal(self.(fn), file_self.(fn))
                    self.(fn) = file_self.(fn);
                end
            end
            is_loaded = true;
        end
        
        function scan_save(self)
            file = ['scan_' self.work_dir '.mat'];
            save(file, 'self')
        end
        
        function scan_add(self, scan_name, scan_values, uedgevar_names, uedgevar_amps)
            %% scan_add(self, scan_name, scan_values, uedgevar_names, uedgevar_amps)
            
            %% set properties
            self.scan.(scan_name).(self.scan_field_value) = unique(sort(scan_values));
            self.scan.(scan_name).(self.scan_field_uedgevar).(self.scan_field_names) = uedgevar_names;
            self.scan.(scan_name).(self.scan_field_uedgevar).(self.scan_field_amps) = uedgevar_amps;
            %% save
            self.scan_save();
        end
        
        function input_diff = scan_find_seed(self)
            %% load the input script
            contents = uedgerun.input_modify(self.script_input);
            contents = cellfun(@uedgerun.input_clean_line, contents, 'UniformOutput', false);
            contents = contents(contains(contents, '='));
            
            matches = cellfun(@(line) regexp(line, '^(.*?)\s*=\s*(.*?)$', 'tokens', 'once'), contents, 'UniformOutput', false);
            keys    = cellfun(@(m) m{1}, matches, 'UniformOutput', false);
            values  = cellfun(@(m) m{2}, matches, 'UniformOutput', false);
            %% find the nearest input_diff
            input_diff = self.scan;
            scan_names = fieldnames(self.scan);
            for i=1:length(scan_names)
                % find the first uedgevar name
                uedgevar_name_tmp = self.scan.(scan_names{i}).(self.scan_field_uedgevar).(self.scan_field_names);
                if ischar(uedgevar_name_tmp)
                    uedgevar_name = uedgevar_name_tmp;
                else
                    uedgevar_name = uedgevar_name_tmp{1};
                end
                % find the key and val in the input
                inds = cellfun(@(key) strcmpi(key, uedgevar_name), keys);
                assert(sum(inds) > 0, ['Can not find "' uedgevar_name '" in "' self.script_input '"!'])
                
                key_input = keys(inds);   key_input = key_input{end};
                val_input = values(inds); val_input = eval(val_input{end});
                % find the nearest value in self.scan and set in output
                scan_name = scan_names{i};                
                input_diff_vals = input_diff.(scan_name).(self.scan_field_value);
                input_diff_amp  = input_diff.(scan_name).(self.scan_field_uedgevar).(self.scan_field_amps);
                ind = findvalue(input_diff_vals * input_diff_amp(1), val_input);
                disp(['      ' scan_name ' = ' num2str(val_input) ' -> ' num2str(input_diff_vals(ind)) ...
                    ' * ' num2str(input_diff_amp(1))])
                input_diff.(scan_name).(self.scan_field_value) = input_diff_vals(ind);
            end
        end
        
        function job = job_struct(self, input_diff, file_init, varargin)
            %% check arguments
            Args.ScriptImage = self.script_image;
            Args.ScriptInput = self.script_input;
            Args.WorkDir = self.work_dir;
            Args.JobID = [];
            Args = parseArgs(varargin, Args);
            %% check input_diff
            uedgerun.input_diff_check(input_diff);
            %% set outputs
            job = struct();               
            job.input_diff = input_diff;
            job.file_init = file_init;
            job.script_input = Args.ScriptInput;
            job.script_image = Args.ScriptImage;
            job.work_dir = Args.WorkDir;
            job.id = Args.JobID;
            %% set jobinfo as cache for massive run
            if length(self.jobinfo) > 1 && ... 
                    self.jobinfo(end).id == job.id
                return
            end
            
            self.jobinfo(end+1).id = job.id;
            self.jobinfo(end).file_init = job.file_init;
            self.jobinfo(end).input_diff = job.input_diff;
        end
        
        function status = job_status(self, input_diff, varargin)
            Args.FailDirName = 'fail';
            Args = parseArgs(varargin, Args);
            
            profile_job_save = uedgerun.generate_file_name(input_diff, 'extension', '.mat');
            profile_job_save = fullfile(self.work_dir, profile_job_save);
            if exist(profile_job_save, 'file')
                status = 'successful';
                return
            end
            
            profile_job_fail = uedgerun.generate_file_name(input_diff, 'extension', '.mat', 'suffix', 'fail');
            profile_job_fail = fullfile(self.work_dir, Args.FailDirName, profile_job_fail);
            if exist(profile_job_fail, 'file')
                status = 'failed';
                return
            end
            
            status = 'pending';
        end
        
        function job_file = job_name_gen(self, jobid)
            job_file = fullfile(self.work_dir, [self.job_file_prefix num2str(jobid) '.mat']);
        end
        
        function job_gen_sequential(self, seq_scan_name, seq_start_value, varargin)
            %% check arguments
            Args.InitialRunFileInit = self.file_init;
            Args.Scan = self.scan;
            Args = parseArgs(varargin, Args);
            
            disp_prefix = [uedgerun.print_prefix ' '];
            %% get seq info
            seq = Args.Scan.(seq_scan_name);
            seq_values = seq.(self.scan_field_value);
            seq_start_value = seq_values(findvalue(seq_values, seq_start_value));
            %% clean
            delete(fullfile(self.work_dir, [self.job_file_prefix '*.mat']))
            %% pre run for using nearest init file
            if contains(lower(Args.InitialRunFileInit), 'near')
                input_diff_initial = self.scan_find_seed();
                jobid = 0;
                job = self.job_struct(input_diff_initial, self.file_init, 'jobid', jobid);
                job_file = self.job_name_gen(jobid);
                disp([disp_prefix 'Generating: "' job_file '"'])
                save(job_file, 'job')
            end
            %% initial run
            % Given the sequential scan name and value, the remain scan
            % parameters are traversed to start using Args.InitialRunFileInit
            
            scan_start = Args.Scan;
            scan_start.(seq_scan_name).(self.scan_field_value) = seq_start_value;
            
            input_diff_list = self.scan_traverse(scan_start);
            for jobid=1:length(input_diff_list)
                job = self.job_struct(input_diff_list{jobid}, Args.InitialRunFileInit, 'jobid', jobid);
                job_file = self.job_name_gen(jobid);
                disp([disp_prefix 'Generating: "' job_file '"'])
                save(job_file, 'job')
            end
            %% sequential run
            scan_remain = rmfield(Args.Scan, seq_scan_name);
            input_diff_remain_list = self.scan_traverse(scan_remain);
            
            seq_values_smaller = flip(seq_values(seq_values < seq_start_value));
            seq_values_larger = seq_values(seq_values > seq_start_value);
            seq_values_cell = {seq_values_smaller, seq_values_larger};
            
            for i=1:length(input_diff_remain_list)
                input_diff = input_diff_remain_list{i};
                input_diff.(seq_scan_name) = seq; 
                for j=1:length(seq_values_cell)
                    seq_values_tmp = seq_values_cell{j};
                    if isempty(seq_values_tmp)
                        continue
                    end
                    
                    jobid = jobid + 1;
                    for k=1:length(seq_values_tmp)
                        if k==1
                            value_init = seq_start_value;
                        else
                            value_init = seq_values_tmp(k-1);
                        end
                        input_diff.(seq_scan_name).(self.scan_field_value) = value_init;
                        file_init_tmp = fullfile(self.work_dir, uedgerun.generate_file_name(input_diff));
                        
                        value = seq_values_tmp(k);
                        input_diff.(seq_scan_name).(self.scan_field_value) = value;
                        
                        job_tmp = self.job_struct(input_diff, file_init_tmp, 'jobid', jobid);
                        if k==1
                            job = job_tmp;
                        else
                            job(k) = job_tmp;
                        end
                    end
                    job_file = self.job_name_gen(jobid);
                    disp([disp_prefix 'Generating: "' job_file '"'])
                    save(job_file, 'job')
                end
            end
            %% save
            self.scan_save()
        end
                
        function job_gen_rerun(self)
            disp_prefix = [uedgerun.print_prefix ' '];
            savedt_file_list = self.profile_files_get('succeed', ...
                'fileextension', uedgerun.file_extension);
            savedt_filenames = {savedt_file_list(:).name};
            input_diff_list = self.scan_traverse(self.scan);
            for jobid=1:length(input_diff_list)
                input_diff = input_diff_list{jobid};
                scan_str = uedgerun.input_diff_gen_str(input_diff);
                if ~any(contains(savedt_filenames, scan_str))
                    continue
                end           
                                
                job = self.job_struct(input_diff, file_init_tmp, 'jobid', jobid);
                job_file = self.job_name_gen(jobid);
                disp([disp_prefix 'Generating: "' job_file '"'])
                save(job_file, 'job')
            end
        end
        
        function job_files = job_get_files(self, varargin)
            %% job_files = job_get_files(self, 'ExcludeRunning', false, 'ExcludeNoInitProfile', false)
            % Excluded cases are used in self.run to avoid constantly adding
            % tasks that are not ready.
            % TODO: still too time consuming for massive task. add argument
            % that max number of job is required
            
            %% check arguments
            Args.ExcludeRunning = false;
            Args.ExcludeNoInitProfile = false;
            Args.MassiveBoostJobNum = 100;
            
            Args = parseArgs(varargin, Args, {'ExcludeRunning', 'ExcludeNoInitProfile'});
            disp_prefix = '      '; % not really need to keep for trimming
            is_parallel = ~isempty(gcp('nocreate')); % TODO: cellfun parallel?
            tic
            %% get all jobs
            job_file_pattern = [self.job_file_prefix, '*.mat'];
            job_files = dir(fullfile(self.work_dir, job_file_pattern));
            if isempty(job_files)
                return
            end
            
            job_paths = [];
            if length(job_files) < Args.MassiveBoostJobNum
                job_paths = abspath(job_files, 0);
                inds = cellfun(@(path) exist(path, 'file') ==2 && ...
                    strcmpi(who(matfile(path)), 'job'), ...
                    job_paths); % consuming time a lot when job number is large, which blocks dispathing
                job_files(~inds) = [];
                job_paths(~inds) = [];
            else
                warning([disp_prefix 'Skipped Job file check for massive run. Make sure no other ""' job_file_pattern '" file inside the work dir.'])
            end
            fprintf(1, [disp_prefix 'Enumeration of job files: '])
            toc
            %% remove running jobs
            if Args.ExcludeRunning
                if isempty(job_paths)
                    job_paths = abspath(job_files, 0);
                end
                
                lock_files = dir(fullfile(self.work_dir, [job_file_pattern '.lock']));
                if ~isempty(lock_files)
                    lock_job_paths = cellfun(@(path) path(1:end-5) , abspath(lock_files, 0), ...
                        'UniformOutput', false);
                    
                    [job_paths, inds] = setdiff(job_paths, lock_job_paths);
                    job_files = job_files(inds);
                end
                fprintf(1, [disp_prefix 'Enumeration of locked files: '])
                toc
            end
            %% remove jobs without init profile
            if Args.ExcludeNoInitProfile
                assert(length(self.jobinfo) >= length(job_files), '"self.jobinfo" is corrupted!')
                %% get jobid by parsing file name                
                jobid_list_cell = cellfun(...
                    @(job_name) str2double(regexp(job_name, [uedgescan.job_file_prefix '(\d+)\.mat'], 'tokens', 'once')), ...
                    {job_files(:).name}, 'UniformOutput', false);
                jobid_list = [jobid_list_cell{:}];
                if length(jobid_list) ~= length(jobid_list_cell)
                    error('Currently does not considering bad job file name, which should not be.')
                end
                %% get init files for job_files                
%                 jobinfo_init_files = cellfun(@self.jobinfo_get_init_profile, ...
%                 jobid_list_cell, 'UniformOutput', false); % very slow for massive run
                jobinfo_ids = [self.jobinfo(:).id];
                [~, inds_map] = ismember(jobid_list, jobinfo_ids); % map self.jobinfo fields to jobid_list
                init_file_list  = {self.jobinfo(inds_map).file_init};
                input_diff_list = {self.jobinfo(inds_map).input_diff};
                %% check if init file exist
                % bad job file
                inds_full = 1:length(job_files);
                inds_near = contains(init_file_list, 'near');
                inds_save = contains(init_file_list, uedgerun.file_save_prefix);
                inds_bad = inds_full(~inds_near & ~inds_save);
                % not exist
                inds_full_save = inds_full(inds_save);
                inds_tmp = cellfun(@(fname) exist(fname, 'file') ~= 2, ...
                    init_file_list(inds_save));
                inds_bad = [inds_bad inds_full_save(inds_tmp)];
                fprintf(1, [disp_prefix 'Enumeration of init files ("savedt*hdf5"): '])
                toc
                % no nearest init file
                inds_full_near = inds_full(inds_near);
                inds_tmp = cellfun(@(inpdiff) ...
                    isempty(self.file_init_find_nearest(self.work_dir, self.scan, inpdiff)), ...
                    input_diff_list(inds_near));
                inds_bad = [inds_bad inds_full_near(inds_tmp)];
                % remove bad
                job_files(inds_bad) = [];
                fprintf(1, [disp_prefix 'Enumeration of init files ("nearest"): '])
                toc
            end
        end
        
        function job_update(self)
            disp_prefix = [uedgerun.print_prefix ' '];
            job_files = self.job_get_files();
            for i=1:length(job_files)
                job_file = abspath(job_files(i));
                jobs = matread(job_file, 'job');
                %% check jobs
                index_pending = [];
                index_done = [];
                for j=1:length(jobs)
                    job = jobs(j);
                    status = self.job_status(job.input_diff);
                    if contains(lower(status), 'pending')
                        index_pending(end+1) = j;
                    end
                    if contains(lower(status), 'successful')
                        index_done(end+1) = j;
                    end
                end
                %% delete done
                if isempty(index_pending)
                    disp([disp_prefix 'Deleting completed: "' job_file '"'])
                    delete(job_file)
                    continue
                end                
            end
        end
        
        function file_profiles = profile_files_get(self, status, varargin)
            Args.FileExtension = '.mat';
            Args = parseArgs(varargin, Args);
            
            succeded_suffix = ['succe*' Args.FileExtension];
            failed_suffix = ['fail*' Args.FileExtension];
            
            ur = uedgerun();
            switch lower(status)
                case {'succeed', 'succeeded', 'successful'}
                    file_profiles = dir(fullfile(self.work_dir, [ur.file_save_prefix '*' succeded_suffix]));
                case {'fail', 'failed'}
                    file_profiles = dir(fullfile(self.work_dir, 'fail', [ur.file_save_prefix '*' failed_suffix]));
                otherwise
                    error(['Unknown status: "' status '"'])
            end
        end
           
        function diary_file = run(self, varargin)
            %% check arguments
            Args.SkipExist = true;
            Args.ID = [];
            Args.SaveJob = true;
            Args.FailDirName = 'fail';
            Args.Debug = false;
            Args.UseParallel = true;
            Args.NumWorkers = maxNumCompThreads;
            Args.LogFileName = 'uedgescan_log.txt';
            Args.PauseTime = 0.5;
            Args = parseArgs(varargin, Args, {'SaveJob', 'SkipExist', 'Debug', 'UseParallel'});
            
            use_parallel = Args.UseParallel;
            num_workers = Args.NumWorkers;
            log_file_name = Args.LogFileName;
            pause_time = Args.PauseTime;
            
            Args = rmfield(Args, 'UseParallel');
            Args = rmfield(Args, 'NumWorkers');
            Args = rmfield(Args, 'LogFileName');
            Args = rmfield(Args, 'PauseTime');
            varargin = struct2vararg(Args);
            disp_prefix = [uedgerun.print_prefix ' '];
            %% open diary
            start_time = now();
            diary_file = uedgerun.get_increment_file_name(log_file_name);
            diary(diary_file)
            diary on
            disp([disp_prefix 'Start time: ' datestr(start_time)])
            %% normal run
            if ~use_parallel
                while(1)
                    job_files = self.job_get_files('ExcludeRunning', 'ExcludeNoInitProfile');
                    if isempty(job_files)
                        disp([disp_prefix 'Exit with no jobs!'])
                        break
                    end

                    for i=1:length(job_files)
                        uedgescan.run_job_file(job_files(i), varargin{:})
                    end
                end
                return
            end
            %% parallel run
            %TODO: use onCleanup to detect use quit
            % uedge instance not closed by cancel(tasks)
            
            %% create parpool
            if isempty(gcp('nocreate'))
                disp([disp_prefix 'Creating parallel pool ...'])
                try
                    pool = parpool(num_workers);
                catch
                    pool = parpool();
                end
            else
                disp([disp_prefix 'Use existing parallel pool!'])
                pool = gcp;
            end
            %% main dispath
            tasks = [];
            logfile_handles = [];
            state = 'init'; % 'init', 'running', 'logging', 'finished'
            while(1)
                %% check if there are new jobs
                diary off
                job_files = self.job_get_files('ExcludeRunning', 'ExcludeNoInitProfile');
                diary on
                task_no = length(job_files);
                
                if task_no == 0 && (strcmpi(state, 'init') || strcmpi(state, 'finished'))
                    disp([disp_prefix 'Exit with no jobs!'])
                    break
                end
                %% create tasks     
                if task_no > 0                    
                    flag_create = false;
                    if strcmpi(state, 'init')
                        index_task = 0;
                        tasks = parallel.FevalFuture.empty(0, task_no);
                        flag_create = true;
                        disp([disp_prefix 'Creating parallel tasks (NO. is ' num2str(task_no) ') ...'])
                        state = 'running';
                    elseif sum( ...
                            strcmpi({tasks.State}, 'running') | ...
                            strcmpi({tasks.State}, 'pending') | ...
                            strcmpi({tasks.State}, 'queued') ...
                            ) < num_workers
                        index_task = length(tasks);
                        flag_create = true;
                        disp([disp_prefix 'Adding parallel tasks (NO. is ' num2str(task_no) ') ...'])
                    end
                    
                    if flag_create
                        for i = 1:task_no
                            % create task
                            job_file = job_files(i);
                            tasks(i+index_task) = parfeval(pool, ...
                                @self.run_job_file, ... % function name
                                1, ... % num of output parameters
                                job_file, varargin{:}); % fun arguments
                            % init log file handle
                            logfile_handles(i+index_task) = -1;
                        end
                    end
                end
                %% collect diary
                if strcmpi(state, 'running')
                    for i=1:length(logfile_handles)
                        % normal update log content
                        file_id = logfile_handles(i);
                        if file_id > 2
                            self.log_print_new(file_id, 'PrintPrefix', uedgerun.print_prefix);
                            continue
                        end
                        % check if task is running
                        if ~strcmpi(tasks(i).State, 'running')
                            continue
                        end
                        % open log file
                        job_file = tasks(i).InputArguments{1};
                        [~, job_file_name] = fileparts(job_file.name);
                        
                        pattern = fullfile(job_file.folder, [job_file_name '_log*.txt']);
                        job_diary = uedgerun.get_latest_file(pattern, start_time);
                        if isempty(job_diary) || ~exist(job_diary, 'file')
                            disp([disp_prefix 'Diary file currently not exist for "' job_file_name '"!'])
                            continue
                        end
                        
                        fh = fopen(job_diary, 'r');
                        if fh < 0
                            disp([disp_prefix 'Diary failed to open for "' job_file_name '"!'])
                            continue
                        end
                        
                        disp([disp_prefix 'Diary opened for "' job_file_name '"!'])
                        logfile_handles(i) = fh;
                    end
                    
                    if all(strcmpi({tasks.State}, 'finished'))
                        state = 'logging';
                    end
                end
                %% final logging
                if strcmpi(state, 'logging')
                    disp([disp_prefix 'Logging the rest ...'])
                    for i=1:length(logfile_handles)
                        file_id = logfile_handles(i);
                        if file_id > 2
                            self.log_print_new(file_id, 'PrintPrefix', uedgerun.print_prefix);
                        end
                    end
                    state = 'finished';
                end
                %% pause
                pause(pause_time)
                diary off
                disp(['Waiting ' num2str(pause_time) ' second ...'])
                diary on
            end
            %% clean
            for i=1:length(logfile_handles)
                file_id = logfile_handles(i);
                if file_id > 2
                    fclose(file_id);
                end
            end
            delete(pool)
            disp([disp_prefix 'End time: ' datestr(now)])
            diary off
        end
        
        function inspect(self, varargin)
            %% check arguments
            Args.Status = {'successful', 'failed', 'pending'};
            Args.Verbose = false;
            Args = parseArgs(varargin, Args, {'Verbose'});
            
            status_filter = Args.Status;
            %% disp job status
            job_files = self.job_get_files();
            for i=1:length(job_files)
                f = job_files(i);
                job = matread(f, 'job');
                disp(['***** ' f.name ' *****'])
                for j=1:length(job)
                    str_save = uedgerun.input_diff_gen_str(job(j).input_diff);
                    [~, profile_init] = fileparts(job(j).file_init);
                    status = self.job_status(job(j).input_diff);
                    if ~isempty(status_filter) && any(contains(lower(status_filter), lower(status(1:4))))
                        if Args.Verbose
                            disp(['"' str_save '" <-- "' profile_init '", "' upper(status) '"'])
                            continue
                        end
                        disp(['"' str_save '": "' upper(status) '"'])
                    end                    
                end
            end
        end
        
        function remain(self, varargin)
            job_list_all = self.scan_traverse(self.scan);
            len_all = length(job_list_all);
            disp(['All cases: ' num2str(len_all)])
            job_list_succeeded = self.profile_files_get('succeeded', ...
                'fileextension', uedgerun.file_extension);
            len_succeeded = length(job_list_succeeded);
            disp(['Succeded cases: ' num2str(len_succeeded)])
            job_list_failed = self.profile_files_get('failed', 'fileextension', '.mat');
            len_failed = length(job_list_failed);
            disp(['Failed cases: ' num2str(len_failed)])
            disp(['Remain cases: ' num2str(len_all - len_succeeded - len_failed)])
        end
             
        function clean(self)
            %% clean run and input files
            uedgerun.clean()
            %% clean not deleted job files
            self.job_update()
            %% clean lock files
            answer = input('Sure to clean job lock files? Please ensure there is no running task! [yes/no]', 's');
            if ~contains(lower(answer), 'y')
                return
            end
            
            lock_files = dir(fullfile(self.work_dir, 'job*.lock'));
            for i=1:length(lock_files)
                lock_file = abspath(lock_files(i));
                disp(['Remove: "' lock_file '"'])
                delete(lock_file)
            end
        end
        
        function time = time_stat(self)
            ur = uedgerun();
            file_list = dir(fullfile(self.work_dir, [ur.file_save_prefix '*.mat']));
            time = zeros(1, length(file_list));
            for i=1:length(file_list)
                f = file_list(i);
                job = matread(f, 'job');
                time(i) = job.elapsed_time;
            end
            uedgedata.figure;
            histogram(time/60);
            xlabel('Elapsed Time [minutes]')
            ylabel('Counts')
            uedgedata.figure_decoration;
        end
    end
end