% Author: Xiang LIU@ASIPP
% E-mail: xliu@ipp.ac.cn
% Created: 2023-12-16
% Version: 0.1.18
% TODO: time out control
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
        file_init % TODO: use relative path, only file name.
        jobinfo % a collection of all job files, used as cache for massive 
                % run, where job_get_files consumes significant time
    end
    
    properties(Access=private)
        tasks;
        log_fid_list;
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
            
            if isstruct(input_diff_list)
                input_diff_list = {input_diff_list};
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
        
        function result = job_file_check(fpath, raise_error)
            if nargin < 2
                raise_error = true;
            end
            
            result = false;
            try
                result = exist(fpath, 'file') ==2 && ...
                    strcmpi(who(matfile(fpath)), 'job');
            catch
            end
            
            if raise_error && ~result
                error('Not a valid job file!')
            end
        end
        
        function pool = parpool_create(num_workers, print_prefix)
            if nargin < 2
                print_prefix = uedgerun.print_prefix;
            end
            
            pool = gcp('nocreate');
            if isempty(pool)
                disp([print_prefix ' Creating parallel pool ...'])
                try
                    pool = parpool(num_workers);
                catch
                    pool = parpool();
                end
            else
                disp([print_prefix ...
                    ' Using current parallel pool (num of workers: ' ...
                    num2str(pool.NumWorkers) ') ...'])
            end
        end
            
        function num = task_remain(tasks)
            num = sum( ...
                strcmpi({tasks.State}, 'running') | ...
                strcmpi({tasks.State}, 'pending') | ...
                strcmpi({tasks.State}, 'queued') ...
                );
        end
        
        function job = run_basic(job, varargin)
            %% job = uedgescan.run_basic(job, 'SkipExist', true, ...
            %      'ID', [], 'SaveJob', true, ...
            %      'FailDirName', 'fail', 'Debug', false)
            % Run a single job task, which does not check job lock.
            
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
            
            if isempty(Args.ID)
                disp_prefix = [uedgerun.print_prefix ' '];
            else
                disp_prefix = [uedgerun.print_prefix '[' num2str(Args.ID) '] '];
            end
            %% skip failed job
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
            %% get init file
            if contains(lower(job.file_init), 'near')
                % find the nearest init file
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
                % find the specified init file
                job_file_init = uedgerun.check_existence({job.file_init, fullfile(work_dir, job.file_init)});
            end
            %% run job
            if isempty(job_file_init) % specified init file not exist
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
        
        function status = run_job_file(fpath_or_dir, varargin)
            %% status = uedgescan.run_job_file(fpath_or_dir, ...
            %      'ID', [], 'SaveJob', true, ...
            %      'FailDirName', 'fail', 'Debug', false)
            % Run a job file, which may contain several job tasks. 
            % Locked jobs will be skipped.
            
            %% check if job file exists
            status = '';
            try
                fpath = abspath(fpath_or_dir);
            catch
                fpath = fpath_or_dir;
            end
            if ~exist(fpath, 'file')
                warning(["File not exist: " fpath])
                return
            end
            
            [job_file_folder, job_file_name] = fileparts(fpath);
            
            jobs = matread(fpath, 'job');
            if ~isempty(jobs(1).id)
                disp_prefix = [uedgerun.print_prefix '[' num2str(jobs(1).id) '] '];
            else
                disp_prefix = [uedgerun.print_prefix ' '];
            end
            %% start diary
            diary_file = fullfile(job_file_folder, [job_file_name '_log.txt']);
            diary_file = uedgerun.get_increment_file_name(diary_file);
            diary(diary_file);
            diary on
            disp([disp_prefix 'Created diary file: "' diary_file '"'])
            %% check if the job is running
            if uedgescan.run_lock(fpath, 'check')
                disp([disp_prefix 'Skip running job: "' job_file_name '.mat"!'])
                status = 'running';
                return
            end
            %% run jobs
            uedgescan.run_lock(fpath, 'lock');
            for i=1:length(jobs)
                if i > 1
                    disp([disp_prefix '-------------------------------------'])
                end
                job_new = uedgescan.run_basic(jobs(i), varargin{:});
                reason = lower(job_new.reason);
                if isempty(reason)
                    continue
                end
                
                if contains(reason, 'kill') % stop running the rest
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
        
        function version()
            getversion({mfilename}, ...
                'versionkeyword', 'version', ....
                'printwarning', false, ...
                'printversion', true, ...
                'printprefix', [uedgerun.print_prefix '[UEDGE-MATLAB]']);
        end
        
        function log_file = log_start(log_file_name, print_prefix)
            %% check arguments
            if nargin < 2
                print_prefix = uedgerun.print_prefix;
            end
            %% create log file and print information
            log_file = uedgerun.get_increment_file_name(log_file_name);
            diary(log_file)
            diary on
            disp('*******************************************************')
            uedgescan.version()
            uedgerun.version()
            disp('*******************************************************')
            disp([print_prefix ' Start time: ' datestr(now())])
            disp([print_prefix ' Log file is "' log_file '"'])
        end
        
        function log_end(print_prefix)
            if nargin < 1
                print_prefix = uedgerun.print_prefix;
            end
            disp([print_prefix ' End time: ' datestr(now)])
            diary off
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
                    if ~isempty(contents)
                        contents = strjoin(contents, '\n');
                        fprintf([contents '\n']);
                    end
                end
            end
        end
        
        function logfile_new = log_trim(logfile_or_folder, varargin)
            %% logfile_new = log_trim(logfile_or_folder, ...
            %    'ExtraPrintPrefix', {})
            
            %% check arguments
            Args.ExtraPrintPrefix = {};
            Args = parseArgs(varargin, Args);
            
            if ischar(Args.ExtraPrintPrefix)
                Args.ExtraPrintPrefix = {Args.ExtraPrintPrefix};
            end
            
            print_prefix_list = {uedgerun.print_prefix};
            print_prefix_list = [print_prefix_list Args.ExtraPrintPrefix];
            %% recursive call
            if exist(logfile_or_folder, 'dir')
                logfiles = dir(fullfile(logfile_or_folder, [uedgescan.job_file_prefix '*.txt']));
                for i=1:length(logfiles)
                    if contains(logfiles(i).name, 'trim')
                        continue
                    end
                    uedgescan.log_trim(abspath(logfiles(i)), varargin{:})
                end
                return
            end
            logfile = logfile_or_folder;
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
        
        function [scan_para_diff, scan_para_common] = savedt_compare(fpath1, fpath2)
            [path1, fname1] = fileparts(fpath1);
            [path2, fname2] = fileparts(fpath2);
            job_file1 = fullfile(path1, [fname1 '.mat']);
            job_file2 = fullfile(path2, [fname2 '.mat']);
            job1 = matread(job_file1, 'job', 1);
            job2 = matread(job_file2, 'job', 1);
            input_diff1 = job1.input_diff;
            input_diff2 = job2.input_diff;
            scan_names = fieldnames(input_diff1);
            scan_para_diff = struct([]);
            scan_para_common = struct([]);
            for i=1:length(scan_names)
                scan_name = scan_names{i};
                val1 = input_diff1.(scan_name).(uedgescan.scan_field_value);
                val2 = input_diff2.(scan_name).(uedgescan.scan_field_value);
                if val1 == val2
                    scan_para_common(1).(scan_name) = val1;
                else
                    scan_para_diff(1).(scan_name) = [val1 val2];
                end
            end
        end
        
        function savedt_files = savedt_find_missing(work_dir)
            savedt_files = dir(fullfile(work_dir, [uedgerun.file_save_prefix '*' uedgerun.file_extension]));
            inds = cellfun(@(path) ...
                exist( strrep(path, uedgerun.file_extension, '.mat'), 'file') ~= 2, ...
                abspath(savedt_files));
            savedt_files(~inds) = [];
        end
        
        function [fail_files, cache] = savedt_find_fail(fail_init, cache)
            %% get fail_files
            if nargin < 2
                cache.jobid = [];
                [fail_dir, ~, ext] = fileparts(fail_init);
                fail_files = dir(fullfile(fail_dir, [uedgerun.file_save_prefix '*' ext]));
                cache.all_fail_files = fail_files;
            else
                fail_files = cache.all_fail_files;
            end
            %% get target jobid
            job = matread(fail_init, 'job');
            jobid = job.id;
            %% slice fail_files
            if nargin < 2
                inds = false(size(fail_files));
                for i=1:length(fail_files)
                    job = matread(abspath(fail_files(i)), 'job');
                    jobid_tmp = job.id;
                    cache.jobid(i) = jobid_tmp;

                    if jobid_tmp ~= jobid
                        continue;
                    end
                    inds(i) = true;
                end
            else
                inds = cache.jobid == jobid;
            end

            fail_files = fail_files(inds);
        end
    end
    
    methods(Access=private)                
        function run_normal(self, varargin)
            disp_prefix = [uedgerun.print_prefix ' '];
            while(1)
                job_files = self.job_file_get_available('MaxRequired', 10);
                if isempty(job_files)
                    disp([disp_prefix 'Exit with no jobs!'])
                    break
                end

                for i=1:length(job_files)
                    uedgescan.run_job_file(job_files(i), varargin{:})
                end
            end
        end
        
        function job_files = task_trim_job_files(self, job_files, task_backup_ratio)
            %% check arguments
            if nargin < 3
                task_backup_ratio = 0.2;
            end
            
            if isempty(job_files)
                return
            end
            %% trim tasks according to parpool size
            num_workers = gcp().NumWorkers;
            max_num_new_tasks = floor( num_workers*(1+task_backup_ratio) );
            if isempty(self.tasks)
                if length(job_files) > max_num_new_tasks
                    job_files(max_num_new_tasks+1:end) = [];
                end
                return
            end
            %% check if there are still many backup tasks
            num_remain_tasks = self.task_remain(self.tasks);
            if num_remain_tasks > max_num_new_tasks
                job_files(:) = [];
                return
            end
            %% trim tasks accoring to num of free backup taks
            num_free_tasks = max_num_new_tasks - num_remain_tasks;
            if length(job_files) > num_free_tasks
                job_files(num_free_tasks+1:end) = [];
            end
        end
        
        function task_enqueue(self, job_files, varargin)
            %% check arguments
            if isempty(job_files)
                return
            end
            %% init self.tasks
            new_task_no = length(job_files);
            if isempty(self.tasks)
                self.tasks = parallel.FevalFuture.empty(0, new_task_no);
                ind_task_tail = 0;
            else
                ind_task_tail = length(self.tasks);
            end
            %% main
            disp([uedgerun.print_prefix ...
                ' Adding tasks (' num2str(new_task_no) ...
                ') into parallel pool (' ...
                num2str(self.task_remain(self.tasks)) ' left) ...'])
            for i = 1:new_task_no
                %% add task
                self.tasks(i+ind_task_tail) = parfeval(...
                    @self.run_job_file, ... % function name
                    1, ... % num of output parameters
                    job_files(i), varargin{:}); % fun arguments
                %% init log file handle
                self.log_fid_list(i+ind_task_tail) = -1;
            end
        end
        
        function task_remove_finished(self)
            %% check and print finished tasks
            inds = [];
            for i=1:length(self.tasks)
                if ~strcmpi(self.tasks(i).State, 'finished')
                    continue
                end
                job_file = self.tasks(i).InputArguments{1};
                [~, job_file_name] = fileparts(job_file.name);
                disp([uedgerun.print_prefix ' Task finished: "' job_file_name '.mat"! ' num2str(self.task_remain(self.tasks)) ' tasks left ...'])
                if ~isempty(self.tasks(i).Error)
                    warning([uedgerun.print_prefix 'This task returned error: ' self.tasks(i).Error.message])
                end
                fid = self.log_fid_list(i);
                if fid > 2
                    fclose(fid);
                end
                inds(end+1) = i;
            end
            %% remove finished tasks
            if ~isempty(inds)
                self.tasks(inds) = [];
                self.log_fid_list(inds) = [];
            end
        end
        
        function task_log_open(self, time_newer, varargin)
            Args.PrintNotExist = true;
            Args = parseArgs(varargin, Args, {'PrintNotExist'});
            disp_prefix = [uedgerun.print_prefix ' '];
            for i=1:length(self.log_fid_list)
                %% check if the log file has already been opened
                if self.log_fid_list(i) ~= -1
                    continue
                end
                %% check if the task is running
                if ~strcmpi(self.tasks(i).State, 'running')
                    continue
                end
                %% get diary file name
                job_file = self.tasks(i).InputArguments{1};
                [~, job_file_name] = fileparts(job_file.name);

                pattern = fullfile(job_file.folder, [job_file_name '_log*.txt']);
                job_diary = uedgerun.get_latest_file(pattern, time_newer);
                if isempty(job_diary) || ~exist(job_diary, 'file')
                    if Args.PrintNotExist
                        disp([disp_prefix 'Diary file currently not exist for "' job_file_name '.mat"!'])
                    end
                    continue
                end
                %% open diary
                fid = fopen(job_diary, 'r');
                if fid < 0
                    disp([disp_prefix 'Diary failed to open for "' job_file_name '.mat"!'])
                    continue
                end

                disp([disp_prefix 'Diary opened for "' job_file_name '.mat"!'])
                self.log_fid_list(i) = fid;
            end
        end
        
        function task_log_collect(self)
            for i=1:length(self.log_fid_list)
                fid = self.log_fid_list(i);
                if fid > 2
                    self.log_print_new(fid, 'PrintPrefix', uedgerun.print_prefix);
                end
            end
        end
        
        function task_log_close(self)
            for i=1:length(self.log_fid_list)
                fid = self.log_fid_list(i);
                if fid > 2
                    fclose(fid);
                end
            end
        end
        
        function run_parallel(self, varargin)
            %TODO: use onCleanup to detect use quit
            % uedge instance not closed by cancel(tasks)
            
            %% check arguments
            Args.SkipExist = true;
            Args.SaveJob = true;
            Args.FailDirName = 'fail';
            Args.Debug = false;
            Args.NumWorkers = maxNumCompThreads;
            Args.TaskBackupRatio = 0.2;
            Args = parseArgs(varargin, Args, {'SaveJob', 'SkipExist', 'Debug'});
            
            num_workers = Args.NumWorkers;
            task_backup_ratio = Args.TaskBackupRatio;
            disp_prefix = [uedgerun.print_prefix ' '];
            %% check if there are any job files
            job_files = self.job_file_get_available('MaxRequired', num_workers);
            if isempty(job_files)
                disp([disp_prefix 'Exit with no jobs!'])
                return
            end
            %% create parpool
            pool = self.parpool_create(num_workers - 1);
            Args = rmfield(Args, 'NumWorkers');
            Args = rmfield(Args, 'TaskBackupRatio');
            %% initial run
            start_time = now();
            varargin = struct2vararg(Args);
            job_files = self.task_trim_job_files(job_files, task_backup_ratio);
            self.task_enqueue(job_files, varargin{:});
            %% dispathing
            state = 'running'; % 'running', 'finished'
            while(1)
                %% get new tasks
                if self.task_remain(self.tasks) < ...
                        floor(gcp().NumWorkers * (1+task_backup_ratio))
                    job_files = self.job_file_get_available(...
                        'MaxRequired', num_workers, ...
                        'PrintPerformance', false);
                else
                    job_files = [];
                end
                %% check if there are new jobs
                if isempty(job_files) && strcmpi(state, 'finished')
                    disp([disp_prefix 'Exit with no jobs!'])
                    break
                end
                %% add new tasks
                job_files = self.task_trim_job_files(job_files, task_backup_ratio);
                self.task_enqueue(job_files, varargin{:});
                %% display task logs
                if strcmpi(state, 'running')
                    diary off
                    self.task_log_open(start_time, 'PrintNotExist', false);
                    self.task_log_collect();
                    diary on
                end
                %% check if all tasks are finished
                if all(strcmpi({self.tasks.State}, 'finished'))
                    state = 'finished';
                else
                    self.task_remove_finished();
                end
                pause(0.1)
            end
            %% clean
            self.task_log_close();
            delete(pool)
        end
        
        function [splits1, splits2, indbits_diff] = savedt_diff_by_job(self, fpath1, fpath2, scan_name_sequence)
            % much more robust than self.savedt_diff_by_name by slow
            [scan_diff, scan_common] = self.savedt_compare(fpath1, fpath2);
            for i=1:length(scan_name_sequence)
                scan_name = scan_name_sequence{i};
                if isfield(scan_diff, scan_name)
                    splits1{i} = num2str(scan_diff.(scan_name)(1), [scan_name '%g']);
                    splits2{i} = num2str(scan_diff.(scan_name)(2), [scan_name '%g']);
                    indbits_diff(i) = true;
                    continue
                end
                assert(isfield(scan_common, scan_name), ['Unknown scan name "' scan_name '"'])
                scan_str = num2str(scan_common.(scan_name), [scan_name '%g']);
                splits1{i} = scan_str;
                splits2{i} = scan_str;
                indbits_diff(i) = false;
            end
        end
    end
    
    methods
        function self = uedgescan(work_dir, script_input, file_init, varargin)
            %% check arguments
            assert(exist(work_dir, 'dir'), ['Woring dir not exist: ' work_dir])
            self.work_dir = work_dir;
            self.scan_load();
            
            if nargin == 1    
                return
            end
            
            Args.ImageScript = '';
            Args = parseArgs(varargin, Args);
            %% set properties
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
        
        function job = job_create(self, input_diff, file_init, varargin)
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
        
        function job_file = job_save(self, job)
            jobid = job.id;
            job_file = fullfile(self.work_dir, [self.job_file_prefix num2str(jobid) '.mat']);
            save(job_file, 'job');
        end
        
        function status = job_status(self, job, varargin)
            Args.FailDirName = 'fail';
            Args = parseArgs(varargin, Args);
            
            input_diff = job.input_diff;
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
            self.jobinfo = [];
            if contains(lower(Args.InitialRunFileInit), 'near')
                disp([disp_prefix 'Finding initial seed closed to "' self.file_init '"'])
                input_diff_initial = self.scan_find_seed();
                jobid = 0;
                job = self.job_create(input_diff_initial, self.file_init, 'jobid', jobid);
                job_file = self.job_save(job);
                disp([disp_prefix 'Generated: "' job_file '"'])
            end
            %% initial run
            % Given the sequential scan name and value, the remain scan
            % parameters are traversed to start using Args.InitialRunFileInit
            
            scan_start = Args.Scan;
            scan_start.(seq_scan_name).(self.scan_field_value) = seq_start_value;
            
            input_diff_list = self.scan_traverse(scan_start);
            for jobid=1:length(input_diff_list)
                job = self.job_create(input_diff_list{jobid}, Args.InitialRunFileInit, 'jobid', jobid);
                job_file = self.job_save(job);
                disp([disp_prefix 'Generated: "' job_file '"'])
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
                        
                        job_tmp = self.job_create(input_diff, file_init_tmp, 'jobid', jobid);
                        if k==1
                            job = job_tmp;
                        else
                            job(k) = job_tmp;
                        end
                    end
                    job_file = self.job_save(job);
                    disp([disp_prefix 'Generated: "' job_file '"'])
                end
            end
            %% save
            self.scan_save()
        end
        
        function job_gen_nearest(self)
            %% clean
            delete(fullfile(self.work_dir, [self.job_file_prefix '*.mat']))
            disp_prefix = [uedgerun.print_prefix ' '];
            %% pre run for using nearest init file
            self.jobinfo = [];
            disp([disp_prefix 'Finding initial seed closed to "' self.file_init '"'])
            input_diff_initial = self.scan_find_seed();
            jobid = 0;
            job = self.job_create(input_diff_initial, self.file_init, 'jobid', jobid);
            job_file = self.job_save(job);
            disp([disp_prefix 'Generated: "' job_file '"'])
            %% set all other jobs' init file as "nearest"
            input_diff_list = self.scan_traverse(self.scan);
            for jobid=1:length(input_diff_list)
                job = self.job_create(input_diff_list{jobid}, 'nearest', 'jobid', jobid);
                job_file = self.job_save(job);
                disp([disp_prefix 'Generated: "' job_file '"'])
            end
            %% save
            self.scan_save()
        end
        
        function job_gen_rerun(self)
            %% prepare
            disp_prefix = [uedgerun.print_prefix ' '];
            savedt_file_list = self.savedt_file_get('succeed', ...
                'fileextension', uedgerun.file_extension);
            savedt_filenames = {savedt_file_list(:).name};
            %% gen job files
            self.jobinfo = [];
            input_diff_list = self.scan_traverse(self.scan);
            for jobid=1:length(input_diff_list)
                input_diff = input_diff_list{jobid};
                scan_str = uedgerun.input_diff_gen_str(input_diff);
                if ~any(contains(savedt_filenames, scan_str))
                    continue
                end           
                                
                job = self.job_create(input_diff, file_init_tmp, 'jobid', jobid);
                job_file = self.job_save(job);
                disp([disp_prefix 'Generated: "' job_file '"'])
            end
            %% save
            self.scan_save()
        end
        
        function job_files = job_file_get(self, varargin)
            %% job_files = self.job_file_get('ExcludeRunning', false, 'JobCheckMax', inf)
            % This func is also a part of dispatching func. Job checking
            % process is time consuming. Set "JobCheckMax" to skip checking
            % when the number of job files is large.
            
            %% check arguments
            Args.ExcludeRunning = false;
            Args.PrintPerformance = true;
            Args.JobCheckMax = inf;
            
            Args = parseArgs(varargin, Args, {'ExcludeRunning', 'PrintPerformance'});
            disp_prefix = '      '; % not really need to keep for trimming
            is_parallel = ~isempty(gcp('nocreate')); % TODO: cellfun parallel?
            %% get all job files
            if Args.PrintPerformance; tic; end
            job_file_pattern = [self.job_file_prefix, '*.mat'];
            job_files = dir(fullfile(self.work_dir, job_file_pattern));
            if isempty(job_files)
                return
            end
            %% validate job file
            if length(job_files) < Args.JobCheckMax
                job_paths = abspath(job_files, 0);
                inds = cellfun(@(path) self.job_file_check(path, 0), ...
                    job_paths, ... % consumes a lot of time when job number 
                    'UniformOutput', false); % is large. set "JobCheckMax" to skip.
                inds = logical([inds{:}]);
                job_files(~inds) = [];
                job_paths(~inds) = [];
            else
                job_paths = [];
            end
            if Args.PrintPerformance
                fprintf(1, [disp_prefix 'Enumeration of job files: '])
                toc
            end
            %% remove running jobs
            if Args.ExcludeRunning
                if Args.PrintPerformance; tic; end
                if isempty(job_paths)
                    job_paths = abspath(job_files, 0);
                end
                
                lock_files = dir(fullfile(self.work_dir, [job_file_pattern '.lock']));
                if ~isempty(lock_files)
                    lock_job_paths = cellfun(@(path) path(1:end-5) , abspath(lock_files, 0), ...
                        'UniformOutput', false);
                    
                    [~, inds] = setdiff(job_paths, lock_job_paths);
                    job_files = job_files(inds);
                end
                if Args.PrintPerformance
                    fprintf(1, [disp_prefix 'Enumeration of locked files: '])
                    toc
                end
            end
        end
        
        function job_files = job_file_get_available(self, varargin)
            %% job_files = self.job_file_available('MaxRequired', inf)
            % This func is also a part of dispatching func. The init file 
            % enumerating part is time consuming. "MaxRequired" is used to
            % boost this process according the task adding need.
            
            %% check arguments
            Args.MaxRequired = inf;
            Args.SliceLength = 10;
            Args.PrintPerformance = true;
            Args = parseArgs(varargin, Args, {'ExcludeRunning', 'PrintPerformance'});
            disp_prefix = '      '; 
            %% get job files
            job_files = self.job_file_get('ExcludeRunning', ...
                'PrintPerformance', Args.PrintPerformance, ...
                'JobCheckMax', 100);
            assert(length(self.jobinfo) >= length(job_files), '"self.jobinfo" is corrupted!')
            %% get jobid by parsing file name                
            job_file_id_list_cell = cellfun(...
                @(job_name) str2double(regexp(job_name, [uedgescan.job_file_prefix '(\d+)\.mat'], 'tokens', 'once')), ...
                {job_files(:).name}, 'UniformOutput', false);
            job_file_id_list = [job_file_id_list_cell{:}];
            assert( length(job_file_id_list) == length(job_file_id_list_cell), ...
                'Currently does not considering bad job file name, which should not be.')
            %% get init files for job_files
            jobinfo_id_list = [self.jobinfo(:).id];
            [~, inds_map] = ismember(job_file_id_list, jobinfo_id_list); % map index from self.jobinfo to job_files
            file_init_list  = {self.jobinfo(inds_map).file_init};
            input_diff_list = {self.jobinfo(inds_map).input_diff};
            %% select job files based on the type of init file
            inds_bad = ~contains(file_init_list, 'near') & ...
                ~contains(file_init_list, uedgerun.file_save_prefix);
            job_files(inds_bad) = [];
            file_init_list(inds_bad) = [];
            input_diff_list(inds_bad) = [];
            inds_job = 1:length(job_files);
            %% init file is specified as "savedt*hdf5"
            if Args.PrintPerformance; tic; end
            indbits_savedt = contains(file_init_list, uedgerun.file_save_prefix);
            inds_job_savedt = inds_job(indbits_savedt);
            
            indbits_inds_job_savedt = cellfun(@(fname) exist(fname, 'file') == 2, ...
                file_init_list(indbits_savedt));
            inds_avail = inds_job_savedt(indbits_inds_job_savedt);
            
            if Args.PrintPerformance
                fprintf(1, [disp_prefix 'Enumeration of init files ("savedt*hdf5"): '])
                toc
            end
            
            if length(inds_avail) >= Args.MaxRequired
                job_files = job_files(inds_avail);
                job_files = job_files(1:Args.MaxRequired);
                return
            end
            %% init file is "nearest"
            if Args.PrintPerformance; tic; end
            indbits_near = contains(file_init_list, 'near');
            inds_job_near = inds_job(indbits_near);
            
            slice_no = 1;
            if length(inds_job_near) > Args.SliceLength
                slice_no = ceil(length(inds_job_near)/Args.SliceLength);
            end
            
            for i=1:slice_no
                if slice_no == 1
                    inds_part = 1:length(inds_job_near);
                else
                    if i == slice_no
                        inds_part = (i-1)*Args.SliceLength:length(inds_job_near);
                    else
                        inds_part = (i-1)*Args.SliceLength + (1:Args.SliceLength);
                    end
                end
                
                inds_job_near_part = inds_job_near(inds_part);
                indbits_inds_job_near_part = cellfun(@(inpdiff) ...
                    ~isempty(self.file_init_find_nearest(self.work_dir, self.scan, inpdiff)), ...
                    input_diff_list(inds_job_near_part));
                
                inds_avail = unique([inds_avail inds_job_near_part(indbits_inds_job_near_part)], 'stable');
                if length(inds_avail) >= Args.MaxRequired
                    break
                end
            end
            if Args.PrintPerformance
                fprintf(1, [disp_prefix 'Enumeration of init files ("nearest"): '])
                toc
            end
            
            job_files = job_files(inds_avail);
            if length(inds_avail) >= Args.MaxRequired
                job_files = job_files(1:Args.MaxRequired);
            end
        end
        
        function job_file_update(self)
            disp_prefix = [uedgerun.print_prefix ' '];
            job_files = self.job_file_get();
            for i=1:length(job_files)
                job_file = abspath(job_files(i));
                jobs = matread(job_file, 'job');
                %% check jobs
                index_pending = [];
                index_done = [];
                for j=1:length(jobs)
                    job = jobs(j);
                    status = self.job_status(job);
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
        
        function file_profiles = savedt_file_get(self, status, varargin)
            Args.FileExtension = '.mat';
            Args = parseArgs(varargin, Args);
            
            succeded_suffix = ['succe*' Args.FileExtension];
            failed_suffix = ['fail*' Args.FileExtension];
            
            switch lower(status)
                case {'succeed', 'succeeded', 'successful'}
                    file_profiles = dir(fullfile(self.work_dir, [uedgerun.file_save_prefix '*' succeded_suffix]));
                case {'fail', 'failed'}
                    file_profiles = dir(fullfile(self.work_dir, 'fail', [uedgerun.file_save_prefix '*' failed_suffix]));
                otherwise
                    error(['Unknown status: "' status '"'])
            end
        end
                
        function diary_file = run(self, varargin)
            %% diary_file = self.run( ...
            %      'SkipExist', true, ...
            %      'SaveJob', true, ...
            %      'FailDirName', 'fail', ...
            %      'Debug', false, ...
            %      'UseParallel', true, ...
            %      'NumWorkers', maxNumCompThreads, ...
            %      'LogFileName', 'uedgescan_log.txt')
            
            %% check arguments
            Args.SkipExist = true;
            Args.SaveJob = true;
            Args.FailDirName = 'fail';
            Args.Debug = false;
            Args.UseParallel = true;
            Args.NumWorkers = maxNumCompThreads;
            Args.LogFileName = 'uedgescan_log.txt';
            Args = parseArgs(varargin, Args, {'SaveJob', 'SkipExist', 'Debug', 'UseParallel'});
            
            use_parallel = Args.UseParallel;
            diary_file = self.log_start(Args.LogFileName);
            disp([uedgerun.print_prefix ' Current dir is "' pwd '"'])
            disp([uedgerun.print_prefix ' Working on "' self.work_dir '"'])
            %% parallel run            
            Args = rmfield(Args, 'UseParallel');
            Args = rmfield(Args, 'LogFileName');
            varargin = struct2vararg(Args);
            
            if use_parallel
                self.run_parallel(varargin{:})
                self.log_end();
                return
            end
            %% normal run
            Args = rmfield(Args, 'NumWorkers');
            varargin = struct2vararg(Args);
            
            self.run_normal(varargin{:});
            self.log_end();
        end
        
        function fail_files = job_file_fail_gather(self, fail_file_or_dir)
            %% single file
            if ischar(fail_file_or_dir)
                fail_files = self.savedt_find_fail(fail_file_or_dir);
                return
            end
            %% multiple files
            files = fail_file_or_dir;
            [fail_files, cache] = self.savedt_find_fail(abspath(files(1)));
            for i=2:length(files)
                fail_files_tmp = self.savedt_find_fail(abspath(files(i)), cache);
                if length(fail_files_tmp) > 1
                    disp(['Found multiple:' files(i).name])
                end
                fail_files(end+1:end+length(fail_files_tmp)) = fail_files_tmp;
            end
        end
        
        function file_init_disp(self, varargin)   
            %% self.file_init_disp('MinDiff', 0, 'Robust', false)   
            
            %% check arguments
            Args.MinDiff = 0;
            Args.Robust = false;
            Args = parseArgs(varargin, Args, 'Robust');
            %% get savedt files
            savedt_files = abspath(...
                dir( fullfile(self.work_dir, [uedgerun.file_save_prefix '*.mat']) ),...
                0);
            %% get delimiters
            delimiters = {'-'};
            scan_names = sort(fieldnames(self.scan));
            delimiters(2:length(scan_names)+1) = scan_names;
            %% main
            for i=1:length(savedt_files)
                %% get files
                fsave = savedt_files{i};
                job = matread(fsave, 'job');
                lines = {'---------------------------------------------------'};
                if isfield(job, 'file_init_used')
                    finit = job.file_init_used;
                else
                    finit = job.file_init;
                end
                %% disp file names
                [~, finit_name] = fileparts(finit);
                [~, fsave_name] = fileparts(fsave);
                if strcmp(finit_name, fsave_name)
                    warning(['Itself: ' fsave_name])
                    continue
                end
                lines{end+1} = ['Init: ' finit_name];
                lines{end+1} = ['Save: ' fsave_name];
                %% disp differences
                % split input_diff str
                scan_str = strsplit(finit_name, '_'); 
                if length(scan_str) < 2
                    continue
                end
                if ~Args.Robust
                    scan_str = scan_str{2};
                    init_splits = strsplit(scan_str, delimiters);
                    scan_str = strsplit(fsave_name, '_'); scan_str = scan_str{2};
                    save_splits = strsplit(scan_str, delimiters);
                end
                % collect diff
                job_file_init = fullfile(self.work_dir, [finit_name '.mat']);
                job_file_save = fullfile(self.work_dir, [fsave_name '.mat']);
                if ~Args.Robust && length(init_splits) == length(save_splits)
                    indbits_diff = cellfun(@(s1, s2) ~strcmp(s1, s2), init_splits, save_splits);
                    if length(delimiters) ~= length(indbits_diff)
                        [init_splits, save_splits, indbits_diff] = ...
                            self.savedt_diff_by_job(...
                            job_file_init, ...
                            job_file_save, ...
                            scan_names);
                    end
                else
                    [init_splits, save_splits, indbits_diff] = ...
                        self.savedt_diff_by_job(...
                        job_file_init, ...
                        job_file_save, ...
                        scan_names);
                end
                if sum(indbits_diff) < Args.MinDiff
                    continue
                end
                % gen diff line
                tmp_deli = delimiters([false indbits_diff]);
                tmp_init = init_splits(indbits_diff);
                tmp_save = save_splits(indbits_diff);
                for j=1:length(tmp_deli)
                    lines{end+1} = ['Diff: ' tmp_deli{j} ...
                        '(' tmp_init{j} ' -> ' tmp_save{j} ')'];
                end
                % print
                for j=1:length(lines)
                    line = lines{j};
                    if contains(line, 'Diff')
                        fprintf(2, [line '\n'])
                    else
                        disp(line)
                    end
                end
                if sum(indbits_diff) > 1
                    fprintf(2, 'Not nearest init is used!\n')
                end
            end
        end
        
        function inspect(self, varargin)
            %% check arguments
            Args.Status = {'successful', 'failed', 'pending'};
            Args.Verbose = false;
            Args = parseArgs(varargin, Args, {'Verbose'});
            
            status_filter = Args.Status;
            %% disp job status
            job_files = self.job_file_get();
            for i=1:length(job_files)
                f = job_files(i);
                job = matread(f, 'job');
                disp(['***** ' f.name ' *****'])
                for j=1:length(job)
                    str_save = uedgerun.input_diff_gen_str(job(j).input_diff);
                    [~, profile_init] = fileparts(job(j).file_init);
                    status = self.job_status(job(j));
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
            job_list_succeeded = self.savedt_file_get('succeeded', ...
                'fileextension', uedgerun.file_extension);
            len_succeeded = length(job_list_succeeded);
            disp(['Succeded cases: ' num2str(len_succeeded)])
            job_list_failed = self.savedt_file_get('failed', 'fileextension', '.mat');
            len_failed = length(job_list_failed);
            disp(['Failed cases: ' num2str(len_failed)])
            disp(['Remain cases: ' num2str(len_all - len_succeeded - len_failed)])
        end
             
        function clean(self)
            %% clean run and input files
            uedgerun.clean()
            %% clean not deleted job files
            self.job_file_update()
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
        
        function varargout = time_stat(self, varargin)
            %% time = self.time_stat('TimeRange', [], ...
            %    'ScanPattern', '*', 'Plot', true, ...
            %    'Print', true)
            
            %% check arguments
            Args.TimeRange = [0 inf];
            Args.ScanPattern = '*';
            Args.Plot = true;
            Args.Print = true;
            Args = parseArgs(varargin, Args, {'Plot', 'Print'});
            %% get elapsed time
            file_list = dir(...
                fullfile(self.work_dir, ...
                [uedgerun.file_save_prefix Args.ScanPattern 'successful.mat']));
            time = zeros(1, length(file_list));
            for i=1:length(file_list)
                f = file_list(i);
                job = matread(f, 'job');
                time(i) = job.elapsed_time;
            end
            %% slice
            time = time/60;
            inds = time < min(Args.TimeRange) | time > max(Args.TimeRange);
            file_list(inds) = [];
            time(inds) = [];
            if nargout == 1
                varargout = {time};
            end
            %% plot
            if Args.Plot
                uedgedata.figure;
                histogram(time);
                xlabel('Elapsed Time [minutes]')
                ylabel('Counts')
                uedgedata.figure_decoration;
            end
            %% print
            for i=1:length(file_list)
                if time(i) > 60; fid = 2; else; fid = 1; end
                fprintf(fid, '%s: %.1f minutes\n', ...
                    uedgerun.extract_input_diff_str(file_list(i).name), time(i));
            end
        end
    end
end