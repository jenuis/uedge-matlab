% Author: Xiang LIU@ASIPP
% E-mail: xliu@ipp.ac.cn
% Created: 2023-12-16
% Version: V 0.1.11
% TODO: make '>>>>> ' a constant across the project
classdef uedgescan < handle    
    properties
        work_dir
        scan = struct()
        script_input
        script_image
        file_init
        job_file_prefix = 'job';
    end
    methods(Static)
        function input_diff_list = scan_traverse(scan, value_field_name)
            %% check arguments
            if nargin < 2
                value_field_name = 'value';
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
                value_field_name = 'value';
            end
            
            for i=1:length(input_diff_list)
                disp(uedgerun.input_diff_gen_str(input_diff_list{i}, 'ValueFieldName', value_field_name))
            end
        end
        
        function job = job_run(job, varargin)
            %% check arguments
            Args.PrintPrefix = '>>>>> ';
            Args.SkipExist = true;
            Args.ID = [];
            Args.SaveJob = true;
            Args.FailDirName = 'fail';
            Args = parseArgs(varargin, Args, {'SaveJob', 'SkipExist'});
            
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
            
            disp_prefix = Args.PrintPrefix;
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
                disp([disp_prefix 'Skip failed job: ' job_file_fail])
                return
            end
            %% check init file
            flag_fail = ~exist(job.file_init, 'file');
            if flag_fail
                file_init_corrected = fullfile(work_dir, job.file_init);
                flag_fail = ~exist(file_init_corrected, 'file');
                if ~flag_fail
                    job.file_init = file_init_corrected;
                end
            end
            %% run job            
            if flag_fail
                %% omit run
                job.status = -1; % no init file
                disp([disp_prefix 'No init file: ' job.file_init])
            else
                %% set input_diff
                ur = uedgerun(job.script_input, job.file_init, 'ImageScript', script_image);

                fnames = fieldnames(job.input_diff);
                for i=1:length(fnames)
                    fn = fnames{i};
                    field = job.input_diff.(fn);
                    ur.input_diff_update(fn, field.value, field.uedgevar.names, field.uedgevar.amps);
                end
                %% check existence
                ur.file_save = fullfile(work_dir, ur.generate_file_name(ur.input_diff));
                if Args.SkipExist && exist(ur.file_save, 'file')
                    disp([disp_prefix 'Skip exist: ' ur.file_save])
                    return
                end
                %% run
                ur.script_run_gen('ID', Args.ID, 'PrintPrefix', Args.PrintPrefix);
                tic
                job.status = ur.run();
                job.elapsed_time = toc;
                disp([disp_prefix 'Total time used: ' num2str(job.elapsed_time/60) ' minutes!'])
            end
            %% save job
            if ~Args.SaveJob
                return
            end
            
            if job.status == 0
                job_file = strrep(ur.file_save, ur.file_extenstion, '.mat');
                run_status = 'successful';
            else
                job_file = job_file_fail;
                run_status = 'failed';
            end
            save(job_file, 'job')
            disp([disp_prefix 'Converged status: ' upper(run_status)])
        end
        
        function time_stat(work_dir)
            ur = uedgerun();
            file_list = dir(fullfile(work_dir, [ur.file_save_prefix '*.mat']));
            time = zeros(length(file_list));
            for i=1:length(file_list)
                f = file_list(i);
                job_file = abspath(f);
                mf = matfile(job_file);
                job = mf.job;
                time(i) = job.elapsed_time;
            end
            uedgedata.figure;
            histogram(time/60);
            xlabel('Elapsed Time [minutes]')
            ylabel('Counts')
            uedgedata.figure_decoration;
        end
    end
    methods(Access=private, Static)
        function status = run_lock(lock_file, mode)
            mode = lower(mode);
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
        
        function run_job_files(work_dir, job_files, index, varargin)
            %% get job file
            fpath = abspath(job_files(index));
            mf = matfile(fpath);
            %% check if the job is running
            if uedgescan.run_lock(fpath, 'check')
                disp(['>>>>> Skipping running job: "' job_files(index).name '"!'])
                pause(0.1)
                return
            end
            %% check and wait init file
            job = mf.job;
            job_file_fail = uedgerun.generate_file_name(job(1).input_diff, 'suffix', 'fail', 'extension', '.mat');
            job_file_fail = fullfile(work_dir, 'fail', job_file_fail);
            tic
            if ~exist(job(1).file_init, 'file') && ~exist(job_file_fail, 'file')
                disp(['[' num2str(job(1).id) ']>>>>> Waiting init file (' num2str(toc) ' s): ' job(1).file_init])
                pause(0.1)
                return
            end
            %% run job
            uedgescan.run_lock(fpath, 'lock');
            for i=1:length(job)
                uedgescan.job_run(job(i), varargin{:});
            end
            %% delete job file
            disp(['>>>>> Completed: "' fpath '"'])
            delete(fpath)
            uedgescan.run_lock(fpath, 'unlock');
        end
    end
    methods
        function self = uedgescan(work_dir, script_input, file_init, varargin)
            %% check arguments
            if nargin == 1
                self.work_dir = work_dir;
                self.scan_load();
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
        end
        
        function scan_load(self)
            file = ['scan_' self.work_dir '.mat'];
            assert(exist(file, 'file'), ['"' file '" not exist!'])
            mf = matfile(file);
            file_self = mf.self;
            fnames = fieldnames(file_self);
            for i=1:length(fnames)
                fn = fnames{i};
                self.(fn) = file_self.(fn);
            end
        end
        
        function scan_save(self)
            file = ['scan_' self.work_dir '.mat'];
            save(file, 'self')
        end
        
        function scan_add(self, scan_name, scan_values, uedgevar_names, uedgevar_amps)
            %% scan_add(self, scan_name, scan_values, uedgevar_names, uedgevar_amps)
            
            %% set properties
            self.scan.(scan_name).value = unique(sort(scan_values));
            self.scan.(scan_name).uedgevar.names = uedgevar_names;
            self.scan.(scan_name).uedgevar.amps = uedgevar_amps;
            %% save
            self.scan_save();
        end
                
        function job = job_struct(self, input_diff, file_init, varargin)
            Args.ScriptImage = self.script_image;
            Args.ScriptInput = self.script_input;
            Args.WorkDir = self.work_dir;
            Args.JobID = [];
            Args = parseArgs(varargin, Args);
            
            uedgerun.input_diff_check(input_diff);
            
            job = struct();               
            job.input_diff = input_diff;
            job.file_init = file_init;
            job.script_input = Args.ScriptInput;
            job.script_image = Args.ScriptImage;
            job.work_dir = Args.WorkDir;
            job.id = Args.JobID;
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
            Args.ValueFieldName = 'value';
            Args = parseArgs(varargin, Args);
            %% get seq info
            seq = self.scan.(seq_scan_name);
            seq_values = seq.(Args.ValueFieldName);
            seq_start_value = seq_values(findvalue(seq_values, seq_start_value));
            %% clean
            delete(fullfile(self.work_dir, [self.job_file_prefix '*.mat']))
            %% run in the 1st place
            scan_start = self.scan;
            scan_start.(seq_scan_name).(Args.ValueFieldName) = seq_start_value;
            
            input_diff_list = self.scan_traverse(scan_start, Args.ValueFieldName);
            for jobid=1:length(input_diff_list)
                job = self.job_struct(input_diff_list{jobid}, self.file_init, 'jobid', jobid);
                job_file = self.job_name_gen(jobid);
                disp(['>>>>> Generating: "' job_file '"'])
                save(job_file, 'job')
            end
            %% run sequence
            scan_remain = rmfield(self.scan, seq_scan_name);
            input_diff_remain_list = self.scan_traverse(scan_remain, Args.ValueFieldName);
            
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
                        input_diff.(seq_scan_name).(Args.ValueFieldName) = value_init;
                        file_init_tmp = fullfile(self.work_dir, uedgerun.generate_file_name(input_diff));
                        
                        value = seq_values_tmp(k);
                        input_diff.(seq_scan_name).(Args.ValueFieldName) = value;
                        
                        job_tmp = self.job_struct(input_diff, file_init_tmp, 'jobid', jobid);
                        if k==1
                            job = job_tmp;
                        else
                            job(k) = job_tmp;
                        end
                    end
                    job_file = self.job_name_gen(jobid);
                    disp(['>>>>> Generating: "' job_file '"'])
                    save(job_file, 'job')
                end
            end
        end
        
        function job_gen_rerun(self)
            savedt_file_list = self.profile_files_get('succeed', 'fileextension', '.hdf5');
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
                disp(['>>>>> Generating: "' job_file '"'])
                save(job_file, 'job')
            end
        end
        
        function job_files = job_get_files(self, no_running)
            %% get all jobs
            job_files = dir(fullfile(self.work_dir, [self.job_file_prefix, '*.mat']));
            if nargin < 2 || ~no_running
                return
            end
            %% remove running jobs
            inds_running = [];
            for i=1:length(job_files)
                fpath = abspath(job_files(i));
                lock_file = strrep(fpath, '.mat', '.mat.lock');
                if exist(lock_file, 'file')
                    inds_running = i;
                end
            end
            job_files(inds_running) = [];
        end
        
        function job_update(self)
            job_files = self.job_get_files();
            for i=1:length(job_files)
                job_file = abspath(job_files(i));
                mf = matfile(job_file);
                jobs = mf.job;
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
                    disp(['>>>>> Deleting completed: "' job_file '"'])
                    delete(job_file)
                    continue
                end                
            end
        end
        
        function job_disp(self, status_filter)
            if nargin < 2
                status_filter = {'successful', 'failed', 'pending'};
            end
            
            job_files = self.job_get_files();
            for i=1:length(job_files)
                f = job_files(i);
                fpath = abspath(f);
                mf = matfile(fpath);
                disp(['***** ' f.name ' *****'])
                job = mf.job;
                for j=1:length(job)
                    str_save = uedgerun.input_diff_gen_str(job(j).input_diff);
                    [~, profile_init] = fileparts(job(j).file_init);
                    status = self.job_status(job(j).input_diff);
                    if ~isempty(status_filter) && any(contains(lower(status_filter), lower(status(1:4))))
                        disp(['"' str_save '" <-- "' profile_init '", "' upper(status) '"'])
                    end                    
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
        
        function job_remain(self, varargin)
            job_list_all = self.scan_traverse(self.scan);
            len_all = length(job_list_all);
            disp(['>>>>> All cases: ' num2str(len_all)])
            job_list_succeeded = self.profile_files_get('succeeded', 'fileextension', '.hdf5');
            len_succeeded = length(job_list_succeeded);
            disp(['>>>>> Succeded cases: ' num2str(len_succeeded)])
            job_list_failed = self.profile_files_get('failed', 'fileextension', '.mat');
            len_failed = length(job_list_failed);
            disp(['>>>>> Failed cases: ' num2str(len_failed)])
            disp(['>>>>> Remain cases: ' num2str(len_all - len_succeeded - len_failed)])
        end
                
        function run(self, varargin)
            %% check arguments
            Args.PrintPrefix = '>>>>> ';
            Args.SkipExist = true;
            Args.ID = [];
            Args.SaveJob = true;
            Args.FailDirName = 'fail';
            Args.UseParallel = true;
            Args = parseArgs(varargin, Args, {'SaveJob', 'SkipExist', 'UseParallel'});
            
            use_parallel = Args.UseParallel;
            dir_work = self.work_dir;
            
            Args = rmfield(Args, 'UseParallel');
            varargin = struct2vararg(Args);
            %% normal run
            if ~use_parallel
                while(1)
                    job_files = self.job_get_files(1);
                    if isempty(job_files)
                        disp('>>>>> Exit with no jobs!')
                        break
                    end

                    for index=1:length(job_files)
                        uedgescan.run_job_files(dir_work, job_files, index, varargin{:})
                    end
                end
                return
            end
            %% parallel run
            while(1)
                job_files = self.job_get_files(1);
                if isempty(job_files)
                    disp('>>>>> Exit with no jobs!')
                    break
                end
                
                parfor index=1:length(job_files)
                    uedgescan.run_job_files(dir_work, job_files, index, varargin{:})
                end
            end
        end
        
        function clean(self)
            uedgerun.clean()
            lock_files = dir(fullfile(self.work_dir, 'job*.lock'));
            for i=1:length(lock_files)
                lock_file = abspath(lock_files(i));
                disp(['Remove: "' lock_file '"'])
                delete(lock_file)
            end
        end
    end
end