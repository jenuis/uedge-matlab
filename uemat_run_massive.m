function uemat_run_massive(work_dir, cluster_name, num_workers, log_file)
    %% check arguments
    if nargin < 4
        log_file = ['pool_' work_dir '_log.txt'];
    end
    
    if nargin < 3
        num_workers = maxNumCompThreads - 1;
        disp([uedgerun.print_prefix ' Using maxNumCompThreads(' num2str(maxNumCompThreads) ')-1'])
    else
        if ischar(num_workers)
            num_workers = str2double(num_workers);
        end
        num_workers = floor(num_workers);
        if isnan(num_workers)
            num_workers = maxNumCompThreads - 1;
            disp([uedgerun.print_prefix ' Changing NumWorkers to maxNumCompThreads(' num2str(maxNumCompThreads) ')-1'])
        end
    end
    
    if nargin < 2
        cluster_name = 'local';
    end
    %% check and close gcp
    p = gcp('nocreate');
    if ~isempty(p)
        delete(p)
    end
    %% open parallel pool
    cluster_profile = parcluster(cluster_name);
    cluster_profile.NumWorkers = num_workers;
    cluster_profile.NumThreads = 1;
    parpool(cluster_profile, num_workers);
    %% call uedgescan.run
    usc = uedgescan(work_dir);
    usc.run('LogFileName', log_file);
end