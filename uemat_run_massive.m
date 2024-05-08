function uemat_run_massive(work_dir, cluster_name_or_num_workers, log_file)
    %% check arguments
    if nargin < 3
        log_file = ['pool_' work_dir '_log.txt'];
    end
    if nargin < 2
        num_workers = maxNumCompThreads - 1;
        disp([uedgerun.print_prefix ' Using maxNumCompThreads(' num2str(num_workers) ')-1'])
    else
        num_workers = str2double(cluster_name_or_num_workers);
        num_workers = floor(num_workers);
        if num_workers >= maxNumCompThreads
            num_workers = maxNumCompThreads - 1;
            disp([uedgerun.print_prefix ' Changing NumWorkers to maxNumCompThreads(' num2str(num_workers) ')-1'])
        end
    end
    %% open parallel pool
    if ~isnan(num_workers)
        parpool(num_workers);
    else
        cluster_name = cluster_name_or_num_workers;
        parpool(cluster_name);
    end
    %% call uedgescan.run
    usc = uedgescan(work_dir);
    usc.run('LogFileName', log_file);
end