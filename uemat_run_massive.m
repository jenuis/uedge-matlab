function uemat_run_massive(work_dir, cluster_name, num_workers, log_file)
    %% check arguments
    if nargin < 4
        log_file = ['pool_' work_dir '_log.txt'];
    end
    
    if nargin < 3
        num_workers = maxNumCompThreads - 1;
        disp([uedgerun.print_prefix ' Using maxNumCompThreads(' num2str(num_workers) ')-1'])
    else
        num_workers = str2double(num_workers);
        num_workers = floor(num_workers);
        if isnan(num_workers) || num_workers >= maxNumCompThreads
            num_workers = maxNumCompThreads - 1;
            disp([uedgerun.print_prefix ' Changing NumWorkers to maxNumCompThreads(' num2str(num_workers) ')-1'])
        end
    end
    
    if nargin < 2
        cluster_name = 'local';
    end
    %% open parallel pool
    cluster = parcluster(cluster_name);
    cluster.NumWorkers = num_workers;
    for i=1:10
        try
            parpool(cluster, num_workers);
            break
        catch ME
            pause(30)
            disp([uedgerun.print_prefix ' Failed to start parpool for ' num2str(i) ' times!'])
            disp(ME.message)
            if i == 10
                disp([uedgerun.print_prefix ' Exit with error: no parallel pool!'])
                return
            end
        end
    end
    %% call uedgescan.run
    usc = uedgescan(work_dir);
    usc.run('LogFileName', log_file);
end