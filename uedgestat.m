% Author: Xiang LIU@ASIPP
% E-mail: xliu@ipp.ac.cn
% Created: 2023-10-11
% Version: V 0.1.11
classdef uedgestat < handle
    properties(Access=protected)
        work_dir
        files
    end

    properties
        scan_paras
        data_1d = struct()
        data_0d = struct()
    end
    
    methods(Static)        
        function phy_name = translate_phy_name(phy_name)
            phy_name = strrep(phy_name, '(', '');
            phy_name = strrep(phy_name, ')', '');
        end
        
        function val = get_struct_value(s, uri)
            uri = uedgedata.check_uri(uri);
            fnames = strsplit(uri(2:end), '/');
            for i=1:length(fnames)
                fn = fnames{i};
                assert(isfield(s, fn), ['field "' fn '" not exist!'])
                s = s.(fn);
            end
            val = s;
        end
        
        function para = scan_parse_file(work_dir, profile_file_name, varargin)
            %% check arguments
            Args.ValueFieldName = 'value';
            Args = parseArgs(varargin, Args);
            
            para = struct();
            %% parse scan para using job record
            job_file = strrep(profile_file_name, uedgerun.file_extension, '.mat');
            job_file = fullfile(work_dir, job_file);
            assert( exist(job_file, 'file'), ['Job file missing: ' job_file])
            job = matread(job_file, 'job');                    
            fnames = fieldnames(job.input_diff);
            for i=1:length(fnames)
                fn = fnames{i};
                para.(fn) = job.input_diff.(fn).(Args.ValueFieldName);
            end
        end
    end
    
    methods
        function self = uedgestat(work_dir)
            addpath_mdslab('divlp');            
            assert(exist(work_dir, 'dir'), '"work_dir" does not exist!')
            self.work_dir = work_dir;
        end
        
        function files = get_files(self)
            %% check and return
            if ~isempty(self.files)
                files = self.files;
                return
            end
            %% list files
            files = dir(fullfile(self.work_dir, [uedgerun.file_save_prefix '*' uedgerun.file_extension]));
            %% filter files
            if isempty(gcp('nocreate'))
                inds = cellfun(@(f) uedgedata.is_uedge_file(f), abspath(files));
            else
                inds = false(size(files));
                parfor i=1:length(files)
                    inds(i) = uedgedata.is_uedge_file(abspath(files(i)));
                end
            end
            files = files(inds);
            %% use relative path
            for i=1:length(files)
                files(i).folder = self.work_dir;
            end
            %% set properties
            self.files = files;
        end
               
        function file_path = get_profile_path(self, varargin)
            %% check arguments
            para_len = length(varargin);
            if para_len > 1
                assert(mod(para_len, 2)==0, 'The input arguments are not paired!')
                filter = struct();
                for i=1:(para_len/2)
                    filter.(varargin{2*i-1}) = varargin{2*i};
                end
            else
                filter = varargin{1};
            end
            %% get by file index
            if isnumeric(filter) && filter > 0 && filter <= length(self.files)
                file_path = abspath(self.files(filter));
                if ~exist(file_path, 'file')
                    file_path = fullfile(self.work_dir, self.files(filter).name);
                end
                return
            end
            %% get by para, recursive call
            assert(isstruct(filter), 'Input "filter" should be a struct!')
            assert(length(self.filter_remain(filter)) < 2, 'Input "filter" can only omit at most one scan parameter!')
            inds = self.filter_select(filter);
            if ~islogical(inds)
                file_path = {};
                for i=1:length(inds)
                    file_path{i} = self.get_profile_path(inds(i));
                end
                return
            end
            %% get by para, single call
            file_path = [];
            for i=1:length(self.files)
                file_name = self.files(i).name;
                para_tmp = self.scan_parse_file(self.work_dir, file_name);
                if isequal(para_tmp, filter)
                    file_path = abspath(self.files(i));
                    if ~exist(file_path, 'file')
                        file_path = fullfile(self.work_dir, file_name);
                    end
                    return
                end
            end
        end
        
        function ud = get_uedge_data(self, file_index)
            ud = uedgedata(self.get_profile_path(file_index));
        end
        
        function job_file = get_job_file(self, file_index)
            file_path = self.get_profile_path(file_index);
            job_file = strrep(file_path, uedgerun.file_extension, '.mat');
        end
        
        function rerun(self, file_index)
            job_file = self.get_job_file(file_index);
            job = matread(job_file, 'job');
            job.file_init = self.get_profile_path(file_index);
            uedgescan.run_basic(job, ...
                'SkipExist', false, ...
                'SaveJob', false, ...
                'FailDirName', '../');
        end
        
        function image_check(self)
            %% check if image file is OK
            savedt_files = self.files;
            if isempty(gcp('nocreate'))
                inds = cellfun(@(f) uedgedata.is_h5(...
                    strrep(f, uedgerun.file_save_prefix, uedgerun.file_image_prefix)), ...
                    abspath(savedt_files));
            else
                inds = false(size(savedt_files));
                parfor i=1:length(savedt_files)
                    inds(i) = uedgedata.is_h5(...
                        strrep(abspath(savedt_files(i)), uedgerun.file_save_prefix, uedgerun.file_image_prefix)...
                        );
                end
            end
            %% repair image file
            inds_repair = 1:length(savedt_files);
            inds_repair = inds_repair(~inds);
            for i=1:length(inds_repair)
                file_index = inds_repair(i);
                disp([uedgerun.print_prefix ' Repairing image for: "' savedt_files(file_index).name '"'])
                self.rerun(file_index);
            end
        end
        
        function scan_parse(self)
            %% check
            if ~isempty(self.scan_paras)
                return
            end
            %% get files
            all_files = self.get_files();
            self.scan_paras = struct();
            workdir = self.work_dir;
            %% get paras from savedt job file
            paras = cell(size(all_files));
            if isempty(gcp('nocreate'))
                for i=1:length(all_files)
                    paras{i} = uedgestat.scan_parse_file(workdir, all_files(i).name);
                end
            else
                parfor i=1:length(all_files)
                    paras{i} = uedgestat.scan_parse_file(workdir, all_files(i).name);
                end
            end
            %% set self.scan_paras
            for i=1:length(paras)
                para = paras{i};
                fnames = fieldnames(para);
                for j=1:length(fnames)
                    var_name = fnames{j};
                    self.scan_paras.(var_name)(i) = para.(var_name);
                end
            end
        end
        
        function para_names = scan_get_names(self)
            para_names = fieldnames(self.scan_paras);
        end
              
        function vals = scan_get_values(self, para_name)
            self.scan_parse();
            vals = self.scan_paras.(para_name);
        end
        
        function para_vals = scan_get_space(self, para_name)
            self.scan_parse();
            
            if nargin < 2
                para_names = self.scan_get_names();
                for i=1:length(para_names)
                    disp(['----- ' para_names{i} ' -----'])
                    disp(self.scan_get_space(para_names{i}))
                end
                return
            end
            
            para_vals = unique(self.scan_get_values(para_name));
        end
        
        function comb = scan_paras_combiner(self, scan_names, convert2cell)
            if nargin < 3
                convert2cell = false;
            end
            if nargin < 2
                scan_names = fieldnames(self.scan_paras);
            end
            
            num_converge = length(self.scan_paras.(scan_names{1}));
            comb = zeros(num_converge, length(scan_names));
            for i=1:num_converge
                comb(i,:) = cellfun(@(name) self.scan_paras.(name)(i), scan_names);
            end
            
            if convert2cell
                comb = mat2cell(comb, ones(1, size(comb, 1)), length(scan_names));
            end
        end
        
        function vals = collect_1d(self, phy_name, poloidal_location)
            %% get file_list
            file_list = self.get_files();
            %% check arguments
            phy_name = lower(phy_name);
            phy_name_attr = self.translate_phy_name(phy_name);
            
            ud = self.get_uedge_data(1);
            phy = ud.read_physical(phy_name, 'warningoff');
            flag_2d = dimnum(phy.data) == 2;
            if ~flag_2d
                flag_1d = dimnum(phy.data) == 1;
                assert(flag_1d, '"phy_name" is not a 1D variable!');
                flag_ot = strcmpi(poloidal_location, 'lo') && phy_name(end) == 'r';
                flag_it = strcmpi(poloidal_location, 'li') && phy_name(end) == 'l';
                assert(flag_ot || flag_it, '"phy_name" is not a valid target variable!');
            end
            %% return if exist
            pol_name = num2str(poloidal_location);
            flag_exist = ...
                isfield(self.data_1d, phy_name_attr) && ...
                isfield(self.data_1d.(phy_name_attr), pol_name) && ...
                ~isempty(self.data_1d.(phy_name_attr).(pol_name));
            if flag_exist
                vals = self.data_1d.(phy_name_attr).(pol_name);
                return
            end
            %% collect
            vals = nan(length(file_list), size(ud.rm,2));
            parfor i=1:length(file_list)
                ud = self.get_uedge_data(i);
                ind_poloidal = ud.get_poloidal_index(poloidal_location);
                phy = ud.read_physical(phy_name, 'warningoff');
                if flag_2d
                    vals(i,:) = phy.data(ind_poloidal, :);
                    continue
                end                
                vals(i,:) = phy.data(:);
            end
            self.data_1d.(phy_name_attr).(pol_name) = vals;
        end

        function vals = collect_0d(self, phy_name)
            %% get file_list
            file_list = self.get_files();
            %% check arguments
            phy_name = lower(phy_name);
            ud = self.get_uedge_data(1);
            phy = ud.read_physical(phy_name);
            assert(dimnum(phy.data) == 0, '"phy_name" is not a 0D variable!');
            %% return if exist
            flag_exist = ...
                isfield(self.data_0d, phy_name) && ...
                ~isempty(self.data_0d.(phy_name));
            if flag_exist
                vals = self.data_0d.(phy_name);
                return
            end
            %% collect
            vals = nan(size(file_list));
            parfor i=1:length(file_list)
                ud = self.get_uedge_data(i);
                phy = ud.read_physical(phy_name);
                vals(i) = phy.data;
            end
            self.data_0d.(phy_name) = vals;            
        end
        
        function vals = extract_0d(self, phy_name, poloidal_location, method)
            %% check arguments
            method = lower(method);
            assert(haselement({'lcfs', 'max', 'int', 'lam', 'lamint'}, method), 'Unknow "method"!')
            phy_name_attr = self.translate_phy_name(phy_name);
            %% get data_1d
            vals_1d = self.collect_1d(phy_name, poloidal_location);
            %% return if exist
            pol_name = num2str(poloidal_location);
            flag_exist = ...
                isfield(self.data_0d, phy_name_attr) && ...
                isfield(self.data_0d.(phy_name_attr), pol_name) && ...
                isfield(self.data_0d.(phy_name_attr).(pol_name), method) && ...
                ~isempty(self.data_0d.(phy_name_attr).(pol_name).(method));
            if flag_exist
                vals = self.data_0d.(phy_name_attr).(pol_name).(method);
                return
            end
            %% collect
            vals = nan(size(self.files));
            switch method
                case 'lcfs'
                    %% Value at LCFS
                    ud = self.get_uedge_data(1);
                    if isempty(ud.lcfs_radial_ind)
                        ind_radial = 9;
                        warning('lcfs_radial_ind = 9 is assumed!')
                    else
                        ind_radial = ud.lcfs_radial_ind;
                    end
                    
                    for i=1:size(vals_1d, 1)
                        vals(i) = vals_1d(i, ind_radial);
                    end
                case 'max'
                    %% Max value
                    for i=1:size(vals_1d, 1)
                        vals(i) = max(vals_1d(i,:));
                    end
                case 'int'
                    %% Integral value
                    ud = self.get_uedge_data(1);
                    x = ud.cal_rrsep(poloidal_location, 'RemoveGhost');
                    for i=1:size(vals_1d, 1)
                        y = vals_1d(i,:);
                        y(ud.get_ghost_inds('radial')) = [];
                        vals(i) = trapz(x, y);
                    end
                case 'lam'
                    %% lambda
                    ud = self.get_uedge_data(1);
                    fit_data.xdata = ud.cal_rrsep(poloidal_location);
                    fit_data.xdata = fit_data.xdata'*1e3; %[m] to [mm]
                    dispstat('','init');
                    tot_no = size(vals_1d, 1);
                    for i=1:tot_no
                        fit_data.ydata = vals_1d(i,:);
                        fit_res = prbfit.eichfit(fit_data, 'mp');
                        vals(i) = fit_res.lam;
                        dispstat(['fitting progress: ' num2str(100*i/tot_no,'%5.1f') '%']);
                    end
                    dispstat('','clean');
                case 'lamint'
                    %% lambda integral
                    ud = self.get_uedge_data(1);
                    xdata = ud.cal_rrsep(poloidal_location, 'RemoveGhost');
                    xdata = xdata'*1e3; %[m] to [mm]
                    for i=1:size(vals_1d, 1)
                        ydata = vals_1d(i,:);
                        ydata(ud.get_ghost_inds('radial')) = [];
                        inds = ydata > 0;
                        if sum(inds) < 3
                            ydata = -ydata;
                            warning(['Inversed: ' phy_name])
                            inds = ydata > 0;
                        end
                        fit_data.xdata = xdata(inds);
                        fit_data.ydata = ydata(inds);
                        vals(i) = prbfit.cal_lambda_int(fit_data, [], 'type', 'raw');
                    end
            end
            
            self.data_0d.(phy_name_attr).(pol_name).(method) = vals;
        end
        
        function extract_0d_batch(self, phy_names, location_names)
            methods = {'lcfs', 'max', 'int', 'lamint'};
            for i=1:length(phy_names)
                for j=1:length(location_names)
                    for k=1:length(methods)
                        disp([uedgerun.print_prefix ' Collecting values for "' phy_names{i} '" on "' upper(location_names{j}) '" using "' methods{k} '" method ...'])
                        self.extract_0d(phy_names{i}, location_names{j}, methods{k});
                    end
                end
            end
        end
        
        function analysis(self, redo)
            %% check arguments
            if nargin < 2
                redo = false;
            end
            
            if ~redo && self.stat_load()
                return
            end
            disp_prefix = [uedgerun.print_prefix ' '];
            %% parsing files
            disp([disp_prefix 'Enumerating UEDGE files ...'])
            self.get_files();
            disp([disp_prefix 'Checking image files ...'])
            self.image_check();
            disp([disp_prefix 'Parsing scanning parameters ...'])
            self.scan_parse();
            %% 2D variables
            phy_names = {'nis(1)', 'tes', 'tis', 'ngs', 'tgs'};
            location_names = {'omp', 'li', 'lo'};
            self.extract_0d_batch(phy_names, location_names);
            % TODO: 'nis', 'ups', 3D or more
            %% 1D variables
            phy_names = {'jl', 'qtl'};
            location_names = {'li'};
            self.extract_0d_batch(phy_names, location_names);
            
            phy_names = {'jr', 'qtr'};
            location_names = {'lo'};
            self.extract_0d_batch(phy_names, location_names);
            %% 0D variables
            phy_names = {'Il', 'Ir'};
            for i=1:length(phy_names)
                self.collect_0d(phy_names{i});
            end
            %% save
            self.stat_save();
        end
        
        function stat_save(self)
            prop_names = fieldnames(self);
            prop_names{end+1} = 'work_dir';
            prop_names{end+1} = 'files';
            res = struct();
            for i=1:length(prop_names)
                p = prop_names{i};
                res.(p) = self.(p);
            end
            save(['stat_' self.work_dir '.mat'], 'res');
        end
        
        function is_loaded = stat_load(self)
            is_loaded = false;
            
            f = ['stat_' self.work_dir '.mat'];
            if ~exist(f, 'file')
                warning(['"' f '" not exist!'])
                return
            end
            
            mat = load(f);
            if ~isfield(mat, 'res')
                warning('not valid uedgestat file!')
                return
            end
            
            prop_names = fieldnames(mat.res);
            for i=1:length(prop_names)
                p = prop_names{i};
                self.(p) = mat.res.(p);
            end
            is_loaded = true;
        end
        
        function slice(self, inds, varargin)
            %% check arguments
            Args.Attributes = [];
            Args = parseArgs(varargin, Args);
            %% set fnames
            fnames = fieldnames(self);
            if ~isempty(Args.Attributes)
                if ischar(Args.Attributes)
                    Args.Attributes = {Args.Attributes};
                end
                fnames = Args.Attributes;
            end
            %% slice
            for i=1:length(fnames)
                fn = fnames{i};
                fn_paths = strsplit(fn, '.');
                %% root attr
                if length(fn_paths) == 1
                    try
                        self.(fn);
                    catch
                        error(['Not a valid attribute of this instance: ' fn])
                    end
                    
                    fnames_sub = fieldnames(self.(fn));
                    attributes_sub = {};
                    for j=1:length(fnames_sub)
                        attributes_sub{end+1} = [fn '.' fnames_sub{j}];
                    end
                    self.slice(inds, 'Attributes', attributes_sub)
                    continue
                end
                %% sub attr
                attr_sub = self;
                for j=1:length(fn_paths)
                    attr_sub = attr_sub.(fn_paths{j});
                end
                
                if isstruct(attr_sub)
                    fnames_sub = fieldnames(attr_sub);
                    attributes_sub = {};
                    for j=1:length(fnames_sub)
                        attributes_sub{end+1} = [fn '.' fnames_sub{j}];
                    end
                    self.slice(inds, 'Attributes', attributes_sub)                    
                    continue
                end
                %% slice
%                 disp([fn ':' num2str(length(attr_sub))])
%                 disp([fn_paths ':' num2str(length(attr_sub))])
                self = setfield(self,fn_paths{:}, attr_sub(inds));
            end
        end
        
        function para_names = filter_remain(self, filter, single_check)
            %% para_names = filter_remain(self, filter, single_check)
            
            %% check arguments
            if nargin < 3
                single_check = false;
            end
            %% get remain sets
            para_names = fieldnames(filter);
            scan_names = self.scan_get_names();
            para_names = setdiff(scan_names, para_names);
            %% check
            if single_check
                assert(length(para_names)==1, 'Input "filter" has not been properly set!')
                para_names = para_names{1};
            end
        end
        
        function inds = filter_select(self, filter)
            %% gen inds
            para_names = fieldnames(filter);
            scan_names = self.scan_get_names();
            inds = true(1, length(self.files));
            for i=1:length(para_names)
                para_name = para_names{i};
                assert(any(contains(scan_names, para_name)), ['invalid filter: "' para_name '"!'])
                para_val = filter.(para_name);
                assert(~isempty(find(self.scan_get_space(para_name)==para_val, 1)), [para_name '=' num2str(para_val) ', not in scan space!'])
                inds = inds & self.scan_get_values(para_name) == para_val;
            end            
            %% sort if only one left
            left_names = self.filter_remain(filter);
            if length(left_names) ~= 1
                return
            end
            scan_name = left_names{1};
            scan_vals = self.scan_get_values(scan_name);
            inds_tmp = 1:length(scan_vals);
            inds = inds_tmp(inds);
            
            [~, inds_tmp] = sort(scan_vals(inds));
            inds = inds(inds_tmp);
        end
                
        function str = filter_gen_str(~, filter)
            filter_para_names = fieldnames(filter);
            str = {};
            for i=1:length(filter_para_names)
                fpn = filter_para_names{i};
                str{i} = [fpn '=' num2str(filter.(fpn))];
            end
        end
        
        function fig = plot_profiles(self, phy_name, poloidal_location, filter, varargin)
            %% fig = plot_profiles(self, phy_name, poloidal_location, filter, 'UseSubplot', 0, 'FontSize', 25)
            
            Args.UseSubplot = false;
            Args.FontSize = 25;
            Args.RemoveGhost = false;
            Args = parseArgs(varargin, Args, {'UseSubplot', 'RemoveGhost'});
            %% get data
            scan_para_name = self.filter_remain(filter);
            if length(scan_para_name) < 1
                fnames = fieldnames(filter);
                filter_multiple = cellfun(@(name) length(filter.(name)) > 1, fnames);
                assert(sum(filter_multiple) == 1, 'Input filter should have only one parameter that has multiple values!')
                scan_para_name = char(fnames(filter_multiple));
                scan_para_vals = self.scan_get_values(scan_para_name);
                para_vals = filter.(scan_para_name);
                filter = rmfield(filter, scan_para_name);
                inds = self.filter_select(filter);
                inds_bit_sel = false(size(inds));
                diff_val = min(diff(sort(scan_para_vals(inds))))/2;
                for i=1:length(inds)
                    ind_tmp = findvaluecrit(para_vals, scan_para_vals(inds(i)), diff_val);
                    inds_bit_sel(i) = ~isnan(ind_tmp);
                end
                inds = inds(inds_bit_sel);
            else
                assert(length(scan_para_name) == 1, ['Input filter should have ' num2str(length(self.scan_get_names())-1) ' fields!'])
                scan_para_name = scan_para_name{1};
                scan_para_vals = self.scan_get_values(scan_para_name);
                inds = self.filter_select(filter);
            end
            
            para_vals = scan_para_vals(inds);
            file_list = self.files(inds);
            %% plot using subplot
            fig = uedgedata.figure;
            if Args.UseSubplot
                row_no = 4;
                col_no = ceil(length(file_list)/4);
                for i=1:length(file_list)
                    subplot(row_no, col_no, i)
                    f = file_list(i);
                    file_path = abspath(f);
                    ud = uedgedata(file_path);
                    ud.plot_1d(phy_name, 'poloidallocationindex', poloidal_location, 'map2omp');
                    title([scan_para_name '=' num2str(para_vals(i))])

                    if i==1
                        s = self.filter_gen_str(filter);
                        title(strjoin(s, ','))
                    end
                end
                set(gca, 'fontsize', Args.FontSize)
                return
            end
            %% plot in single figure
            legend_str = {};
            for i=1:length(file_list)
                ud = uedgedata(abspath(file_list(i)));
                if haselement({'jr','jl','qtr','qtl'},lower(phy_name))
                    ud.plot_1d_target(phy_name, 'map2omp', 'RemoveGhost', Args.RemoveGhost);
                else
                    ud.plot_1d_profile(phy_name, poloidal_location, 'map2omp', 'RemoveGhost', Args.RemoveGhost);
                end
                hold on
                legend_str{end+1} = [scan_para_name '=' num2str(para_vals(i))];
            end
            legend(legend_str)
            s = self.filter_gen_str(filter);
            title(strjoin(s, ','))
            set(gca, 'fontsize', Args.FontSize)
        end
        
        function fig = plot_0d(self, phy_uri_x, phy_uri_y, filter, varargin)
            if nargin < 4
                filter = struct();
            end
            
            Args.FontSize = 25;
            Args.PlotFilter = 1;
            Args.RemainParaIndex = 1;
            Args = parseArgs(varargin, Args, {'PlotFilter'}); 
            
            scan_para_name = self.filter_remain(filter);
            %% recursive call
            if length(scan_para_name) > 1
                assert(length(scan_para_name) == 2, '"filter" should leave at most two scans!')
                scan_name = scan_para_name{Args.RemainParaIndex};
                scan_vals = self.scan_get_space(scan_name);
                fig = uedgedata.figure;
                hold on
                legend_str = {};
                for i=1:length(scan_vals)
                    filter_tmp = filter;
                    filter_tmp.(scan_name) = scan_vals(i);
                    self.plot_0d(phy_uri_x, phy_uri_y, filter_tmp, varargin{:});
                    legend_str{end+1} = [scan_name '=' num2str(scan_vals(i))];
                end  
                legend(legend_str)
                return
            end
            %% single call    
            assert(length(scan_para_name) == 1, '"filter" should leave only one scan!')
            x = self.get_struct_value(self.data_0d, phy_uri_x);
            y = self.get_struct_value(self.data_0d, phy_uri_y);
            inds = self.filter_select(filter);
            
            fig = uedgedata.figure;
            plot(x(inds), y(inds), '*:', 'linewidth', 2.5)
            xlabel(phy_uri_x)
            ylabel(phy_uri_y)
            set(gca, 'fontsize', Args.FontSize)
        end
        
        function fig = plot_scan(self, phy_name, poloidal_location, method, filter, varargin)
            %% fig = plot_scan(self, phy_name, poloidal_location, method, filter, 'FontSize', 25, 'PlotFilter', 1)
            
            assert(haselement({'lcfs', 'max', 'int'}, method), '"method" should be in {"lcfs", "max", "int"}')
            Args.FontSize = 25;
            Args.PlotFilter = 1;
            Args = parseArgs(varargin, Args, {'PlotFilter'});
            
            font_size = Args.FontSize;
            %% recursive call check
            fnames = fieldnames(filter);
            flag_recursive = 0;
            for i=1:length(fnames)
                n = fnames{i};
                v = filter.(n);
                if length(v) > 1
                    flag_recursive = 1;
                    break
                end
            end
            %% recursive call
            if flag_recursive
                Args_tmp = Args;
                Args_tmp.PlotFilter = 0;
                varargin_tmp = struct2vararg(Args_tmp);
                filter_tmp = filter;
                legend_str = {};
                for i=1:length(v)
                    filter_tmp.(n) = v(i);
                    self.plot_scan(phy_name, poloidal_location, method, filter_tmp, varargin_tmp{:})
                    hold on
                    legend_str{end+1} = [n '=' num2str(v(i))];
                end
                legend(legend_str)
                filter_tmp = rmfield(filter_tmp, n);
                s = self.filter_gen_str(filter_tmp);
                title(strjoin(s,','))
                return
            end
            %% extract scan para vals
            scan_para_name = self.filter_remain(filter, 1);
            scan_para_vals = self.scan_get_values(scan_para_name);
            %% extract physical vals
            phy_vals = self.extract_0d(phy_name, poloidal_location, method);
            %% plot
            inds = self.filter_select(filter);
            fig = uedgedata.figure;
            plot(scan_para_vals(inds), phy_vals(inds), '*--', 'linewidth', 2.5, 'markersize',8)
            xlabel([scan_para_name '-scan'])
            ylabel([phy_name '_{' method '}'])
            set(gca, 'fontsize', font_size)
            
            if Args.PlotFilter
                s = self.filter_gen_str(filter);
                title(strjoin(s,','))
            end
        end
        
        function animate_2d(self, phy_name, filter, varargin)
            %% animate_2d(self, phy_name, filter, 'PauseTime', 0.1)
            
            Args.PauseTime = 0.1;
            Args = parseArgs(varargin, Args);
            
            scan_name = self.filter_remain(filter, 1);
            scan_vals = self.scan_get_values(scan_name);
            inds = self.filter_select(filter);
            
            assert(length(inds) > 1, 'Not enough slice to animate!')
            for i=1:length(inds)
                index = inds(i);
                ud = self.get_uedge_data(index);
                ud.plot_2d(phy_name);
                title ( [scan_name '=' num2str(scan_vals(index))] )
                pause(Args.PauseTime);
            end
        end
        
        function fig = contour_profiles(self, phy_name, poloidal_location, filter, varargin)
            %% fig = contour_profiles(self, phy_name, poloidal_location, filter, 'UseSubplot', 0, 'FontSize', 25)
            
            fig = figure;
            self.plot_profiles(phy_name, poloidal_location, filter, varargin{:});
            lines = getlines;
            X = zeros(size(lines));
            C = [];
            for i=1:length(lines)
                l = lines(i);
                sp = strsplit(l.DisplayName, '=');
                X(i) = str2double(sp{end});
                C(:,i) = l.YData;
            end
            Y = l.XData;
            
            x_label = sp{1};
            a = findobj(fig, 'type', 'axe'); 
            y_label = get(get(a, 'xlabel'), 'string');
            c_label = get(get(a, 'ylabel'), 'string');
            
            hold off
            contourfjet(X, Y, C);
            hold on
            hline(0)
            xlabel(x_label)
            ylabel(y_label)
            
            c = colorbar;
            c.Label.String = c_label;
        end
        
        function fig = plot_fail(self, num_paras_max)
            %% check arguments
            if nargin < 2
                num_paras_max = 4;
            end
            assert(num_paras_max <= 4, '"num_paras_max" should be less than 5!')
            %% get fail cases
            usc = uedgescan(self.work_dir);
            [all_cases, scan_names] = usc.scan_combiner();
            converged_cases = self.scan_paras_combiner(scan_names);
            fail_cases = setdiff(all_cases, converged_cases, 'rows');
            %% get scan spaces
            scan_spaces = {};
            inds_trim = [];
            for i=1:length(scan_names)
                scan_spaces{i} = self.scan_get_space(scan_names{i});
                if length(scan_spaces{i}) == 1
                    inds_trim(end+1) = i;
                end
            end
            %% trim
            scan_names_trim = scan_names(inds_trim);
            scan_names(inds_trim) = [];
            
            scan_spaces_trim = scan_spaces(inds_trim);
            scan_spaces(inds_trim) = [];
            
            fail_cases(:,inds_trim) = [];
            %% sort
            scan_space_lens = cellfun(@(x) length(x), scan_spaces);
            [~, inds_sort] = sort(scan_space_lens);
            
            scan_names = scan_names(inds_sort);
            scan_spaces = scan_spaces(inds_sort);
            fail_cases = fail_cases(:, inds_sort);
            %% select scan value if length(scan_names) is larger than 4
            inds_sel = true(size(fail_cases, 1), 1);
            scan_names_sel_len = length(scan_names) - num_paras_max;
            for i=1:scan_names_sel_len                
                scan_space = scan_spaces{i};
                val_sel = input(['Select value for "' scan_names{i} '":\n' mat2str(scan_space) '\n']);
                val = scan_space(findvalue(scan_spaces{i}, val_sel));
                
                scan_names_trim{end+1} = scan_names{i};
                scan_spaces_trim{end+1} = val;
                
                inds_sel = inds_sel & fail_cases(:,i) == val;
            end
            
            fail_cases = fail_cases(inds_sel, :);
            
            scan_names(1:scan_names_sel_len) = [];
            scan_spaces(1:scan_names_sel_len) = [];
            fail_cases(:, 1:scan_names_sel_len) = [];
            %% plot
            assert(length(scan_names) >= 2, 'The number of the remain scan parameters is less than 2!')
            marker_size = 8;
            
            fig = uedgedata.figure;
            title_str = {};
            for i=1:length(scan_names_trim)
                title_str{i} = num2str(scan_spaces_trim{i}, [scan_names_trim{i} '=%g']);
            end
            title_str = strjoin(title_str, ',');
            
            X = fail_cases(:, 1);
            Y = fail_cases(:, 2);
            if length(scan_names) == 2
                plot(X, Y, 'o', 'markersize', marker_size, 'markerfacecolor', 'k');
                title(title_str);
                xlabel(scan_names{1});
                ylabel(scan_names{2});
                uedgedata.figure_decoration;
                return
            end
            
            Z = fail_cases(:, 3);
            if length(scan_names) == 3
                scatter(X, Y, marker_size*10, Z, 'filled', 'markeredgecolor', 'b');
                title(title_str);
                xlabel(scan_names{1});
                ylabel(scan_names{2});
                ax = colorbar();
                colormap('hot');
                xlabel(ax, scan_names{3});
                uedgedata.figure_decoration;
                return
            end
            
            C = fail_cases(:, 4);
            scatter3(X, Y, Z, marker_size*5, C, 'filled');
            title(title_str);
            xlabel(scan_names{1});
            ylabel(scan_names{2});
            zlabel(scan_names{3});
            ax = colorbar();
            colormap('hot');
            xlabel(ax, scan_names{4});
            uedgedata.figure_decoration;
        end
    end
end