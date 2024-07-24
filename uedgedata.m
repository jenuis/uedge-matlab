% Author: Xiang LIU@ASIPP
% E-mail: xliu@ipp.ac.cn
% Created: 2023-10-11
% Version: V 0.1.10
classdef uedgedata < handle
    properties(Access=private)
        depdencies = {...
            'repopathctrl.m', ...
            'addpath_mdslab.m', ...
            }
        
        phy_prefix = 'bbb.'
        
        bdry_ind_radial_start
        bdry_ind_radial_end
        bdry_ind_poloidal_start
        bdry_ind_poloidal_end
        
        grid_no_poloidal
        grid_no_radial
    end
    
    properties
        file_savedt
        file_image
        file_mesh = 'mesh.hdf5'
        file_job
        ur % a uedgerun instance
        
        nx
        ny
        rm
        zm
        omp_poloidal_ind
        lcfs_radial_ind
        X
        Y
        
        geometry
        rbdry
        zbdry
        zmid
        yylb
        yyrb
        psinormc
    end
    
    methods(Static)
        function value = permute(value, warning_off)
            %% check arguments
            if nargin < 2
                warning_off = false;
            end
            %% get dim
            dim = dimnum(value, warning_off);
            
            if dim < 2
                return
            end
            %% permute
            value = permute(value, dim:-1:1);
            if ~warning_off && dimnum(value, warning_off) < dim
                warning('Input value is automatically shrinked!')
            end
        end  
        
        function uri = check_uri(uri)
            if any(contains(uri, '.'))
                uri = strrep(uri, '.', '/');
            end
            
            if uri(1) ~= '/'
                uri = ['/' uri];
            end
        end
        
        function flag = is_h5(file)
            try
                h5info(file);
                flag = true;
            catch
                flag = false;
            end
        end
        
        function flag = is_uedge_file(file_savedt)
            try
                code_name = h5readatt(file_savedt, '/bbb', 'code');
                flag = strcmpi(code_name, 'uedge');
                h5info(file_savedt,'/bbb/nis');
            catch
                flag = false;
            end
        end
        
        function flag = is_phy_struct(phy)
            flag = isstruct(phy);
            if flag == false
                return
            end
            
            field_names = {'name', 'data', 'unit'};
            for i=1:length(field_names)
                flag = flag && isfield(phy, field_names{i});
            end
        end
        
        function fig = figure(varargin)
            Args.FigurePosition = '';
            Args = parseArgs(varargin, Args);
            
            current_fig = get(groot, 'currentfigure');
            if isempty(current_fig)
                fig = figure;
            else
                fig = current_fig;
            end
            set(fig, 'color', 'w')
            
            if ~isempty(Args.FigurePosition)
                setfigposition(Args.FigurePosition);
            end
        end
        
        function figure_decoration(varargin)
            Args.FontSize = 25;
            Args = parseArgs(varargin, Args);
            set(gca, 'fontsize', Args.FontSize);
        end
    end
    
    methods(Access=private)
        function check_dependency(self)
            for i=1:length(self.depdencies)
                dep = self.depdencies{i};
                assert(exist(dep)==2, ['Can not find dependency: "' dep '"'])
            end
        end

        function patch_decoration(self, fig)
            if nargin < 2
                fig = gcf;
            end
            
            figure(fig)
            xlim([min(self.rm, [], 'all') max(self.rm, [], 'all')])
            ylim([min(self.zm, [], 'all') max(self.zm, [], 'all')])
            daspect([1 1 1])
            box on
            grid on
            xlabel('R [m]')
            ylabel('Z [m]')
        end
        
        function compare_experiment_profiles(self, rrsep_exp, ne_exp, te_exp, ti_exp, d, chie, chii)
            %% check argument
            if nargin < 6
                d = [];
                chie = [];
                chii = [];
            end
            %% plot density
            legend_str = {'n_{e,exp}', 'n_{e,UEDGE}'};
            subplot(1, 3, 1)
            hold on
            if ~isempty(d)
                yyaxis right
                plot(d.r*1e3, d.value, 'linewidth', 2)
                legend_str{end+1} = 'D';
                yyaxis left
            end
            plot(rrsep_exp*1e3, ne_exp, 'linewidth', 2)
            self.plot_1d_profile('nis(1)', 'omp', 'xaxisunit', 'mm');
            xlabel('R-R_{sep}')
            ylabel([])
            legend(legend_str)
            %% plot Te
            legend_str = {'T_{e,exp}', 'T_{e,UEDGE}'};
            subplot(1, 3, 2)
            hold on
            if ~isempty(chie)
                yyaxis right
                plot(chie.r*1e3, chie.value, 'linewidth', 2)
                legend_str{end+1} = '\chi_e';
                yyaxis left
            end
            plot(rrsep_exp*1e3, te_exp, 'linewidth', 2)
            self.plot_1d_profile('tes', 'omp', 'xaxisunit', 'mm');
            xlabel('R-R_{sep}')
            ylabel([])
            legend(legend_str)
            %% plot Ti
            legend_str = {'T_{i,exp}', 'T_{i,UEDGE}'};
            subplot(1, 3, 3)
            hold on
            if ~isempty(chii)
                yyaxis right
                plot(chii.r*1e3, chii.value, 'linewidth', 2)
                legend_str{end+1} = '\chi_i';
                yyaxis left
            end
            plot(rrsep_exp*1e3, ti_exp, 'linewidth', 2)
            self.plot_1d_profile('tis', 'omp', 'xaxisunit', 'mm');
            xlabel('R-R_{sep}')
            ylabel([])
            legend(legend_str)
        end
    end
    
    methods
        function self = uedgedata(file_savedt, file_image)
            %% load dependencies   
            self.check_dependency();
            addpath_mdslab('divlp');
            %% check argument
            if nargin == 0
                return
            end
            % check input profile         
            assert(exist(file_savedt,'file'), 'Input profile does not exist!')
            assert(self.is_uedge_file(file_savedt), 'Input profile is not an UEDGE file!')
            self.file_savedt = file_savedt;
            % set uedgerun instance
            [folder, name] = fileparts(file_savedt);
            file_input = [strrep(name, uedgerun.file_save_prefix, 'rd_in') '.py'];
            self.ur = uedgerun(...
                {file_input, fullfile(folder, 'rd_in_new.py'), fullfile(folder, 'rd_in.py')}, ...
                file_savedt ...
                );
            % set other propeties
            self.file_mesh = fullfile(folder, [self.ur.file_mesh_prefix, self.ur.file_extension]);
            self.file_job = self.ur.check_existence(strrep(self.file_savedt, self.ur.file_extension, '.mat'));  
            %% set mesh variables
            self.nx = double(self.savedt_read('com.nx'));
            self.ny = double(self.savedt_read('com.ny'));
            self.rm = self.savedt_read('com.rm');
            self.zm = self.savedt_read('com.zm');
            
            self.set_bdry_index();
            self.set_patch_XY();
            self.find_midplane_index();
            %% check image file
            if nargin == 2 && exist(file_image, 'file')
                self.file_image = file_image;
            end
            self.image_check();             
            %% load extra mesh
            if self.mesh_load()
                self.find_lcfs_index();
            end
        end
         
        function is_loaded = mesh_load(self, mesh_file)
            %% check arguments
            if nargin < 2
                mesh_file = self.file_mesh;
            end
            
            is_loaded = false;
            if ~exist(mesh_file, 'file')
                warning('Can not find mesh file, try call "uedgedata.mesh_save()" to generate!')
%                 self.file_mesh = []; % keep to call mesh_save()
                return
            end
            %% load
            var_list = {'geometry', 'rbdry', 'zbdry', 'yylb', 'yyrb', 'zmid', 'psinormc'};
            file_savedt_store = self.file_savedt;
            self.file_savedt = mesh_file;
            for i=1:length(var_list)
                n = var_list{i};
                v = self.savedt_read(['com.' n]);
                self.(n) = v;
            end
            self.file_savedt = file_savedt_store;
            if length(self.geometry) == 1 && iscellstr(self.geometry)
                self.geometry = strtrim(self.geometry{1});
            end
            is_loaded= true;
        end
        
        function mesh_save(self)
            %% gen image script
            script_image_old = self.ur.script_image;
            self.ur.script_image = self.ur.generate_uuid_file('prefix', 'uedgeimage');
            contents = {
                'var_list = [', ...
                '"com.geometry",', ...
                '"com.rbdry", # radial coordinates of last closed flux surface from EFIT', ...
                '"com.zbdry", # vertical coordinates of last closed flux surface from EFIT', ...
                '"com.yylb", #radial dist. on left boundary from septrx.', ...
                '"com.yyrb", #radial dist. on right boundary from septrx.', ...
                '"com.zmid", # vertical height of midplane above bottom of EFIT mesh', ...
                '"com.psinormc", # norm spi ', ...
                ']', ...
                '', ...
                ['hdf5_save("' self.file_mesh '", var_list)']};
            contents = strjoin(contents, '\n');
            self.ur.script_save(self.ur.script_image, contents);
            %% run
            self.ur.file_save = self.file_savedt;
            self.ur.script_run_gen('CheckInputDiffLeft', 'Dt', 1e-9);
            self.ur.run();
            %% clean
            delete(self.ur.script_image);
            self.ur.script_image = script_image_old;
            self.ur.file_save = [];
        end

        function flag = image_check(self)
            %% check self.file_image
            flag = false;
            if ~isempty(self.file_image) && exist(self.file_image, 'file')
                flag = true;
                return
            end
            %% check according to self.file_savedt            
            f_image = strrep(self.file_savedt, self.ur.file_save_prefix, self.ur.file_image_prefix);
            if strcmpi(self.file_savedt, f_image)
                return
            end
            
            if exist(f_image, 'file')
                self.file_image = f_image;
                flag = true;
                return
            end
        end
               
        function set_bdry_index(self)
            %% get all grid size
            self.grid_no_poloidal = size(self.rm, 1);
            self.grid_no_radial = size(self.rm, 2);
            %% get ghost no
            ghost_no_poloidal = (self.grid_no_poloidal-self.nx)/2;
            ghost_no_radial = (self.grid_no_radial-self.ny)/2;
            %% set bdry ind
            self.bdry_ind_radial_start = 1 + ghost_no_radial;
            self.bdry_ind_radial_end = self.ny + ghost_no_radial;
            self.bdry_ind_poloidal_start = 1 + ghost_no_poloidal;
            self.bdry_ind_poloidal_end = self.nx + ghost_no_poloidal;  
        end
        
        function inds = get_ghost_inds(self, direction)
            %% poloidal direction
            if strcmpi(direction, 'x') || ...
                    strcmpi(direction, 'poloidal')
                inds = [1:self.bdry_ind_poloidal_start-1 self.bdry_ind_poloidal_end+1:self.grid_no_poloidal];
                return
            end
            %% radial direction
            if strcmpi(direction, 'y') || ...
                    strcmpi(direction, 'radial')
                inds = [1:self.bdry_ind_radial_start-1 self.bdry_ind_radial_end+1:self.grid_no_radial];
                return
            end
        end
        
        function set_patch_XY(self)
            %% check
            if ~isempty(self.X) && ~isempty(self.Y)
                return
            end
            %% cal and set self.X and self.Y
            self.X = zeros(4, self.grid_no_radial*self.grid_no_poloidal); % 4*N arrays
            self.Y = self.X; % 4*N arrays
            for ix = 1:self.grid_no_poloidal
                for iy = 1:self.grid_no_radial
                    v = zeros(4, 2);
                    v(1, :) = [self.rm(ix, iy, 2), self.zm(ix, iy, 2)];
                    v(2, :) = [self.rm(ix, iy, 3), self.zm(ix, iy, 3)];
                    v(3, :) = [self.rm(ix, iy, 5), self.zm(ix, iy, 5)];
                    v(4, :) = [self.rm(ix, iy, 4), self.zm(ix, iy, 4)];
                    
                    ind = (ix-1)*self.grid_no_radial+iy;
                    self.X(:, ind) = v(:,1);
                    self.Y(:, ind) = v(:,2);
                end
            end
        end
        
        function index = find_midplane_index(self)
            %% check and return
            if ~isempty(self.omp_poloidal_ind)
                index = self.omp_poloidal_ind;
                return
            end
            %% get r_max for all poloidal locations
            r_max = zeros(1, self.nx);
            for i=self.bdry_ind_poloidal_start:self.bdry_ind_poloidal_end
                indices = self.get_radial_indices(i);
                r = self.X(:, indices);
                r_max(i) = max(r, [], 'all');
            end
            [~, index] = max(r_max);
            %% set omp_poloidal_ind
            self.omp_poloidal_ind = index;
        end
                   
        function r_lcfs = get_div_rrsep(self, target_position, varargin)
            %% check arguments
            Args.RemoveGhost = false;
            Args = parseArgs(varargin, Args, {'RemoveGhost'});
            %% choose target r_lcfs accoording to image file
            switch lower(target_position)
                case {'li', 'lb'}
                    r_lcfs = self.yylb;
                case {'lo', 'rb'}
                    r_lcfs = self.yyrb;
                otherwise
                    error('"target_position" should be in {"ui","uo","li","lo","lb","rb"}')
            end
            %% remove ghost
            if ~isempty(r_lcfs) && Args.RemoveGhost
                r_lcfs = r_lcfs(self.bdry_ind_radial_start:self.bdry_ind_radial_end);
            end
        end
         
        function poloidal_index = get_poloidal_index(self, poloidal_location)
            %% already a poloidal index, return
            if isnumeric(poloidal_location)
                 assert(poloidal_location >= 1 && poloidal_location <= self.grid_no_poloidal, '"poloidal_location" is numerical, but out of range!')
                 poloidal_index = poloidal_location;
                return
            end
            %% a string representation
            assert(ischar(poloidal_location), '"poloidal_location" should be a string!')
            switch lower(poloidal_location)
                case {'inner_divertor', 'left_boundary', 'lb', 'li'}
                    poloidal_index = self.bdry_ind_poloidal_start;
                case {'outer_divertor', 'right_boundary', 'rb', 'lo'}
                    poloidal_index = self.bdry_ind_poloidal_end;
                case {'outboard_midplane', 'omp'}
                    poloidal_index = self.omp_poloidal_ind;
                otherwise
                    error('Unkown argument "poloidal_location"!')
            end
        end
        
        function indices = get_radial_indices(self, poloidal_location)
            poloidal_index = self.get_poloidal_index(poloidal_location);
            indices = (poloidal_index-1)*self.grid_no_radial + (1:self.grid_no_radial);
        end
                
        function index = find_lcfs_index(self)
            %% check and return
            if ~isempty(self.lcfs_radial_ind)
                index = self.lcfs_radial_ind;
                return
            end
            %% get r_lcfs
            r_lcfs = self.get_div_rrsep('rb');
            assert(~isempty(r_lcfs), '"mesh" file not loaded!')
            %% find lcfs index
            index = findvalue(r_lcfs, 0);
            self.lcfs_radial_ind = index;
        end
             
        function info = savedt_read_info(self, uri)
            %% check arguments
            if nargin < 2
                uri = [];
            end
            %% read and return
            if isempty(uri)
                info = h5info(self.file_savedt);
                return
            end
            
            info = h5info(self.file_savedt, self.check_uri(uri));
        end
        
        function attr_val = savedt_read_attr(self, var_uri, attr_name)
            var_uri = self.check_uri(var_uri);
            attr_val = h5readatt(self.file_savedt, var_uri, attr_name);
        end
        
        function groups = savedt_list_groups(self, uri)
            %% check arguments
            groups = {};
            
            if nargin < 2
                uri = [];
            end
            %% read info
            info = self.savedt_read_info(uri);
            
            if isempty(info.Groups)
                return
            end
            %% recursive call
            this_groups = {info.Groups(:).Name};
            for i=1:length(this_groups)
                next_groups = self.savedt_list_groups(this_groups{i});
                if isempty(next_groups)
                    groups{end+1} = this_groups{i};
                    continue
                end
                groups(end+1:end+length(next_groups)) = next_groups;
            end
        end
        
        function vars = savedt_list_variables(self)
            %% get groups
            vars = {};
            groups = self.savedt_list_groups();
            if isempty(groups)
                return
            end
            %% get variables
            for i=1:length(groups)
                info = self.savedt_read_info(groups{i});
                datasets = info.Datasets;
                if isempty(datasets)
                    continue
                end
                dataset_names = {info.Datasets(:).Name};
                for j=1:length(dataset_names)
                    vars{end+1} = [groups{i} '/' dataset_names{j}];
                end
            end
        end
        
        function val = savedt_read(self, var_uri, varargin)
            Args.WarningOff = false;
            Args = parseArgs(varargin, Args, {'WarningOff'});
            
            var_uri = self.check_uri(var_uri);
            val = h5read(self.file_savedt, var_uri);
            if isnumeric(val)
                val = self.permute(val, Args.WarningOff);
            end
        end
        
        function phy = savedt_read_physical(self, phy_name, varargin)
            ind = regexp(phy_name, '\(\d+\)', 'once');
            if ~isempty(ind)
                dim = str2double(phy_name(ind+1));
                phy_name = phy_name(1:end-ind+1);
            end
            
            phy_uri = [self.phy_prefix phy_name];
            
            if ~isempty(ind)
                dims = self.savedt_size(phy_name);
                assert(length(dims) > 2, ['"' phy_name '" has no third dimension, please remove the index!'])
                assert(dims(3) >= dim, ['The max index of dim-3 of "' phy_name '" is ' num2str(dims(3)) '!'])
            end

            val = self.savedt_read(phy_uri, varargin{:});    
            if ~isempty(ind)
                val = val(:,:,dim);
            end
            
            phy_unit = self.savedt_read_attr(phy_uri, 'units');
            
            if strcmpi(phy_unit, 'j')
                phy_unit = 'eV';
                val = val/self.savedt_read([self.phy_prefix 'ev']);
            end
            
            if strcmpi(phy_unit, 'm^-3')
                phy_unit = '10^{19} m^{-3}';
                val = val*1e-19;
            end
            
            phy.name = phy_name;
            phy.data = val;
            phy.unit = phy_unit;
               
            if ~isempty(ind)
                phy.dim = dim;
            end
        end
        
        function val = savedt_max(self, phy_name)
            phy = self.savedt_read_physical(phy_name);
            val = max(phy.data, [], 'all');
        end
        
        function dims = savedt_size(self, phy_name)
            dims = [];
            uri = [self.phy_prefix phy_name];
            info = self.savedt_read_info(uri);
            if isempty(info)
                return
            end
            
            dims = info.Dataspace.Size;
            dims = flip(dims); % in line with permute
        end
        
        function info = image_read_info(self)
            assert(self.image_check(), ['image file does not exist for "' self.file_savedt '"'])
            info = h5info(self.file_image);
        end

        function vars = image_list_variables(self)
            info = self.image_read_info();
            vars = {};
            for i=1:length(info.Datasets)
                vars{end+1} = info.Datasets(i).Name;
            end
        end

        function phy = image_read_physical(self, phy_name)
            %% check arguments
            assert(self.image_check(), "image file does not exist!")
            image_variables = self.image_list_variables();
            [flag, ind] = haselement(lower(image_variables), lower(phy_name));
            
            Isat_names = {'Il', 'Ir'};
            [flag2, ind2] = haselement(lower(Isat_names), lower(phy_name));
            %% cal
            if ~flag && flag2
                phy_name = Isat_names{ind2};
                js = self.read_physical(strrep(phy_name, 'I', 'j'));
                phy.name = phy_name;
                phy.data = sum(js.data);
                phy.unit = js.unit;
                return
            end
            %% read
            assert(flag, 'Not valid "phy_name"')
            phy_name = image_variables{ind};
            val = h5read(self.file_image, self.check_uri(phy_name));
            val = self.permute(val);
            unit = '';
            try
                unit = h5readatt(self.file_image, self.check_uri(phy_name), 'units');
            catch
            end

            phy.name = phy_name;
            phy.data = val;
            phy.unit = unit;
        end
        
        function phy = read_physical(self, phy_name, varargin)
            %% savedt variables
            phy_name = lower(phy_name);
            variables = self.savedt_list_variables();
            phy_name_real = phy_name;
            if contains(phy_name_real, '(')
                phy_name_real = phy_name_real(1:end-3);
            end
            flag = haselement(lower(variables), self.check_uri([self.phy_prefix phy_name_real]));
            if flag
                phy = self.savedt_read_physical(phy_name, varargin{:});
                return
            end
            %% image variables
            variables = self.image_list_variables();
            variables{end+1} = 'Il';
            variables{end+1} = 'Ir';
            
            phy_name_org = phy_name;
            if haselement(lower(variables), 'jsatr')
                switch lower(phy_name)
                    case 'jl'
                        phy_name = 'jsatl';
                    case 'jr'
                        phy_name = 'jsatr';
                end
            end
            
            [flag, ind] = haselement(lower(variables), phy_name);
            assert(flag, ['"' phy_name '" is none of a savedt and image variables!']);
            phy_name = variables{ind};
            phy = self.image_read_physical(phy_name);
            
            if ~strcmpi(phy_name_org, phy_name)
                phy.name = phy_name_org;
            end
        end

        function fig = plot_bdry(self, varargin)
            Args.ZMidZero = 0;
            Args = parseArgs(varargin, Args, {'ZMidZero'});
            
            fig = self.figure();
            if isempty(self.rbdry) || isempty(self.zbdry)
                return
            end
            
            r = self.rbdry;
            z = self.zbdry;
            if Args.ZMidZero
                z = z - self.zmid;
            end
            
            plot(r, z, 'r', 'linewidth', 1.5)
        end
        
        function plot_wall(self)
            file_wall = self.ur.check_existence({'wall_position.mat', '../wall_position.mat'}, true);
            wall_r = matread(file_wall, 'wall_r');
            wall_z = matread(file_wall, 'wall_z');
            plot(wall_r, wall_z, 'b-', 'linewidth', 2.5)
            axis equal
            xlim([min(wall_r) max(wall_r)])
            ylim([min(wall_z) max(wall_z)])
        end
        
        function fig = plot_mesh(self, varargin)
            Args.LineColor = 'k';
            Args.PlotBdry = 1;
            Args.ZMidZero = 0;
            Args = parseArgs(varargin, Args, {'ZMidZero', 'PlotBdry'});
            
            self.set_patch_XY();
            
            r = self.X;
            z = self.Y;
            if Args.ZMidZero
                z = z - self.zmid;
            end
            
            fig = self.figure();
            patch(r, z, '', 'FaceColor', 'None', 'EdgeColor', Args.LineColor)
            hold on
            
            if Args.PlotBdry
                self.plot_bdry('ZMidZero', Args.ZMidZero);
            end
            
            self.patch_decoration();
            self.figure_decoration();
            
            try
                self.plot_wall()
            catch
            end
        end
        
        function val = serialize(self, val)
            %% serialize only 2D data for plotting
            val_size = size(val);
            assert(length(val_size) == 2 && val_size(1) >= self.nx && val_size(2) >= self.ny, 'Not a valid 2D value passed in!')
            val = reshape(val', numel(val), 1);% N*1 arrays            
        end 
        
        function fig = plot_2d(self, phy_name)
            phy = self.read_physical(phy_name);
            assert(dimnum(phy.data) == 2, '"phy_name" is not a 2D variable!')
            C = self.serialize(phy.data);
            
            fig = self.figure();
            patch(self.X, self.Y, C, 'EdgeColor', 'none')

            self.patch_decoration();
            colormap(jet);
            cb = colorbar;
            ylabel(cb, [phy.name ' [' phy.unit ']'])
            self.figure_decoration();
            set(gca,'colorscale','log')
        end
       
        function vals = extract_1d(self, phy_name_or_phy, poloidal_location, varargin)
            Args.RemoveGhost = 0;
            Args = parseArgs(varargin, Args, {'RemoveGhost'});
            
            if ischar(phy_name_or_phy)
                phy = self.read_physical(phy_name_or_phy);
            else
                assert(self.is_phy_struct(phy_name_or_phy), 'Not a valid physical variable struct!');
                phy = phy_name_or_phy;
            end
            
            assert(dimnum(phy.data) == 2, 'Not a 2D variable!')
            
            ind = self.get_poloidal_index(poloidal_location);
            vals = phy.data(ind, :);
            
            if Args.RemoveGhost
                vals = vals(self.bdry_ind_radial_start:self.bdry_ind_radial_end);
            end
        end
        
        function [r_lcfs, label] = cal_rrsep(self, poloidal_location, varargin)
            Args.XaxisUnit = 'm';
            Args.RemoveGhost = false;
            Args = parseArgs(varargin, Args, {'RemoveGhost'});
            assert(haselement({'m', 'mm'}, Args.XaxisUnit), '"Args.XaxisUnit" should be in {"m", "mm"}')
            
            indices = self.get_radial_indices(poloidal_location);
            
            x = self.X(:, indices);
            y = self.Y(:, indices);
            
            x = mean(x);
            y = mean(y);
            %TODO: use correct one, test against yylb and yyrb
            
            [~, seg_len] = arclength(x([1 1:end]), y([1 1:end]));
            r_lcfs = cumsum(seg_len);
            assert(~isempty(self.lcfs_radial_ind), '"lcfs_radial_ind" is empty!')
            r_lcfs = r_lcfs - r_lcfs(self.lcfs_radial_ind);
            
            if strcmpi(Args.XaxisUnit, 'mm')
                r_lcfs = r_lcfs * 1e3;
            end
            
            if Args.RemoveGhost
                r_lcfs = r_lcfs(self.bdry_ind_radial_start:self.bdry_ind_radial_end);
            end
            
            poloidal_location = num2str(poloidal_location, '%i');
            label = ['R-R_{LCFS,' upper(poloidal_location) '} [' Args.XaxisUnit ']'];
        end
        
        function [fit_res, fit_data] = cal_lambda(self, phy_name, varargin)
            Args.RemoveGhost = true;
            Args.PlotResult = true;
            Args.ZeroBg = 1;
            Args.MsParallel = 0;
            Args.MsStartNo = 4;
            Args.FitBdry = [1000 100 100 50 100; 0 0 0 -50 -100];
            Args = parseArgs(varargin, Args, {'RemoveGhost', 'PlotResult', 'ZeroBg'});
            
            assert(haselement({'jr','jl','qtr','qtl'},lower(phy_name)), '"phy_name" should be target particle/heat flux!')
            phy = self.read_physical(phy_name);
            
            fit_data.xdata = self.cal_rrsep('OMP', 'RemoveGhost', Args.RemoveGhost);
            fit_data.ydata = abs(phy.data);
            fit_data.xtype = 'rrsep';
            fit_data.ytype = phy_name;
            
            fit_data.xdata = fit_data.xdata(:)*1e3;
            fit_data.ydata = fit_data.ydata(:);
            if Args.RemoveGhost
                fit_data.ydata = fit_data.ydata(self.bdry_ind_radial_start:self.bdry_ind_radial_end);
            end
            
            plot_res = Args.PlotResult;
            Args = rmfield(Args, 'RemoveGhost');
            Args = rmfield(Args, 'PlotResult');
            
            varargin = struct2vararg(Args);
            fit_res = prbfit.eichfit(fit_data, varargin{:});
            
            if plot_res
                prbfit.fitdata_plot(fit_data, fit_res);
            end
        end
        
        function fig = plot_1d_common(self, phy, poloidal_location_x, varargin)
            %% check arguments
            Args.XaxisUnit = 'm';
            Args.RemoveGhost = false;
            Args.PlotLCFS = true;
            Args = parseArgs(varargin, Args, {'RemoveGhost', 'PlotLCFS'});
            assert(self.is_phy_struct(phy), '"phy" is not a valid physical variable struct!')
            assert(dimnum(phy.data) == 1, '"phy.data" should be 1D')
            %% get xaxis
            [x, x_label] = self.cal_rrsep(poloidal_location_x, ...
                'XaxisUnit', Args.XaxisUnit, ...
                'RemoveGhost', Args.RemoveGhost);
            %% get yaxis
            y = phy.data;
            if Args.RemoveGhost
                y = y(self.bdry_ind_radial_start:self.bdry_ind_radial_end);
            end
            %% plot
            fig = self.figure();
            plot(x, y, 'linewidth', 2);
            
            hold on            
            if Args.PlotLCFS
                vline(0)
            end
            
            y_label = phy.name;
            if ~isempty(phy.unit)
                y_label = [y_label ' [' phy.unit ']'];
            end
            
            xlabel(x_label);
            ylabel(y_label);
            self.figure_decoration();
        end
        
        function fig = plot_1d_profile(self, phy_name, poloidal_location, varargin)
            %% fig = plot_1d_profile(self, phy_name, poloidal_location,  'Map2OMP', false, 'XaxisUnit', 'm', 'RemoveGhost', false, 'PlotLCFS', true)
            
            %% check arguments
            Args.Map2OMP = false;
            Args.XaxisUnit = 'm';
            Args.RemoveGhost = false;
            Args.PlotLCFS = true;
            Args = parseArgs(varargin, Args, {'Map2OMP', 'RemoveGhost', 'PlotLCFS'});
            %% get data
            phy = self.read_physical(phy_name);
            phy.data = self.extract_1d(phy, poloidal_location);
            if ischar(poloidal_location)
                phy.name = [phy.name '_{' poloidal_location '}'];
            else
                phy.name = [phy.name '_{' num2str(poloidal_location, 'pol-ind=%i') '}'];
            end
            
            if Args.Map2OMP
                poloidal_location = 'omp';
            end
            %% plot
            Args = rmfield(Args, 'Map2OMP');
            varargin = struct2vararg(Args);
            fig = self.plot_1d_common(phy, poloidal_location, varargin{:});            
        end
        
        function fig = plot_1d_target(self, phy_name, varargin)
            %% fig = plot_1d_target(self, phy_name, 'Map2OMP', false, 'XaxisUnit', 'm', 'RemoveGhost', false, 'PlotLCFS', true)
            
            %% check arguments
            Args.Map2OMP = false;
            Args.XaxisUnit = 'm';
            Args.RemoveGhost = false;
            Args.PlotLCFS = true;
            Args = parseArgs(varargin, Args, {'Map2OMP', 'RemoveGhost', 'PlotLCFS'});
            assert(haselement({'jr','jl','qtr','qtl'},lower(phy_name)), '"phy_name" should be target particle/heat flux!')
            %% get data
            phy = self.read_physical(phy_name);
            
            poloidal_location = [phy_name(end) 'b'];
            if Args.Map2OMP
                poloidal_location = 'omp';
            end
            %% plot
            Args = rmfield(Args, 'Map2OMP');
            varargin = struct2vararg(Args);
            fig = self.plot_1d_common(phy, poloidal_location, varargin{:});
        end
        
        function fig = plot_1d_patch(self, phy_name, poloidal_location, varargin)
            %% fig = plot_1d_patch(self, phy_name, poloidal_location, 'Map2OMP', false, 'XaxisUnit', 'm')
            
            %% check arguments
            Args.Map2OMP = false;
            Args.XaxisUnit = 'm';
            Args = parseArgs(varargin, Args, {'Map2OMP'});
            %% read data
            phy = self.read_physical(phy_name);
            %% plot
            fig = self.figure();
            indices = self.get_radial_indices(poloidal_location);
            val = self.serialize(phy.data); % 2D data
            v = val(indices);
              
            if Args.Map2OMP
                indices = self.get_radial_indices('omp');
            end

            x = self.X(:, indices);
            y = self.Y(:, indices);

            patch(x, y, v, 'EdgeColor', 'none')

            ylabel([phy.name ' [' phy.unit ']'])
        end
        
        function fig = view(self, varargin)
            %% check arguments
            Args.RemoveGhost = true;
            Args.Target = {'lo', 'li'};
            Args.DetachTemperature = 3;
            Args = parseargs(Args, varargin{:});
            %% create figure
            figure;
            fig = self.figure;
            setfigposition;
            
            row_no = 2;
            col_no = 4;
            sub_no = 1;
            %% 1D OMP
            subplot(row_no, col_no, sub_no)
            yyaxis left
            self.plot_1d_profile('nis(1)', 'omp', 'RemoveGhost', Args.RemoveGhost);
            yyaxis right
            self.plot_1d_profile('ngs', 'omp', 'RemoveGhost', Args.RemoveGhost);
            try
                textbp(self.file_savedt, 'interpreter', 'none')
            catch
            end
            
            sub_no = sub_no + 1;
            subplot(row_no, col_no, sub_no)
            self.plot_1d_profile('tes', 'omp', 'RemoveGhost', Args.RemoveGhost);
            self.plot_1d_profile('tis', 'omp', 'RemoveGhost', Args.RemoveGhost);
            ylabel('T_{omp} [eV]')
            legend('T_e', 'T_i')            
            %% 1D target
            sub_no = sub_no + 1;
            subplot(row_no, col_no, sub_no)
            yyaxis left
            self.plot_1d_profile('nis(1)', Args.Target, 'RemoveGhost', Args.RemoveGhost);
            yyaxis right
            self.plot_1d_profile('ngs', Args.Target, 'RemoveGhost', Args.RemoveGhost);
            
            sub_no = sub_no + 1;
            subplot(row_no, col_no, sub_no)
            self.plot_1d_profile('tes', Args.Target, 'RemoveGhost', Args.RemoveGhost);
            self.plot_1d_profile('tis', Args.Target, 'RemoveGhost', Args.RemoveGhost);
            ylabel(['T_{' upper(Args.Target) '} [eV]'])
            legend('T_e', 'T_i')
%             plot(0, Args.DetachTemperature*1.01, 'w')
            l = getlines(); l = l(1);
            if all(l.YData <= Args.DetachTemperature)
                textbp('Detached', 'fontsize', 25, 'color', 'r')
            end
            hline(Args.DetachTemperature);
            vline(0);
            
            tag = strrep(lower(Args.Target(end)), 'o', 'r');
            tag = strrep(tag, 'i', 'l');
            
            sub_no = sub_no + 1;
            subplot(row_no, col_no, sub_no)
            try
                self.cal_lambda(['j' tag], 'RemoveGhost', Args.RemoveGhost);
            catch ME
                textbp(ME.message, 'color', 'r')
            end
            
            sub_no = sub_no + 1;
            subplot(row_no, col_no, sub_no)
            try
                self.cal_lambda(['qt' tag], 'RemoveGhost', Args.RemoveGhost);
            catch ME
                textbp(ME.message, 'color', 'r')
            end
            %% 2D
            sub_no = sub_no + 1;
            subplot(row_no, col_no, sub_no)
            self.plot_2d('nis(1)');
            
            sub_no = sub_no + 1;
            subplot(row_no, col_no, sub_no)
            self.plot_2d('tes');
        end
        
        function fig = plot_experiment(self, experiment_file)
            %% check argument
            if nargin < 2
                experiment_file = 'experiment_profiles.mat';
            end
            %% load data
            rrsep_exp = matread(experiment_file, 'r', true);
            ne_exp = matread(experiment_file, 'n', true);
            te_exp = matread(experiment_file, 'te', true);
            ti_exp = matread(experiment_file, 'ti', true);
            %% plot
            fig = figure;
            self.compare_experiment_profiles(rrsep_exp, ne_exp, te_exp, ti_exp);
            self.figure_decoration;
            setfigposition;
        end
                
        function status = rerun(self, varargin)
            %% check argument
            Args.FileSavedt = self.file_savedt;
            Args = parseArgs(varargin, Args);
            %% check init file
            file_init_max_time = strrep(self.file_savedt, uedgerun.file_extension, ['_max-iter' uedgerun.file_extension]);
            if exist(file_init_max_time, 'file') == 2
                disp([uedgerun.print_prefix ' Restart with "max-iter" file ...'])
                self.ur.file_init = file_init_max_time;
            end
            %% gen script
            self.ur.file_save = Args.FileSavedt;
            self.ur.script_run_gen('CheckInputDiffLeft', 'Dt', 1e-9);
            %% run
            status = self.ur.run();
            %% clean
            self.ur.file_init = self.file_savedt;
            self.ur.file_save = [];
            %% set self.file_image and rename image file
            if status ~= 0 || exist(self.ur.script_image, 'file') ~= 2
                return
            end
            % image file named according to savedt file
            file_image_new = strrep(self.file_savedt, uedgerun.file_save_prefix, uedgerun.file_image_prefix);
            if exist(file_image_new, 'file')
                self.file_image = file_image_new;
                return
            end
            % image file with default name
            case_dir = fileparts(self.file_savedt);
            file_image_new = fullfile(case_dir, [uedgerun.file_image_prefix, uedgerun.file_extension]);
            assert(exist(file_image_new, 'file'), 'Image file not generated!')
            if isempty(self.file_image)
                self.file_image = file_image_new;
                return
            end
            % image file has already been specified
            if ~strcmpi(file_image_new, self.file_image)
                movefile(file_image_new, self.file_image)
            end
        end
        
        function status = retry(self)
            %% set rundt parameter
            max_time_org = self.ur.rundt_max_time;
            self.ur.rundt_max_time = 1e-7;
            %% gen script
            file_retry = strrep(self.file_savedt, uedgerun.file_extension, ['_retry' uedgerun.file_extension]);
            self.ur.file_save = file_retry;
            self.ur.script_run_gen('BCIncludeDt', 'Dt', 1e-10);
            %% run
            status = self.ur.run();
            assert(status == 103, 'Should return "max time reached!"')
            %% cal rerun
            self.ur.rundt_max_time = max_time_org;
            self.ur.file_init = strrep(file_retry, uedgerun.file_extension, ['_max-time' uedgerun.file_extension]);
            self.rerun();
            %% clean
            self.ur.file_init = self.file_savedt;
            self.ur.file_save = [];
        end
        
        function tune_transport_coeff(self, varargin)
            %% check arguments
            Args.ExperimentProfile = '';
            Args.TransportCoeff = [];
            Args.TransportCoeffSavePath = 'transp_coeffs.mat';
            Args.OverWrite = false;
            Args = parseArgs(varargin, Args, {'OverWrite'});
            
            experiment_file = Args.ExperimentProfile;
            transp_coeff = Args.TransportCoeff;
            
            assert(~isempty(transp_coeff), '"TransportCoeff" is empty!')
            assert(isstruct(transp_coeff), '"TransportCoeff" should be a struct!')
            
            field_names = {'d', 'chie', 'chii'};
            for i=1:length(field_names)
                n = field_names{i};
                assert(fieldexist(transp_coeff, n), ['"TransportCoeff" has no attribute: "' n '"!'])
                v = transp_coeff.(n);
                assert(fieldexist(v, 'r'), ['"TransportCoeff" has no attribute: "' n '.r"!'])
                assert(fieldexist(v, 'value'), ['"TransportCoeff" has no attribute: "' n '.value"!'])
                assert(length(v.r) == length(v.value), ['"' n '.r" has different length from "' n '.value"!'])
            end
            %% load experiment profiles
            rrsep_exp = NaN;
            ne_exp = NaN;
            te_exp = NaN;
            ti_exp = NaN;
            if ~isempty(experiment_file) && exist(experiment_file, 'file') == 2
                rrsep_exp = matread(experiment_file, 'r', true);
                assert(sum(rrsep_exp > 0) > 1 && sum(rrsep_exp < 0) > 1, '"r" should be "R - Rsep"!')
                ne_exp = matread(experiment_file, 'n', true);
                te_exp = matread(experiment_file, 'te', true);
                ti_exp = matread(experiment_file, 'ti', true);
            end
            %% collect variables
            gridr  = self.cal_rrsep('omp');
            d      = transp_coeff.d;
            chie   = transp_coeff.chie;
            chii   = transp_coeff.chii;
            %% plot modified transport coefficients
            fh = figure;
            self.compare_experiment_profiles(rrsep_exp, ne_exp, te_exp, ti_exp, d, chie, chii);
            setfigposition;
            [~,~,b]=ginput(2);
            assert(sum(b) == 2, 'User Quit!')
            close(fh)
            %% save transport coefficients
            save(Args.TransportCoeffSavePath, 'gridr', '-v7')
            mf = matfile(Args.TransportCoeffSavePath, 'writable', true);
            mf.d = d;
            mf.chie = chie;
            mf.chii = chii;
            %% modify input file
            uedgerun.input_exist_parameter(self.ur.script_input, 'transp_coeff_path', true);
            contents = uedgerun.input_modify(self.ur.script_input, 'transp_coeff_path', ['"' fullfile(pwd, Args.TransportCoeffSavePath) '"']);
            uedgerun.script_save(self.ur.script_input, contents);
            %% run
            file_savedt_new = strrep(self.file_savedt, uedgerun.file_extension, ['_new' uedgerun.file_extension]);
            status = self.rerun('FileSavedt', file_savedt_new);
            assert(status == 0, 'Failed to converge for the input "TransportCoeff"!')
            %% plot
            ud = uedgedata(file_savedt_new);
            figure;
            setfigposition;
            ud.compare_experiment_profiles(rrsep_exp, ne_exp, te_exp, ti_exp, d, chie, chii);
            ud.view();
            %% overwrite?
            if Args.OverWrite
                movefile(ud.file_savedt, self.file_savedt);
                movefile(ud.file_image, self.file_image);
            end
        end
    end
end