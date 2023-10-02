function [eloc, labels, theta, radius, indices] = readlocs( filename, varargin ); 
 
 if nargin < 1
     help readlocs;
     return;
 end;

 chanformat(1).type         = 'polhemus';
 chanformat(1).typestring   = 'Polhemus native .elp file';
 chanformat(1).description  = [ 'Polhemus native coordinate file containing scanned electrode positions. ' ...
                                'User must select the direction ' ...
                                'for the nose after importing the data file.' ];
 chanformat(1).importformat = 'readelp() function';
 chanformat(2).type         = 'besa';
 chanformat(2).typestring   = 'BESA spherical .elp file';
 chanformat(2).description  = [ 'BESA spherical coordinate file. Note that BESA spherical coordinates ' ...
                                'are different from Matlab spherical coordinates' ];
 chanformat(2).skipline     = -1;
 chanformat(2).importformat = { 'type' 'labels' 'sph_theta_besa' 'sph_phi_besa' 'sph_radius' };
 chanformat(3).type         = 'xyz';
 chanformat(3).typestring   = 'Matlab .xyz file';
 chanformat(3).description  = [ 'Standard 3-D cartesian coordinate files with electrode labels in ' ...
                                'the first column and X, Y, and Z coordinates in columns 2, 3, and 4' ];
 chanformat(3).importformat = { 'channum' '-Y' 'X' 'Z' 'labels'};
 chanformat(4).type         = 'sfp';
 chanformat(4).typestring   = 'BESA or EGI 3-D cartesian .sfp file';
 chanformat(4).description  = [ 'Standard BESA 3-D cartesian coordinate files with electrode labels in ' ...
                                'the first column and X, Y, and Z coordinates in columns 2, 3, and 4.' ...
                                'Coordinates are re-oriented to fit the EEGLAB standard of having the ' ...
                                'nose along the +X axis.' ];
 chanformat(4).importformat = { 'labels' '-Y' 'X' 'Z' };
 chanformat(4).skipline     = -1;
 chanformat(5).type         = 'loc';
 chanformat(5).typestring   = 'EEGLAB polar .loc file';
 chanformat(5).description  = [ 'EEGLAB polar .loc file' ];
 chanformat(5).importformat = { 'channum' 'theta' 'radius' 'labels' };
 chanformat(6).type         = 'sph';
 chanformat(6).typestring   = 'Matlab .sph spherical file';
 chanformat(6).description  = [ 'Standard 3-D spherical coordinate files in Matlab format' ];
 chanformat(6).importformat = { 'channum' 'sph_theta' 'sph_phi' 'labels' };
 chanformat(7).type         = 'asc';
 chanformat(7).typestring   = 'Neuroscan polar .asc file';
 chanformat(7).description  = [ 'Neuroscan polar .asc file, automatically recentered to fit EEGLAB standard' ...
                                'of having ''Cz'' at (0,0).' ];
 chanformat(7).importformat = 'readneurolocs';
 chanformat(8).type         = 'dat';
 chanformat(8).typestring   = 'Neuroscan 3-D .dat file';
 chanformat(8).description  = [ 'Neuroscan 3-D cartesian .dat file. Coordinates are re-oriented to fit ' ...
                                'the EEGLAB standard of having the nose along the +X axis.' ];
 chanformat(8).importformat = 'readneurolocs';
 chanformat(9).type         = 'elc';
 chanformat(9).typestring   = 'ASA .elc 3-D file';
 chanformat(9).description  = [ 'ASA .elc 3-D coordinate file containing scanned electrode positions. ' ...
                                'User must select the direction ' ...
                                'for the nose after importing the data file.' ];
 chanformat(9).importformat = 'readeetraklocs';
 chanformat(10).type         = 'chanedit';
 chanformat(10).typestring   = 'EEGLAB complete 3-D file';
 chanformat(10).description  = [ 'EEGLAB file containing polar, cartesian 3-D, and spherical 3-D ' ...
                                'electrode locations.' ];
 chanformat(10).importformat = { 'channum' 'labels'  'theta' 'radius' 'X' 'Y' 'Z' 'sph_theta' 'sph_phi' ...
                                'sph_radius' };
 chanformat(10).skipline     = 1;
 chanformat(11).type         = 'custom';
 chanformat(11).typestring   = 'Custom file format';
 chanformat(11).description  = 'Custom ASCII file format where user can define content for each file columns.';
 chanformat(11).importformat = '';
 % ---------------------------------------------------------------------------------------------------
 % ----- ADD MORE FORMATS HERE -----------------------------------------------------------------------
 % ---------------------------------------------------------------------------------------------------
 
 listcolformat = { 'labels' 'channum' 'theta' 'radius' 'sph_theta' 'sph_phi' ...
       'sph_radius' 'sph_theta_besa' 'sph_phi_besa' 'gain' 'calib' 'type' ...
       'X' 'Y' 'Z' '-X' '-Y' '-Z' 'custom1' 'custom2' 'custom3' 'custom4' 'ignore' 'not def' };
 
 if isstr(filename) & strcmp(filename, 'getinfos')
    eloc = chanformat;
    labels = listcolformat;
    return;
 end;
 
 g = finputcheck( varargin, ...
    { 'filetype'       'string'  {}                 '';
      'importmode'  'string'  { 'eeglab' 'native' } 'eeglab';
      'defaultelp'  'string'  { 'besa'   'polhemus' } 'polhemus';
      'skiplines'   'integer' [0 Inf]             [];
      'elecind'     'integer' [1 Inf]            [];
      'format'       'cell'     []                    {} }, 'readlocs');
 if isstr(g), error(g); end;  
 
 if isstr(filename)

    if strcmpi(g.filetype, 'autodetect'), g.filetype = ''; end;
    g.filetype = strtok(g.filetype);
    periods = find(filename == '.');
    fileextension = filename(periods(end)+1:end);
    g.filetype = lower(g.filetype);
    if isempty(g.filetype)
        switch lower(fileextension),
         case {'loc' 'locs' }, g.filetype = 'loc';
         case 'xyz', g.filetype = 'xyz'; 
           fprintf( [ 'WARNING: Matlab Cartesian coord. file extension (".xyz") detected.\n' ... 
                   'If importing EGI Cartesian coords, force type "sfp" instead.\n'] );
         case 'sph', g.filetype = 'sph';
         case 'ced', g.filetype = 'chanedit';
         case 'elp', g.filetype = g.defaultelp;
         case 'asc', g.filetype = 'asc';
         case 'dat', g.filetype = 'dat';
         case 'elc', g.filetype = 'elc';
         case 'eps', g.filetype = 'besa';
         case 'sfp', g.filetype = 'sfp';
         otherwise, g.filetype =  ''; 
        end;
        fprintf('readlocs(): ''%s'' format assumed from file extension\n', g.filetype); 
    else 
        if strcmpi(g.filetype, 'locs'),  g.filetype = 'loc'; end
        if strcmpi(g.filetype, 'eloc'),  g.filetype = 'loc'; end
    end;
    if ~isempty(g.filetype) & ~strcmpi(g.filetype, 'custom') ...
            & ~strcmpi(g.filetype, 'asc') & ~strcmpi(g.filetype, 'elc') & ~strcmpi(g.filetype, 'dat')
       indexformat = strmatch(lower(g.filetype), { chanformat.type }, 'exact');
       g.format = chanformat(indexformat).importformat;
       if isempty(g.skiplines)
          g.skiplines = chanformat(indexformat).skipline;
       end;
       if isempty(g.filetype) 
          error( ['readlocs() error: The filetype cannot be detected from the \n' ...
                  '                  file extension, and custom format not specified']);
       end;
    end;
    
    if strcmp(g.filetype, 'asc') | strcmp(g.filetype, 'dat')
        eloc = readneurolocs( filename );
        eloc = rmfield(eloc, 'sph_theta'); % for the conversion below
        eloc = rmfield(eloc, 'sph_theta_besa'); % for the conversion below
    elseif strcmp(g.filetype, 'elc')
        eloc = readeetraklocs( filename );
        %eloc = read_asa_elc( filename ); % from fieldtrip
        %eloc = struct('labels', eloc.label, 'X', mattocell(eloc.pnt(:,1)'), 'Y', ...
        %                        mattocell(eloc.pnt(:,2)'), 'Z', mattocell(eloc.pnt(:,3)'));
        eloc = convertlocs(eloc, 'cart2all');
        eloc = rmfield(eloc, 'sph_theta'); % for the conversion below
        eloc = rmfield(eloc, 'sph_theta_besa'); % for the conversion below
    elseif strcmp(lower(g.filetype(1:end-1)), 'polhemus') | ...
            strcmp(g.filetype, 'polhemus')
        try, 
            [eloc labels X Y Z]= readelp( filename );
            if strcmp(g.filetype, 'polhemusy')
                tmp = X; X = Y; Y = tmp;
            end;
            for index = 1:length( eloc )
                eloc(index).X = X(index);
                eloc(index).Y = Y(index);    
                eloc(index).Z = Z(index);    
            end;
        catch, 
            disp('readlocs(): Could not read Polhemus coords. Trying to read BESA .elp file.');
            [eloc, labels, theta, radius, indices] = readlocs( filename, 'defaultelp', 'besa', varargin{:} );
        end;
    else      

        if isempty(g.skiplines), g.skiplines = 0; end;
        array = load_file_or_array( filename, max(g.skiplines,0));
        if size(array,2) < length(g.format)
            fprintf(['readlocs() warning: Fewer columns in the input than expected.\n' ...
                     '                    See >> help readlocs\n']);
        elseif size(array,2) > length(g.format)
            fprintf(['readlocs() warning: More columns in the input than expected.\n' ...
                     '                    See >> help readlocs\n']);
        end;
        if g.skiplines == -1
            if isempty(array{1,2})
                disp('BESA header detected, skipping three lines...');
                array = load_file_or_array( filename, -2);
            end;
        end;

        indexbeg = 1;
        while isempty(array{indexbeg,1}) | ...
                (isstr(array{indexbeg,1}) & array{indexbeg,1}(1) == '%' )
            indexbeg = indexbeg+1;
        end;
        array = array(indexbeg:end,:);

        for indexcol = 1:min(size(array,2), length(g.format))
            [str mult] = checkformat(g.format{indexcol});
            for indexrow = 1:size( array, 1)
                if mult ~= 1
                    eval ( [ 'eloc(indexrow).'  str '= -array{indexrow, indexcol};' ]);
                else
                    eval ( [ 'eloc(indexrow).'  str '= array{indexrow, indexcol};' ]);
                end;
            end;
        end;
    end;
    if isfield(eloc, 'sph_theta_besa')
        if isnumeric(eloc(1).type)
            disp('BESA format detected ( Theta | Phi )');
            for index = 1:length(eloc)
                eloc(index).sph_phi_besa   = eloc(index).labels;
                eloc(index).sph_theta_besa = eloc(index).type;
                eloc(index).labels         = '';
                eloc(index).type           = '';
            end;
            eloc = rmfield(eloc, 'labels');
        elseif isnumeric(eloc(1).labels)
            disp('BESA format detected ( Elec | Theta | Phi )');
            for index = 1:length(eloc)
                eloc(index).sph_phi_besa   = eloc(index).sph_theta_besa;
                eloc(index).sph_theta_besa = eloc(index).labels;
                eloc(index).labels         = eloc(index).type;
                eloc(index).type           = '';
                eloc(index).radius         = 1;
            end;           
        else
            disp('BESA format detected ( Type | Elec | Theta | Phi | Radius )');           
        end;
        
        try
            eloc = convertlocs(eloc, 'sphbesa2all');
            eloc = convertlocs(eloc, 'topo2all'); % problem with some EGI files (not BESA files)
        catch, disp('Warning: coordinate conversion failed'); end;
        fprintf('Readlocs: BESA spherical coords. converted, now deleting BESA fields\n');   
        fprintf('          to avoid confusion (these fields can be exported, though)\n');   
        eloc = rmfield(eloc, 'sph_phi_besa');
        eloc = rmfield(eloc, 'sph_theta_besa');

    elseif isfield(eloc, 'sph_theta')
        try
            eloc = convertlocs(eloc, 'sph2all');  
        catch, disp('Warning: coordinate conversion failed'); end;
    elseif isfield(eloc, 'X')
        try
            eloc = convertlocs(eloc, 'cart2all');  
        catch, disp('Warning: coordinate conversion failed'); end;
    else 
        try
            eloc = convertlocs(eloc, 'topo2all');  
        catch, disp('Warning: coordinate conversion failed'); end;
    end;
    

    if ~isfield(eloc, 'labels')
        fprintf('readlocs(): Inserting electrode labels automatically.\n');
        for index = 1:length(eloc)
            eloc(index).labels = [ 'E' int2str(index) ];
        end;
    else 
        for index = 1:length(eloc)
            if isstr(eloc(index).labels)
                tmpdots = find( eloc(index).labels == '.' );
                eloc(index).labels(tmpdots) = [];
            end;
        end;
    end;
    if isfield(eloc, 'channum')
        if ~isnumeric(eloc(1).channum)
            error('Channel numbers must be numeric');
        end;
        allchannum = [ eloc.channum ];
        if any( sort(allchannum) ~= allchannum )
            fprintf('readlocs(): Re-sorting channel numbers based on ''channum'' column indices\n');
            [tmp newindices] = sort(allchannum);
            eloc = eloc(newindices);
        end;
        eloc = rmfield(eloc, 'channum');      
    end;
 else
     if isstruct(filename)
         eloc = filename;
     else
         disp('readlocs(): input variable must be a string or a structure');
     end;        
 end;
 if ~isempty(g.elecind)
     eloc = eloc(g.elecind);
 end;
 if nargout > 2
     tmptheta          = { eloc.theta }; % check which channels have (polar) coordinates set
     indices           = find(~cellfun('isempty', tmptheta));
     indbad            = find(cellfun('isempty', tmptheta));
     tmptheta(indbad)  = { NaN };
     theta             = [ tmptheta{:} ];
 end;
 if nargout > 3
     tmprad            = { eloc.radius };
     tmprad(indbad)    = { NaN };
     radius            = [ tmprad{:} ];
 end;
 %tmpnum = find(~cellfun('isclass', { eloc.labels }, 'char'));
 %disp('Converting channel labels to string');
 for index = 1:length(eloc)
     if ~isstr(eloc(index).labels)
         eloc(index).labels = int2str(eloc(index).labels);
     end;
 end;
 labels = { eloc.labels };
 if isfield(eloc, 'ignore')
     eloc = rmfield(eloc, 'ignore');
 end;
 fidnames = { 'nz' 'lpa' 'rpa' };
 for index = 1:length(fidnames)
     ind = strmatch(fidnames{index}, lower(labels), 'exact');
     if ~isempty(ind), eloc(ind).type = 'FID'; end;
 end;
 
 return;
 function array = load_file_or_array( varname, skiplines );
      if isempty(skiplines),
        skiplines = 0;
     end;
     if exist( varname ) == 2
         array = loadtxt(varname,'verbose','off','skipline',skiplines);
     else % variable in the global workspace
          % --------------------------
          try, array = evalin('base', varname);
          catch, error('readlocs(): cannot find the named file or variable, check syntax');
          end;
     end;     
 return;
 function [str, mult] = checkformat(str)
     mult = 1;
     if strcmpi(str, 'labels'),         str = lower(str); return; end;
     if strcmpi(str, 'channum'),        str = lower(str); return; end;
     if strcmpi(str, 'theta'),          str = lower(str); return; end;
     if strcmpi(str, 'radius'),         str = lower(str); return; end;
     if strcmpi(str, 'ignore'),         str = lower(str); return; end;
     if strcmpi(str, 'sph_theta'),      str = lower(str); return; end;
     if strcmpi(str, 'sph_phi'),        str = lower(str); return; end;
     if strcmpi(str, 'sph_radius'),     str = lower(str); return; end;
     if strcmpi(str, 'sph_theta_besa'), str = lower(str); return; end;
     if strcmpi(str, 'sph_phi_besa'),   str = lower(str); return; end;
     if strcmpi(str, 'gain'),           str = lower(str); return; end;
     if strcmpi(str, 'calib'),          str = lower(str); return; end;
     if strcmpi(str, 'type') ,          str = lower(str); return; end;
     if strcmpi(str, 'X'),              str = upper(str); return; end;
     if strcmpi(str, 'Y'),              str = upper(str); return; end;
     if strcmpi(str, 'Z'),              str = upper(str); return; end;
     if strcmpi(str, '-X'),             str = upper(str(2:end)); mult = -1; return; end;
     if strcmpi(str, '-Y'),             str = upper(str(2:end)); mult = -1; return; end;
     if strcmpi(str, '-Z'),             str = upper(str(2:end)); mult = -1; return; end;
     if strcmpi(str, 'custom1'), return; end;
     if strcmpi(str, 'custom2'), return; end;
     if strcmpi(str, 'custom3'), return; end;
     if strcmpi(str, 'custom4'), return; end;
     error(['readlocs(): undefined field ''' str '''']);
