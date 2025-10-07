%% Master script — Macrocirculation build/run/post → Micro build → CSV handoff → Launch Micro & NGVU
% Run from repo root. 

%% -------- CONFIG --------
% Macrocirculation (build & run)
DO_BUILD      = true;      % true: configure & compile Macrocirculation; false: just run (if RUN_MACRO)
CLEAN_BUILD   = false;     % true: delete Macrocirculation/build before build
RUN_MACRO     = true;      % true: run solver and then python postprocess; false: skip run/postprocess
USE_MPI       = false;     % true: run via mpirun
MPI_PROCS     = 1;         % number of ranks if USE_MPI = true

T_END         = 20;                       % --t-end for C++ and Python
OUTPUT_DIR    = './output';               % relative to executable dir
MESH_FILENAME = '37-vessels.json';        % Macrocirculation/data/1d-meshes
SET_INLET     = false;                    % pass inlet explicitly?
INLET_NAME    = 'cw_in';                  % value for --inlet-name if SET_INLET

EXE_NAME      = 'MacrocirculationNonlinear1DSolver';
PROJ_DIR      = 'Macrocirculation';
BUILD_DIRNAME = 'build';
BIN_SUBPATH   = fullfile('bin','macrocirculation');

% Python post-processing options
PY_SCRIPT     = 'write_flow.py';
PY_VESSELS    = '36';                     % passed to --vessels
PY_FILEPATH   = fullfile(OUTPUT_DIR, 'abstract_vessels.json');
PY_OUTPUT     = fullfile(OUTPUT_DIR, 'pica_flow_heart_period.csv');

% Microcirculation build (DUNE)
DO_BUILD_MICRO = true;  % if true: verify compilers==9, clean artifacts, run ./build.sh
MICRO_DIR = fullfile('Microcirculation','Code_Upload','DUNE'); % relative to repo root
MICRO_DUNE_MODULES_TO_CLEAN = { ...
  'dune-common','dune-foamgrid','dune-functions','dune-geometry', ...
  'dune-grid','dune-istl','dune-localfunctions','dune-typetree'};

% Microcirculation run control
RUN_MICRO_IN_BACKGROUND = true;   % true: run as batch (&), false: run blocking in foreground
MICRO_ARG = 'ArtificialNetwork';  % required argument for dune-angiogenesis
MICRO_STDOUT_FILE = 'micro_stdout.log';   % where to redirect micro stdout/stderr

% NVU vessel CSV drop directory (where NGVU writes NVU_Vessel_*.csv, and Micro reads them)
NVU_DROP_DIRNAME = 'NVU_Vessels';  % created inside dune-angiogenesis/build-cmake/src

% NGVU driver script
NGVU_DIR     = 'NGVU';                         % folder at repo root
NGVU_MAIN    = 'quadripartite_NVC_main.m';     % main driver file inside NGVU

%% -------- PATHS --------
thisFile = mfilename('fullpath');
if isempty(thisFile), repo_root = pwd; else, repo_root = fileparts(thisFile); end

proj_dir  = fullfile(repo_root, PROJ_DIR);
build_dir = fullfile(proj_dir, BUILD_DIRNAME);
exe_dir   = fullfile(build_dir, BIN_SUBPATH);
exe_path  = fullfile(exe_dir, EXE_NAME);

mesh_abs  = fullfile(proj_dir, 'data', '1d-meshes', MESH_FILENAME);
out_abs   = fullfile(exe_dir, OUTPUT_DIR);
macro_csv_expected = fullfile(exe_dir, PY_OUTPUT); % expected CSV path

micro_root     = fullfile(repo_root, MICRO_DIR);
angio_build    = fullfile(micro_root, 'dune-angiogenesis', 'build-cmake');
micro_src_dir  = fullfile(angio_build, 'src');
micro_exe_path = fullfile(micro_src_dir, 'dune-angiogenesis');
micro_csv_sink_dir = micro_src_dir;  % destination for pica_flow_heart_period.csv
micro_csv_sink     = fullfile(micro_csv_sink_dir, 'pica_flow_heart_period.csv');

% NVU vessel drop dir (Micro will read; NGVU will write here)
nvu_drop_dir  = fullfile(micro_src_dir, NVU_DROP_DIRNAME);

% NGVU driver path
ngvu_main_path = fullfile(repo_root, NGVU_DIR, NGVU_MAIN);

%% -------- CHECKS --------
if ~isfolder(proj_dir), error('Folder not found: %s', proj_dir); end
if ~isfile(mesh_abs),   error('Mesh file not found: %s', mesh_abs); end

%% -------- BUILD (optional, Macrocirculation) --------
if DO_BUILD
    if CLEAN_BUILD && isfolder(build_dir)
        fprintf('[CLEAN] %s\n', build_dir);
        rmdir(build_dir, 's');
    end
    if ~isfolder(build_dir), mkdir(build_dir); end

    prev = pwd; c0 = onCleanup(@() cd(prev)); 
    cd(build_dir);

    fprintf('[BUILD] cmake ..\n');
    [s1,o1] = system('cmake ..'); fprintf('%s', o1);
    if s1 ~= 0, error('CMake configuration failed.'); end

    fprintf('[BUILD] make -j\n');
    [s2,o2] = system('make -j'); fprintf('%s', o2);
    if s2 ~= 0, error('Compilation failed.'); end
else
    if RUN_MACRO && ~isfile(exe_path)
        error('Executable not found:\n  %s\nSet DO_BUILD = true to compile first.', exe_path);
    end
end

%% -------- RUN SOLVER + PY POST (Macrocirculation) --------
if RUN_MACRO
    if ~isfolder(exe_dir), error('Executable directory not found: %s', exe_dir); end
    if ~isfolder(out_abs), mkdir(out_abs); end

    args = sprintf('--mesh-file "%s" --t-end %g --output-directory "%s"', ...
                   mesh_abs, T_END, OUTPUT_DIR);
    if SET_INLET && ~isempty(INLET_NAME)
        args = sprintf('%s --inlet-name %s', args, INLET_NAME);
    end

    if USE_MPI
        if MPI_PROCS < 1, MPI_PROCS = 1; end
        cmd = sprintf('mpirun -np %d "%s" %s', MPI_PROCS, exe_path, args);
    else
        cmd = sprintf('"%s" %s', exe_path, args);
    end

    fprintf('[RUN] cd %s\n', exe_dir);
    fprintf('[RUN] %s\n', cmd);
    prev = pwd; c1 = onCleanup(@() cd(prev)); 
    cd(exe_dir);

    [st_run,out_run] = system(cmd); fprintf('%s', out_run);
    if st_run ~= 0, error('Solver exited with status %d.', st_run); end
    fprintf('[DONE] Solver finished. Output dir: %s\n', out_abs);

    % Python post-process
    py_path = fullfile(exe_dir, PY_SCRIPT);
    if ~isfile(py_path), error('Python post-processing script not found: %s', py_path); end

    py_cmd = sprintf('python3 "%s" --vessels %s --filepath "%s" --output "%s" --t-end %g', ...
                     py_path, PY_VESSELS, PY_FILEPATH, PY_OUTPUT, T_END);

    fprintf('[POST] %s\n', py_cmd);
    [st_py, out_py] = system(py_cmd); fprintf('%s', out_py);
    if st_py ~= 0, error('Post-processing script exited with status %d.', st_py); end

    fprintf('[DONE] Post-processing complete. CSV: %s\n', macro_csv_expected);
else
    fprintf('[SKIP] RUN_MACRO = false → skipping solver and python postprocess.\n');
end

%% -------- MICRO CIRCULATION BUILD (if enabled) --------
if DO_BUILD_MICRO
    % 1) Verify compilers are version 9
    fprintf('\n[CHECK] Verifying gcc/g++/gfortran are version 9...\n');
    compilers = {'gcc','g++','gfortran'};
    for k = 1:numel(compilers)
        cmdv = sprintf('%s -dumpfullversion -dumpversion', compilers{k});
        [sV, outV] = system(cmdv);
        if sV ~= 0 || isempty(strtrim(outV))
            [~, outV] = system([compilers{k} ' --version']);
        end
        tok = regexp(outV, '(\d+)(\.\d+)*', 'tokens', 'once');
        if isempty(tok)
            error('Could not parse %s version from: %s', compilers{k}, strtrim(outV));
        end
        major = str2double(tok{1});
        if isnan(major) || major ~= 9
            fprintf(2, ['[ERROR] %s major version is not 9 (detected: %s).\n' ...
                        'Install and select v9, e.g.:\n' ...
                        '  sudo apt-get update && sudo apt-get install gcc-9 g++-9 gfortran-9\n' ...
                        '  sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 90\n' ...
                        '  sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-9 90\n' ...
                        '  sudo update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-9 90\n' ...
                        'Then select them:\n' ...
                        '  sudo update-alternatives --config gcc\n' ...
                        '  sudo update-alternatives --config g++\n' ...
                        '  sudo update-alternatives --config gfortran\n'], ...
                        upper(compilers{k}), strtrim(outV));
            error('Please switch %s to version 9 and rerun.', compilers{k});
        else
            fprintf('[OK] %s version: %s', compilers{k}, strtrim(outV));
            if outV(end) ~= newline, fprintf('\n'); end
        end
    end

    % 2) Clean specified build artifacts
    if ~isfolder(micro_root)
        error('Microcirculation DUNE directory not found: %s', micro_root);
    end

    fprintf('[CLEAN] Removing build-cmake in selected DUNE modules...\n');
    for i = 1:numel(MICRO_DUNE_MODULES_TO_CLEAN)
        modPath = fullfile(micro_root, MICRO_DUNE_MODULES_TO_CLEAN{i}, 'build-cmake');
        if isfolder(modPath)
            fprintf('  - %s\n', modPath);
            rmdir(modPath, 's');
        else
            fprintf('  - (skip) %s (not found)\n', modPath);
        end
    end

    fprintf('[CLEAN] Cleaning dune-angiogenesis/build-cmake files and dirs...\n');
    if isfolder(angio_build)
        filesToDelete = {'CMakeCache.txt','CMakeDoxyfile.in','CMakeDoxygenDefaults.cmake','CMkaeDoxygenDefaults.cmake'};
        for i = 1:numel(filesToDelete)
            fpath = fullfile(angio_build, filesToDelete{i});
            if isfile(fpath)
                fprintf('  - file %s\n', fpath);
                delete(fpath);
            else
                fprintf('  - (skip) %s (not found)\n', fpath);
            end
        end
        dirsToRemove = {'CMakeFiles','cmake'};
        for i = 1:numel(dirsToRemove)
            dpath = fullfile(angio_build, dirsToRemove{i});
            if isfolder(dpath)
                fprintf('  - dir %s\n', dpath);
                rmdir(dpath, 's');
            else
                fprintf('  - (skip) %s (not found)\n', dpath);
            end
        end
    else
        fprintf('  - (skip) %s (not found)\n', angio_build);
    end

    % 3) Run ./build.sh in MICRO_DIR
    fprintf('[BUILD] Running build.sh in: %s\n', micro_root);
    prev = pwd; c2 = onCleanup(@() cd(prev)); 
    cd(micro_root);

    if ~isfile(fullfile(micro_root,'build.sh'))
        error('build.sh not found in %s', micro_root);
    end
    system('chmod +x build.sh');

    [sB, oB] = system('./build.sh'); fprintf('%s', oB);
    if sB ~= 0, error('Microcirculation build.sh exited with status %d.', sB); end

    fprintf('[DONE] Microcirculation build completed.\n');
else
    fprintf('[SKIP] DO_BUILD_MICRO = false → skipping microcirculation build.\n');
end

%% -------- CSV HANDOFF: locate CSV and copy into micro/src --------
fprintf('\n[CSV] Looking for pica_flow_heart_period.csv...\n');

% 1) Preferred location: expected macrocirculation output dir
found_csv = '';
if isfile(macro_csv_expected)
    found_csv = macro_csv_expected;
    fprintf('[CSV] Found at expected location: %s\n', found_csv);
else
    % 2) Recursive search under Macrocirculation/build/bin/macrocirculation
    if isfolder(exe_dir)
        pattern = fullfile(exe_dir, '**', 'pica_flow_heart_period.csv');
        matches = dir(pattern);
        if ~isempty(matches)
            % pick the most recent match
            [~, idx] = max([matches.datenum]);
            found_csv = fullfile(matches(idx).folder, matches(idx).name);
            fprintf('[CSV] Found via recursive search: %s\n', found_csv);
        end
    end

    % 3) If still not found, check micro sink dir (maybe it already exists there)
    if isempty(found_csv)
        if isfile(micro_csv_sink)
            fprintf('[CSV] Already present in micro sink: %s\n', micro_csv_sink);
            found_csv = micro_csv_sink; % treat as found to avoid copy
        end
    end
end

% Act based on what we found
if isempty(found_csv)
    fprintf(2, ['[CSV] pica_flow_heart_period.csv not found in expected macrocirculation ' ...
                'output, nor under %s, nor in %s.\n' ...
                '→ Please run the macrocirculation simulation to generate it (set RUN_MACRO = true).\n'], ...
                exe_dir, micro_csv_sink_dir);
    error('Missing pica_flow_heart_period.csv; cannot proceed to Microcirculation run.');
else
    % If the found file is not already at the sink, copy it
    if ~strcmp(found_csv, micro_csv_sink)
        if ~isfolder(micro_csv_sink_dir)
            fprintf('[CSV] Creating destination dir: %s\n', micro_csv_sink_dir);
            mkdir(micro_csv_sink_dir);
        end
        fprintf('[CSV] Copying to micro sink: %s → %s\n', found_csv, micro_csv_sink);
        copyfile(found_csv, micro_csv_sink, 'f');
    end
    fprintf('[CSV] Ready: %s\n', micro_csv_sink);
end

%% -------- PREP NVU VESSEL DROP DIR (Micro <— NGVU) --------
if ~isfolder(micro_src_dir), error('Micro src dir not found: %s', micro_src_dir); end
if ~isfolder(nvu_drop_dir)
    fprintf('[NVU] Creating drop directory: %s\n', nvu_drop_dir);
    mkdir(nvu_drop_dir);
else
    % Clean older NVU_Vessel_*.csv files
    oldFiles = dir(fullfile(nvu_drop_dir, 'NVU_Vessel_*.csv'));
    if ~isempty(oldFiles)
        fprintf('[NVU] Cleaning %d old NVU_Vessel CSV files in %s\n', numel(oldFiles), nvu_drop_dir);
        for k = 1:numel(oldFiles)
            delete(fullfile(oldFiles(k).folder, oldFiles(k).name));
        end
    end
end

%% -------- START MICROCIRCULATION (dune-angiogenesis ArtificialNetwork) --------
if ~isfile(micro_exe_path)
    error('Micro executable not found: %s', micro_exe_path);
end

fprintf('\n[MICRO] Launching dune-angiogenesis %s in %s\n', MICRO_ARG, micro_src_dir);
prev = pwd; c3 = onCleanup(@() cd(prev)); 
cd(micro_src_dir);

if RUN_MICRO_IN_BACKGROUND
    % Start in background via bash -lc to capture PID
    cmd = sprintf('bash -lc ''nohup "./dune-angiogenesis" %s > %s 2>&1 < /dev/null & echo $!''', ...
                  MICRO_ARG, MICRO_STDOUT_FILE);
    [s_bg, out_bg] = system(cmd);
    if s_bg ~= 0
        error('[MICRO] Failed to start micro in background. Status=%d', s_bg);
    end
    micro_pid = str2double(strtrim(out_bg));
    if isnan(micro_pid) || micro_pid <= 0
        fprintf(2, '[MICRO] Warning: could not parse PID from output: %s\n', out_bg);
    else
        fprintf('[MICRO] Started in background. PID=%d, log=%s\n', micro_pid, fullfile(micro_src_dir, MICRO_STDOUT_FILE));
    end
else
    % Foreground (blocking) run
    cmd = sprintf('"%s" %s', micro_exe_path, MICRO_ARG);
    fprintf('[MICRO] %s\n', cmd);
    [s_fg, out_fg] = system(cmd); fprintf('%s', out_fg);
    if s_fg ~= 0, error('[MICRO] Foreground micro exited with status %d.', s_fg); end
    fprintf('[MICRO] Foreground micro finished.\n');
end

%% -------- LAUNCH NGVU MERGED DRIVER (writes NVU_Vessel_*.csv into nvu_drop_dir) --------
% The NGVU code already writes to Micro src (build-cmake/src) by default in your driver.
% If you need to enforce the exact folder, set a global there to nvu_drop_dir.
fprintf('\n[NGVU] Running merged Quadripartite+NVC driver: %s\n', ngvu_main_path);
if ~isfile(ngvu_main_path)
    error('NGVU driver not found: %s', ngvu_main_path);
end

% Add NGVU to path temporarily and run
addpath(fullfile(repo_root, NGVU_DIR));
try
    run(ngvu_main_path);
catch ME
    rmpath(fullfile(repo_root, NGVU_DIR));
    rethrow(ME);
end
rmpath(fullfile(repo_root, NGVU_DIR));

fprintf('\n[MASTER] NGVU driver finished. Master script done.\n');
