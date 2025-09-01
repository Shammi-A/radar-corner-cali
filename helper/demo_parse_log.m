function out = demo_parse_log(rootDir, testSel)
% DEMO_PARSE_LOG  Pick test (T1/T2/...), find latest .txt, parse, and list associated .bin files.
% Returns: struct with test info, log path, bin files, RadarParams, Raw.

    if nargin < 1 || isempty(rootDir) || ~isfolder(rootDir)
        fprintf(2,'Root not found → pick folder that contains T1/T2/... \n');
        picked = uigetdir(pwd, 'Select calibration root (contains T1, T2, ...)');
        if isequal(picked,0), error('No folder selected'); end
        rootDir = picked;
    end

    % ---- discover tests ----
    tests = listTestFolders(rootDir);
    if isempty(tests), error('No "T#" or "Test #" folders under: %s', rootDir); end

    % ---- choose test ----
    if nargin < 2 || isempty(testSel)
        fprintf('Discovered tests:\n');
        for i=1:numel(tests), fprintf('  [%2d] %s\n', i, tests(i).name); end
        k = input(sprintf('Select test [1..%d] (default %d): ', numel(tests), numel(tests)));
        if isempty(k), k = numel(tests); end
        testIdx = k;
    else
        testIdx = matchTestSelection(tests, testSel);
        if isnan(testIdx), error('Could not match "%s" to available tests.', string(testSel)); end
    end

    testFolder = tests(testIdx).folder;
    testName   = tests(testIdx).name;

    % ---- latest .txt in this test (recursive) ----
    txtFiles = rdir(fullfile(testFolder, '**', '*.txt'));
    if isempty(txtFiles), error('No .txt logs inside: %s', testFolder); end
    [~, j] = max([txtFiles.datenum]);
    logPath = fullfile(txtFiles(j).folder, txtFiles(j).name);

    % ---- parse log (you can pass preferred profileId here if needed) ----
    % opts.profileId = 0;  % uncomment to force profile 0
    [RadarParams, Raw] = mmws_parse_log(logPath);

    % ---- find StartRecord base & associated BINs ----
    [binBasePath, baseStem] = parseStartRecordBase(logPath);
    binTbl = findAssociatedBins(binBasePath, baseStem, testFolder);

    % ---- summary ----
    fprintf('\n=== Test: %s ===\n', testName);
    fprintf('Folder : %s\n', testFolder);
    fprintf('Log    : %s\n', logPath);
    if ~isempty(binTbl)
        disp(binTbl(:, {'rawIdx','sizeMB','path'}));
    else
        warning('No associated *_Raw_*.bin files found.');
    end
    mmws_print_summary(RadarParams);

    out = struct('testName',testName, 'testFolder',testFolder, 'logPath',logPath, ...
                 'binBasePath',binBasePath, 'binFiles',binTbl, ...
                 'RadarParams',RadarParams, 'Raw',Raw);
end

% ----------------- local helpers -----------------
function tests = listTestFolders(rootDir)
    D = dir(rootDir); D = D([D.isdir]);
    names = {D.name};
    keep  = false(size(names));
    for i=1:numel(names)
        n = names{i};
        if any(strcmp(n,{'.','..'})), continue; end
        if ~isempty(regexpi(n,'^T\d+$')) || ~isempty(regexpi(n,'^Test\s*\d+$'))
            keep(i)=true;
        end
    end
    D = D(keep);
    for i=1:numel(D), D(i).folder = fullfile(rootDir, D(i).name); end
    idx = arrayfun(@(d) extractTestNumber(d.name), D);
    [~,ord] = sort(idx);
    tests = struct('name',{D(ord).name}, 'folder',{D(ord).folder});
end

function k = matchTestSelection(tests, sel)
    if isnumeric(sel), want = double(sel);
    else, m = regexp(char(string(sel)),'(\d+)','tokens','once');
         if isempty(m), k=NaN; return; end
         want = str2double(m{1});
    end
    nums = arrayfun(@(t) extractTestNumber(tests(t).name), 1:numel(tests));
    k = find(nums==want,1,'first'); if isempty(k), k=NaN; end
end

function n = extractTestNumber(name)
    m = regexp(name,'(\d+)','tokens','once');
    if isempty(m), n=NaN; else, n=str2double(m{1}); end
end

function files = rdir(pat)
    files = dir(pat); files = files(~[files.isdir]);  % supports ** in modern MATLAB
end

function [basePath, baseStem] = parseStartRecordBase(logPath)
    txt = fileread(logPath);
    basePath = ''; baseStem='';
    m = regexp(txt, 'CaptureCardConfig_StartRecord\("([^"]+)\.bin"\s*,\s*\d+\)', 'tokens', 'once');
    if isempty(m), return; end
    basePath = m{1};
    [~, baseStem, ~] = fileparts([m{1} '.bin']);
end

function T = findAssociatedBins(basePath, baseStem, fallbackRoot)
    T = table();
    if ~isempty(basePath)
        [prefDir, stem] = fileparts(basePath);
        if isempty(baseStem), baseStem = stem; end
        F = dir(fullfile(prefDir, sprintf('%s_Raw_*.bin', baseStem)));
        if isempty(F), F = dir(fullfile(prefDir, sprintf('%s*.bin', baseStem))); end
        T = makeBinTable(F);
    end
    if isempty(T)
        F = dir(fullfile(fallbackRoot, '**', '*_Raw_*.bin'));
        if isempty(F), F = dir(fullfile(fallbackRoot, '**', '*.bin')); end
        T = makeBinTable(F);
    end
    if ~isempty(T), T = sortrows(T, {'rawIdx','datenum'}); end
end

function T = makeBinTable(F)
    if isempty(F), T=table(); return; end
    n = numel(F);
    path  = strings(n,1); sizeMB=zeros(n,1); rawIdx=NaN(n,1); dn=zeros(n,1);
    for i=1:n
        path(i) = string(fullfile(F(i).folder, F(i).name));
        sizeMB(i) = F(i).bytes/(1024*1024);
        dn(i) = F(i).datenum;
        m = regexp(F(i).name,'_Raw_(\d+)\.bin$','tokens','once');
        if ~isempty(m), rawIdx(i) = str2double(m{1}); end
    end
    T = table(path,sizeMB,rawIdx,dn,'VariableNames',{'path','sizeMB','rawIdx','datenum'});
end
