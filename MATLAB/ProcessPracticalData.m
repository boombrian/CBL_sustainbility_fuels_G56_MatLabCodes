clear; 
clc;

% Specify folder
outputFolder = 'Data';

% If the folder doesnt exist, then create it
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
    fprintf('Created output folder: %s\n', outputFolder);
end


% Load slowDAQ data
fprintf('\n Select your slowDAQ data file...\n');
[slowFileName, slowFilePath] = uigetfile('*.txt', 'Select the slowDAQ File');

% Stop script if file selection is closed
if isequal(slowFileName, 0)
    error('File selection canceled by user.');
end

% Combine folder and file name to create the full path and read it
slowFullPath = fullfile(slowFilePath, slowFileName);
slowDAQ = readmatrix(slowFullPath);

% Remove any rows that contain invalid numbers
slowDAQ = slowDAQ(all(isfinite(slowDAQ), 2), :);

% Check that slowDAQ has 5 columns
[~, numSlowCols] = size(slowDAQ);
if numSlowCols ~= 5
    error('slowDAQ must have 5 columns, the loaded file has %d columns.', numSlowCols);
end


% Load fastDAQ data - Same procedure as with the slowDAQ data
fprintf(' Select your fastDAQ data file...\n');
[fastFileName, fastFilePath] = uigetfile('*.txt', 'Select the fastDAQ File');

if isequal(fastFileName, 0)
    error('File selection canceled by user.');
end

fastFullPath = fullfile(fastFilePath, fastFileName);
fastDAQ = readmatrix(fastFullPath);

fastDAQ = fastDAQ(all(isfinite(fastDAQ), 2), :);

% Check that fastDAQ has 3 columns
[numFastRows, numFastCols] = size(fastDAQ);
if numFastCols ~= 3
    error('fastDAQ must have 3 columns, the loaded file has %d columns.', numFastCols);
end


% Clean up the cycles and crank angles
pointsPerCycle = 3600;
expectedCA = (-360 : 0.2 : 359.8)'; % Created crank angle list we will replace with at the end

% Verify the total number of rows is a multiple of cycles completed
totalCycles = numFastRows / pointsPerCycle;
if mod(totalCycles, 1) ~= 0
    error('fastDAQ row count (%d) is not an exact multiple of 3,600 points.', numFastRows);
end

% Iterate through each cycle
for i = 1:totalCycles
    % Calculate row indices for the current cycle
    startIndex = (i - 1) * pointsPerCycle + 1;
    endIndex   = i * pointsPerCycle;
    
    % Extract the current cycle data
    cycleData = fastDAQ(startIndex:endIndex, :);
    
    % Map the cycle data to the range we use [-360, +360)
    cycleData(:, 1) = mod(cycleData(:, 1) + 360, 720) - 360;
    
    % Sort the rows by crank angle chronologically
    cycleData = sortrows(cycleData, 1);
    
    % Put our list of crank angles we created in column 1
    cycleData(:, 1) = expectedCA;
    
    % Put cleaned cycle data back into fastDAQ
    fastDAQ(startIndex:endIndex, :) = cycleData;
end


% Export cleaned data
userAnswer = inputdlg('Enter filename:', ...
                      'Save Cleaned Data', [1 60], {''});

% Stop script if canceled
if isempty(userAnswer) || isempty(strtrim(userAnswer{1}))
    error('Export canceled or no filename was entered.');
end

% Change name so that it is not invalid
runID = matlab.lang.makeValidName(userAnswer{1});

% Make the output text file paths
slowOutFile = fullfile(outputFolder, [runID '_slowdaq.txt']);
fastOutFile = fullfile(outputFolder, [runID '_fastdaq.txt']);

% Save the text files
writematrix(slowDAQ, slowOutFile, 'Delimiter', 'space');
writematrix(fastDAQ, fastOutFile, 'Delimiter', 'space');

fprintf('\nData processing complete!\n');
fprintf('Saved: %s\n', slowOutFile);
fprintf('Saved: %s\n', fastOutFile);