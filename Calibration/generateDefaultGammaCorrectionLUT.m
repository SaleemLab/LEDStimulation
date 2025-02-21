% Define LUT parameters
LUT_size = 1041;             % Total number of elements (0 to 1040)
LUT_values = 0:(LUT_size-1);   % Create an array [0 1 2 ... 1040]
outputFileName = 'defaultLUT.txt';

% Open file for writing
fileID = fopen(outputFileName, 'w');
if fileID == -1
    error('Cannot open output file.');
end

% Write the array declaration header
fprintf(fileID, 'const uint16_t defaultLUT[%d] = {\n', LUT_size);

% Write each LUT value (10 values per line for readability)
for i = 1:length(LUT_values)
    fprintf(fileID, ' %d', LUT_values(i));
    if i < length(LUT_values)
        fprintf(fileID, ',');
    end
    % Insert a newline every 10 values
    if mod(i, 10) == 0
        fprintf(fileID, '\n');
    end
end

% Close the array declaration
fprintf(fileID, '\n};\n');
fclose(fileID);

disp(['Default LUT saved to ', outputFileName]);