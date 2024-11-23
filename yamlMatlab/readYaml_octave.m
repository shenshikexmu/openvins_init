function results = readYaml_octave(filePath)
% Lloyd Russell 2017
% Simple little function to read simple little YAML format parameter files
% from  https://github.com/llerussell/ReadYAML

% read file line by line
fid = fopen(filePath, 'r');
data = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
fclose(fid);

% remove empty lines
data = deblank(data{1});
data(cellfun('isempty', data)) = [];

% prepare final results structure
results = [];

% parse the contents (line by line)
for i = 1:numel(data)
    
    % extract this line
    thisLine = data{i};
    
    % ignore if this line is a comment, start or end of document
    if strcmpi(thisLine(1), '#') || strcmpi(thisLine(1:3), '---') || strcmpi(thisLine(1:3), '...')
        continue
    end
    
    % find the seperator between key and value
    sepIndex = find(thisLine==':', 1, 'first');
    
    % get the key name (remove whitespace)
    key = strtrim(thisLine(1:sepIndex-1));  
    
    % get the value, ignoring any comments (remove whitespace)
    value = strsplit(thisLine(sepIndex+1:end), '#');
    value = strtrim(value{1}); 
    
    % attempt to convert value to numeric type
    [convertedValue, success] = str2num(value);
    if success
        value = convertedValue;
    end
    
    % store the key and value in the results
    results.(key) = value;
end