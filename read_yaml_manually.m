function sensor_data = read_yaml_manually(filename)
    % Read the entire YAML file into a string
    fid = fopen(filename, 'r');
    yaml_content = fread(fid, '*char')';
    fclose(fid);
    
    % Extract the fields manually using regular expressions
    sensor_data.sensor_type = extract_value(yaml_content, 'sensor_type');
    sensor_data.comment = extract_value(yaml_content, 'comment');
    
    % Extract the transformation matrix T_BS
    T_BS_str = extract_value(yaml_content, 'data: \[(.*?)\]', 'match');
    T_BS_data = str2num(T_BS_str);
    sensor_data.T_BS = reshape(T_BS_data, 4, 4)';
    
    % Extract the camera specs
    sensor_data.rate_hz = str2num(extract_value(yaml_content, 'rate_hz'));
    sensor_data.resolution = str2num(extract_value(yaml_content, 'resolution: \[(.*?)\]', 'match'){1});
    sensor_data.camera_model = extract_value(yaml_content, 'camera_model');
    sensor_data.intrinsics = str2num(extract_value(yaml_content, 'intrinsics: \[(.*?)\]', 'match'){1});
    sensor_data.distortion_model = extract_value(yaml_content, 'distortion_model');
    sensor_data.distortion_coefficients = str2num(extract_value(yaml_content, 'distortion_coefficients: \[(.*?)\]', 'match'){1});
end

% Helper function to extract values from YAML content using regular expressions
function value = extract_value(yaml_content, pattern, option = 'once')
    tokens = regexp(yaml_content, pattern, 'tokens', option);
    if iscell(tokens) && !isempty(tokens)
        value = strtrim(tokens{1}{1});
    else
        value = '';
    end
end
