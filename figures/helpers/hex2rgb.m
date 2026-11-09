function rgb = hex2rgb(hexStr)
    hexStr = char(hexStr);
    if hexStr(1) == '#', hexStr = hexStr(2:end); end
    rgb = reshape(sscanf(hexStr, '%2x'), 3, []).'/255;
end