function out = DataType(sf, bpp)
switch sf
    case 1
        switch bpp
            case 8
                out = 'uint8';
            case 16
                out = 'uint16';
            case 32
                out = 'uint32';
        end
    case 2
        switch bpp
            case 8
                out = 'int8';
            case 16
                out = 'int16';
            case 32
                out = 'int32';
        end
    case 3
        switch bpp
            case 32
                out = 'single';
            case 64
                out = 'double';
        end
end
end