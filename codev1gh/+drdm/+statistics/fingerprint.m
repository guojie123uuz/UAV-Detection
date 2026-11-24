function value = fingerprint(description)
%FINGERPRINT 计算配置结构的 SHA-256 标识。
%
% 输入：description：数值、结构体或元胞组成的配置。
% 输出：value：十六进制标识字符串。
%
% 结构字段排序后序列化；复数数组分别记录实部和虚部。

bytes = unicode2native(jsonencode(json_value(description)),'UTF-8');
digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(bytes);
value = lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end

function out = json_value(value)
out = value;
if isnumeric(value) && ~isreal(value)
    out = struct('complex_array',true,'shape',size(value),'real',real(value),'imaginary',imag(value));
elseif isstruct(value)
    out = orderfields(value);
    names = fieldnames(out);
    for element = 1:numel(out)
        for k = 1:numel(names)
            out(element).(names{k}) = json_value(out(element).(names{k}));
        end
    end
elseif iscell(value)
    for k = 1:numel(value)
        out{k} = json_value(value{k});
    end
end
end
