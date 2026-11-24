function value = code_signature()
%CODE_SIGNATURE 计算当前函数库源码与 MATLAB 版本标识。
%
% 输入：无。
% 输出：value：SHA-256 处理版本标识。
%
% 相对文件名排序后计算，便于在不同目录核对同一版本。

packageRoot = fileparts(fileparts(mfilename('fullpath')));
files = dir(fullfile(packageRoot,'**','*.m'));
records = repmat(struct('path','','text',''),1,numel(files));
for k = 1:numel(files)
    name = fullfile(files(k).folder,files(k).name);
    records(k) = struct('path',strrep(name(numel(packageRoot)+2:end),'\','/'),'text',fileread(name));
end
[~,order] = sort({records.path});
value = drdm.statistics.fingerprint(struct('matlab',version,'source',records(order)));
end
