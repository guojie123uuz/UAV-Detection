function w = window_vector(name,count)
%WINDOW_VECTOR 生成矩形窗、Hann 窗或 Hamming 窗。
%
% 输入：name：窗名称；count：样本数。
% 输出：w：count×1 实数窗向量。
%
% Hann 和 Hamming 采用对称窗定义。

validateattributes(count,{'numeric'},{'scalar','integer','positive'});
switch lower(char(name))
    case 'rectangular'
        w = ones(count,1);
    case 'hann'
        if count==1
            w = 1;
        else
            w = 0.5-0.5*cos(2*pi*(0:count-1)'/(count-1));
        end
    case 'hamming'
        if count==1
            w = 1;
        else
            w = 0.54-0.46*cos(2*pi*(0:count-1)'/(count-1));
        end
    otherwise
        error('drdm:Window','Unknown observation window: %s.',char(name));
end
end
