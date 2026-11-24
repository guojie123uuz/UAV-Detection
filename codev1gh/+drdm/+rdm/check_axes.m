function [delay,doppler] = check_axes(delay,doppler)
%CHECK_AXES 统一时延和多普勒坐标的排列。
%
% 输入：delay、doppler：严格递增的实数坐标向量。
% 输出：delay：行向量；doppler：列向量。
%
% 所有 RDM 矩阵按“多普勒行、时延列”存储。

validateattributes(delay,{'numeric'},{'vector','nonempty','finite','real'});
validateattributes(doppler,{'numeric'},{'vector','nonempty','finite','real'});
delay = double(delay(:).');
doppler = double(doppler(:));
assert(all(diff(delay)>0) && all(diff(doppler)>0), ...
    'drdm:Axes','Axes must be unique and strictly increasing.');
end
