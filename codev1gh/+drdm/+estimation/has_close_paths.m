function close = has_close_paths(paths,context)
%HAS_CLOSE_PATHS 检查尺度化坐标中过于接近的路径。
%
% 输入：paths、context：路径和参数尺度。
% 输出：close：是否存在距离小于 duplicate_radius 的径对。
%
% 与候选去重使用相同的时延和多普勒尺度。

close = false;
for i = 1:numel(paths)
    for j = i+1:numel(paths)
        distance = hypot((paths(i).tau_s-paths(j).tau_s)/context.delay_scale_s, ...
            (paths(i).nu_hz-paths(j).nu_hz)/context.doppler_scale_hz);
        if distance<context.options.duplicate_radius
            close = true;
            return;
        end
    end
end
end
