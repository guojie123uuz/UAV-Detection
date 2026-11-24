function y = apply_path(source,path,options)
%APPLY_PATH 施加有符号分数时延、多普勒和复增益。
%
% 输入：source：带前后保护样本的源；path：tau_s、nu_hz、gain_complex；options：插值参数。
% 输出：y：CPI 内 N×1 复数路径分量。
%
% 分数时延采用 Kaiser 加窗 sinc 插值，多普勒相位以 CPI 中心为零点。

if nargin<3
    options = struct('radius',64,'kaiser_beta',10);
end
validateattributes(source.iq,{'numeric'},{'column','nonempty','finite'});
validateattributes(source.fs_hz,{'numeric'},{'scalar','real','positive','finite'});
validateattributes(source.cpi_first_sample,{'numeric'},{'scalar','integer','positive'});
validateattributes(source.cpi_num_samples,{'numeric'},{'scalar','integer','positive'});
validateattributes(path.tau_s,{'numeric'},{'scalar','real','finite'});
validateattributes(path.nu_hz,{'numeric'},{'scalar','real','finite'});
validateattributes(path.gain_complex,{'numeric'},{'scalar','finite'});
validateattributes(options.radius,{'numeric'},{'scalar','integer','>=',4});
validateattributes(options.kaiser_beta,{'numeric'},{'scalar','real','finite','positive'});
% 将有符号时延分为整数移位与 [0,1) 范围的分数部分。
delay = path.tau_s*source.fs_hz;
lag = floor(delay);
fraction = delay-lag;
if abs(delay-round(delay))<1e-10
    lag = round(delay);
    fraction = 0;
end
indices = source.cpi_first_sample+(0:source.cpi_num_samples-1)-lag;
radius = options.radius*(fraction~=0);
assert(indices(1)-radius>=1 && indices(end)+radius<=numel(source.iq), ...
    'drdm:SourceGuard','Insufficient source history/future for this signed delay.');
% 仅提取输出 CPI 所需的局部支撑，减少对保护区的重复卷积。
first = indices(1)-radius;
last = indices(end)+radius;
x = double(source.iq(first:last));
localIndices = radius+(1:source.cpi_num_samples);
if fraction~=0
    k = (-radius:radius)';
    z = k-fraction;
    h = sin(pi*z)./(pi*z);
    h = h.*besseli(0,options.kaiser_beta*sqrt(max(0,1-(k/radius).^2)))/ ...
        besseli(0,options.kaiser_beta);
    h = h/sum(h);
    x = conv(x,h,'same');
end
% 多普勒相位参考 CPI 中心，使复增益在不同拟合步骤中保持同一约定。
t = ((0:source.cpi_num_samples-1)'-(source.cpi_num_samples-1)/2)/source.fs_hz;
y = path.gain_complex*x(localIndices).*exp(1i*2*pi*path.nu_hz*t);
end
