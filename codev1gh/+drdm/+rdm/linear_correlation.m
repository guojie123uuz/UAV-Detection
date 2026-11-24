function values = linear_correlation(y,p,fs,lags,doppler,referenceOffset,method)
%LINEAR_CORRELATION 计算带多普勒补偿的线性互相关。
%
% 输入：y、p：列向量；fs：采样率；lags：整数时延；doppler：Hz；referenceOffset：参考偏移；method：fft/direct。
% 输出：values：多普勒×时延的复相关矩阵。
%
% C(l,f)=Σ y[n]conj(p[n+offset-l])exp(-j2πft[n])/fs，t 以 y 中心为原点。

validateattributes(y,{'numeric'},{'column','nonempty','finite'});
validateattributes(p,{'numeric'},{'column','nonempty','finite'});
validateattributes(fs,{'numeric'},{'scalar','positive','real','finite'});
validateattributes(referenceOffset,{'numeric'},{'scalar','integer','finite'});
[lags,doppler] = drdm.rdm.check_axes(lags,doppler);
assert(all(lags==round(lags)),'drdm:SampleLag','IQ correlation needs integer sample lags.');
assert(all(abs(doppler)<=fs/2),'drdm:DopplerAlias','Doppler exceeds sample Nyquist range.');
n = numel(y);
m = numel(p);
t = ((0:n-1)'-(n-1)/2)/fs;
values = complex(zeros(numel(doppler),numel(lags)));
switch lower(char(method))
    case 'fft'
        lengthFFT = 2^nextpow2(n+m-1);
        spectrum = conj(fft(p,lengthFFT));
        relativeLags = lags-referenceOffset;
        supported = relativeLags>=-(m-1) & relativeLags<=n-1;
        indices = mod(relativeLags(supported),lengthFFT)+1;
        for f = 1:numel(doppler)
            c = ifft(fft(y.*exp(-1i*2*pi*doppler(f)*t),lengthFFT).*spectrum)/fs;
            values(f,supported) = c(indices).';
        end
    case 'direct'
        % 直接求和逐时延处理，避免构造完整的时延×多普勒×时间张量。
        for k = 1:numel(lags)
            first = max(1,1-referenceOffset+lags(k));
            last = min(n,m-referenceOffset+lags(k));
            indices = (first:last)';
            product = y(indices).*conj(p(indices+referenceOffset-lags(k)));
            for f = 1:numel(doppler)
                values(f,k) = sum(product.*exp(-1i*2*pi*doppler(f)*t(indices)))/fs;
            end
        end
    otherwise
        error('drdm:CorrelationMethod','Method must be fft or direct.');
end
end
