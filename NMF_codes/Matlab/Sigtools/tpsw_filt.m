function output = tpsw_filt(x, M, N, alpha)

tpsw_coeffs = zeros(N-1, 1);
tpsw_coeffs(M:end) = 1;
tpsw_coeffs = 1/(2*(N-M))*[tpsw_coeffs(end:-1:1);0;tpsw_coeffs];

freqComponent = fftfilt(tpsw_coeffs, x);
background = x;
background(x > alpha*freqComponent) = freqComponent(x > alpha*freqComponent);
background = fftfilt(tpsw_coeffs, background);

output = x - background;
output(output < 0) = eps;