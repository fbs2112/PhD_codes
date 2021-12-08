function wav_output = wav_rec_eval(signal)

for i = 1:8
    c(:,i) = idwt(signal(:,(i-1)*2 + 1), signal(:,i*2), 'dmey');
end

for i = 1:4
    b(:,i) = idwt(c(:,(i-1)*2 + 1), c(:,i*2), 'dmey');
end

for i = 1:2
    a(:,i) = idwt(b(:,(i-1)*2 + 1), b(:,i*2), 'dmey');
end

wav_output = idwt(a(:,1), a(:,2), 'dmey');