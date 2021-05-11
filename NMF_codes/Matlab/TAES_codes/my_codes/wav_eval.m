function wav_output = wav_eval(signal)

%------1 stage-------------------------
[a1,a2] = dwt(signal, 'dmey');

%------2 stage-------------------------
[b1,b2] = dwt(a1, 'dmey');
[b3,b4] = dwt(a2, 'dmey');

%------3 stage-------------------------
[c1,c2] = dwt(b1, 'dmey');
[c3,c4] = dwt(b2, 'dmey');
[c5,c6] = dwt(b3, 'dmey');
[c7,c8] = dwt(b4, 'dmey');

%------4 stage-------------------------

[d1,d2] = dwt(c1, 'dmey');
[d3,d4] = dwt(c2, 'dmey');

[d5,d6] = dwt(c3, 'dmey');
[d7,d8] = dwt(c4, 'dmey');

[d9,d10] = dwt(c5, 'dmey');
[d11,d12] = dwt(c6, 'dmey');

[d13,d14] = dwt(c7, 'dmey');
[d15,d16] = dwt(c8, 'dmey');

wav_output(:,1) = d1;
wav_output(:,2) = d2;
wav_output(:,3) = d3;
wav_output(:,4) = d4;
wav_output(:,5) = d5;
wav_output(:,6) = d6;
wav_output(:,7) = d7;
wav_output(:,8) = d8;
wav_output(:,9) = d9;
wav_output(:,10) = d10;
wav_output(:,11) = d11;
wav_output(:,12) = d12;
wav_output(:,13) = d13;
wav_output(:,14) = d14;
wav_output(:,15) = d15;
wav_output(:,16) = d16;

wav_output = abs(wav_output);