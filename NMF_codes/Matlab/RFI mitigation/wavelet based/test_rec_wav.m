clear;
clc;
close all;

signal = randn(10000, 1) + 1j*randn(10000, 1);

wav_out = wav_dec_eval(signal);

[a,d] = dwt(signal, 'dmey');

xrec = idwt(a, d, 'dmey');

signalHat = wav_rec_eval(wav_out);

diffVec = signal - signalHat(1:length(signal));
errorRec = (diffVec)'*(diffVec)/length(diffVec);