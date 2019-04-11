function [output] = detection_eval(y, threshold)

output = false(length(y), 1);
output(y > threshold) = true;
