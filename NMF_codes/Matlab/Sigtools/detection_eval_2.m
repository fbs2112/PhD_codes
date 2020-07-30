function [output] = detection_eval_2(y, threshold)


if any(y > threshold)
    output = true;
else
    output = false;
end

