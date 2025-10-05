% ---------- utilities ----------
function out = negCheck(input, width)
% Convert an extensive quantity 'input' [uM*m] into concentration [uM]
% using compartment width 'width' [m]. Prevent negatives underflow.
    if (width > 0) && ((input / width) > 0)
        out = input / width;
    else
        out = 1e-180;
    end
end