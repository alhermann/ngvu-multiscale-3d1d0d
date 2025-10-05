function [out]=Hill(dezelfde,constante,power)
    out = (dezelfde^power)/(dezelfde^power + constante^power);
end