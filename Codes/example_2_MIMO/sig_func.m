function out = sig_func(signal,power)

    out = (abs(signal).^power).*sign(signal);

end