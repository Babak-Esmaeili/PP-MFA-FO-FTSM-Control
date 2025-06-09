function out = sat_func(s,sigma)

    out = (s./sigma).*(abs(s)<=sigma) + (sign(s)).*(abs(s)>sigma);

end