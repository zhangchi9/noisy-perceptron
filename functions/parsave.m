function parsave(varargin)

for i = 1:nargin
    s = inputname(i);
    v = genvarname(s);
    eval([v '= varargin{i};']);
end
save(filename)