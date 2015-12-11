% Matrix to string in WKT (Well-Known-Text) format
% 2015, Manuel Claeys Bouuaert

function strWKT=mat2strWKT(mat,n)
    if nargin < 2
        n = 15;
    end
    str = mat2str(mat,n);
    str = strrep(str, ';', ',');
    str = strrep(str, '[', '(');
    strWKT = strrep(str, ']', ')');
end