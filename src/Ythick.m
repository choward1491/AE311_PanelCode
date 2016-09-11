%Project functions
function output = Ythick(x)
output = (.08/.2)*(.29690*sqrt(x)-.126.*x-.35160.*x.^2+.28430.*x.^3-.10150.*x.^4);
end
