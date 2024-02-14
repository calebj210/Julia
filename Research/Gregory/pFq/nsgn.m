% Modified sign function to return -1 when z = 0
function result = nsgn(z)
if z == 0
    result = -1;
else
    result = sign(z);
end
end