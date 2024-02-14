%% getPathIndices function
function indices = getPathIndices(z, g)
indices = getLinearIndices(getIdx(z, g), g);
end