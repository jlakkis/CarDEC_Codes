sct_norm = function(array, cellnames, genenames){
    rownames(array) = cellnames
    colnames(array) = genenames
    array = t(array)
    array = sctransform::vst(array)$y
    return(t(array))
}