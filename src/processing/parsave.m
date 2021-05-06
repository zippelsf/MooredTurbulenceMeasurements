
function [] = parsave(fname, C)
    save(fname,  '-struct', 'C', '-append')
end
