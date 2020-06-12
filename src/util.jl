
# Converts a MAT file to an HDF5 file
function mat_to_h5(fn, out)
    matdata = matread(fn)
    h5open(out, "w") do file
        for (k, v) in matdata
            write(file, k, size(v, 2) == 1 ? vec(v) : v)
        end
    end
end
