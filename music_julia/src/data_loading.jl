Get_SC_Meta_Info(file) = DataFrame(readdlm(file, '\t', header=true))
Get_SC_Count_Info(file) = DataFrame(readdlm(file, '\t', header=true, footerskip=1))
Get_bulk_Meta_Info(file) = DataFrame(readdlm(file, '\t', header=true))
Get_bulk_Count_Info(file) = DataFrame(readdlm(file, '\t', header=true, footerskip=1))