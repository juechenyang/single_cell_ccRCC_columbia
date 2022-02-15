source("call_libraries.R")
all_raw_file_names = list.files("./Michigan_raw/", full.names = T)
a = Read10X_h5(all_raw_file_names[[1]])
sa = CreateSeuratObject(a)
