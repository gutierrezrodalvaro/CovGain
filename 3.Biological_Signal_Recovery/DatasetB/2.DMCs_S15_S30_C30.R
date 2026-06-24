# With this script We determine which CpG is each methodology clasiffiying as DMCs
# Data used to do the Figure 3
# Alvaro Gutierrez Rodriguez
 
library(methylKit)

load('meth_objs_min5.Rdata')

#  "CTM15x_meth" "CTM30x_meth" "STD15x_meth" "STD30x_meth"

# STD15x_meth
S15_diff <- calculateDiffMeth(STD15x_meth, mc.cores = 2)
S15_DMCs <- getMethylDiff(S15_diff, qvalue = 0.05, difference = 15)

# CTM15x_meth
C15_diff <- calculateDiffMeth(CTM15x_meth, mc.cores = 2)
C15_DMCs <- getMethylDiff(C15_diff, qvalue = 0.05, difference = 15)

# CTM30x_meth
C30_diff <- calculateDiffMeth(CTM30x_meth, mc.cores = 2)
C30_DMCs <- getMethylDiff(C30_diff, qvalue = 0.05, difference = 15)

# CTM15x_meth
S30_diff <- calculateDiffMeth(STD30x_meth, mc.cores = 2)
S30_DMCs <- getMethylDiff(S30_diff, qvalue = 0.05, difference = 15)

save(S15_diff, C15_diff, C30_diff, S30_diff, S15_DMCs, C15_DMCs, C30_DMCs, S30_DMCs, file = 'DMCs_F5.Rdata')


