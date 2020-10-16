cbaywet = read.table('cbaywet')
dim(cbaywet)
cbaywet = cbaywet[,-1]
write.table(cbaywet, 'cbayJulSep0019', row.names = F, col.names = F)
