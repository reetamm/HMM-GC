cbaywet = read.table('datafiles/cbaywet')
dim(cbaywet)
cbaywet = cbaywet[,-1]
write.table(cbaywet, 'datafiles/cbayJulSep0019', row.names = F, col.names = F)
