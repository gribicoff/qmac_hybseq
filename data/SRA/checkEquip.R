# check whether sequencing equipment is biased by latitude
library(openxlsx)

datt <- read.delim('SRA_metadata-hybseqv2.txt', header = T)
datx <- read.xlsx('Plant.1.0-Hyb-Seqv2.xlsx', startRow=13)

datx$lat <- substr(datx$lat_lon,1, 8) |> as.numeric()
names(datx) <- gsub('*', '', names(datx), fixed = T)
row.names(datx) <- datx$sample_name

datt$latitude <- datx[datt$sample_name, 'lat']
datt$species <- datx[datt$sample_name, 'organism']

pdf('latitudeVinstrument.pdf', 11, 8.5)
layout(matrix(1:2, 1))
boxplot(latitude ~ instrument_model, 
    data = datt, main = 'All taxa')
boxplot(latitude ~ instrument_model, 
    data = datt[datt$species == 'Quercus macrocarpa', ], main = 'Q. macrocarpa only')
dev.off()

sink(file = 't.tests.txt')

print('QUERCUS MACROCARPA ONLY:')
t.test(latitude ~ instrument_model, 
    data = datt[datt$species == 'Quercus macrocarpa', ]) |> print()

print('ALL TAXA:')
t.test(latitude ~ instrument_model, 
    data = datt) |> print()

sink(file = NULL)
