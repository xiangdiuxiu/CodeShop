tar -zcvf mvpart.tar.gz mvpart/
R CMD REMOVE mvpart
R CMD INSTALL mvpart.tar.gz
