-- how many solver iterations to run.  right now it's diverging, so ... it's off
maxiter = 100

initCond = 'flat'
--initCond = 'stellar_schwarzschild'
--initCond = 'stellar_kerr_newman'

--solver = 'conjgrad'
--solver = 'conjres'
--solver = 'gmres'
solver = 'jfnk'

size = 4

-- how big to construct the grid
bodyRadii = 2

outputFilename = 'out.txt'
