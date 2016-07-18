maxiter = 0

--initCond = 'flat'
--initCond = 'stellar_schwarzschild'
initCond = 'stellar_kerr_newman'

--solver = 'conjgrad'
--solver = 'conjres'
--solver = 'gmres'
solver = 'jfnk'

size = 64

-- how big to construct the grid
bodyRadii = 2

outputFilename = 'out.txt'
