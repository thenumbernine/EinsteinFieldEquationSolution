-- how many solver iterations to run.
-- right now all linear solvers fail (maybe because I'm adjusting b mid-step?)
-- the jfnk will run one or two iterations, but always converge to prims=0
-- and the efe error looks worse with a flat solution than it does with a Schwarzschild stellar solution 
-- ... I think fixing that last one will help the results converge.
maxiter = 0

-- this holds things like radius of the grid and stress-energy prim information 
body = 'earth'
--body = 'sun'
--body = 'empty'

--initCond = 'flat'
initCond = 'stellar_schwarzschild'
--initCond = 'stellar_kerr_newman'
--initCond = 'rotating_EM_field'

--solver = 'conjgrad'
--solver = 'conjres'
--solver = 'gmres'
solver = 'jfnk'

size = 64

-- how much bigger is our grid than our spherical body 
bodyRadii = 2

outputFilename = 'out.txt'
