-- grid size
--size = 1
--size = 4
size = 8
--size = {16, 1, 1}
--size = 64
-- 10*8^3 = 5120

-- body specifies the radius of the problem, and the initial stress energy primitives
--body = 'Null'
body = 'earth'
--body = 'sun'
--body = 'EMUniformField'
--body = 'em_line'

-- this says how much bigger is our grid than our body radius
bodyRadii = 2

-- initCond specifies the inital metric primitives
initCond = 'flat'
--initCond = 'stellar_schwarzschild'
--initCond = 'stellar_kerr_newman'
--initCond = 'EMUniformField'
--initCond = 'em_line'

--solver = 'conjgrad'
--solver = 'conjres'
--solver = 'gmres'
solver = 'jfnk'

-- how many solver iterations to run.
-- right now all linear solvers fail (maybe because I'm adjusting b mid-step?)
-- the jfnk will run one or two iterations, but always converge to prims=0
-- and the efe error looks worse with a flat solution than it does with a Schwarzschild stellar solution 
-- ... I think fixing that last one will help the results converge.
maxiter = 1000

outputFilename = 'out.txt'
