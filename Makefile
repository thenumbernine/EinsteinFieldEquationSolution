DIST_FILENAME=EFESoln
DIST_TYPE=app
include ../Common/Base.mk
include ../Tensor/Include.mk
include ../Solvers/Include.mk
include ../LuaCxx/Include.mk
include ../Parallel/Include.mk
CXXFLAGS_linux+=-pthread
LDFLAGS_linux+=-pthread
