include Path.mk

#MKLROOT=/opt/intel/oneapi/mkl/2025.0
Mesh_INCDIRS=-I $(Hdf5IncludePath) \
		-I $(GmshIncludePath) \
		-I $(EigenIncludePath) \
		-I /usr/local/include \
        -I $(UmfpackIncludePath) \
		-I ./include -I/usr/include/opencascade

Mesh_Lib_= -Xcompiler=-fopenmp \
	-lpthread \
	-lm \
	-lcudadevrt \
	-lcudart_static \
	-lrt \
	-lpthread \
	-ldl \
	-lgmsh \
	-lumfpack \
	-lamd \
	-lhdf5_cpp \
	-lhdf5 \
	-lsz \
	-lz \
	-ldl \
	-L$(GmshLibraryPath) \
	-L$(Hdf5LibraryPath) \
 	-L$(UmfpackLibraryPath) \
	-L/usr/lib/x86_64-linux-gnu \
	-lTKernel -lTKBO -lTKBRep -lTKGeomBase -lTKG2d -lTKG3d -lTKMath -lTKBool -lTKTopAlgo -lTKSTL

all: main

.PHONY: all main

main: ./main.cu
	$(NVCC) -DUSE_DOUBLES ./main.cu -o ./main $(Mesh_INCDIRS) $(Mesh_Lib_) \
	-Xcudafe --display_error_number -Xcudafe --diag_suppress=3057 \
	-Xcudafe --diag_suppress=1301 \
	-Xcudafe --diag_suppress=3059 \
	-Xcudafe --diag_suppress=3060 \
	-Xcudafe --diag_suppress=3058 \
	-Xcudafe --diag_suppress=3056 \
	-arch=sm_60 -std=c++17 -rdc=true --extended-lambda


clean:
	rm -rf main