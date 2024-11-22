#include "./include/DFNGen.cuh"

int main(int argc, char *argv[])
{
    time_t t;
    time(&t);

    cuDFNsys::DFNGen myDFN;
    myDFN.RandomSeed = (unsigned long)t;

    myDFN.NumFractures = atoi(argv[1]);
    double FractureCenterLimit[6] = {0, 100, 0, 100, 0, 100};

    myDFN.FractureLocationDistributionFunctionPointer = &cuDFNsys::FractureLocationGenerationFunctionUniform; // built-in function
    myDFN.GenerateFractureLocations(FractureCenterLimit);

    double NormalRange[6] = {-1, 1, -1, 1, -1, 1};
    myDFN.FractureNormalDistributionFunctionPointer = &cuDFNsys::FractureNormalVectorGenerationFunctionUniform; // built-in function
    myDFN.GenerateFractureNormals(NormalRange);

    double FracSizeUniformRange[2] = {1, 20};
    myDFN.FractureVerticesFunctionPointer = &cuDFNsys::FractureVerticesGenerationFunctionRegularTriangleUniformRadius; // built-in function
    myDFN.GenerateFractureVertices(FracSizeUniformRange);

    // define aperture and conductivity
    // for example, b = 1e-10 * raidusOfCircumscribed ^ 0.3
    // and k =b ^ 3 / 12
    for (int i = 0; i < myDFN.NumFractures; ++i)
    {
        myDFN.Aperture[i] = pow(1e-10 * Double3Norm(myDFN.FractureVertices[myDFN.VerticesIndex[i].x] - myDFN.Location[i]), 0.3);
        myDFN.Conductivity[i] = pow(myDFN.Aperture[i], 3.) / 12.;
    }

    // for (int i = 0; i < myDFN.NumFractures; ++i)
    //     for (int j = 0; j < 3; ++j)
    //         std::cout << myDFN.FractureVertices[i * 3 + j] << (j == 2 ? "\n" : " ");
    // for (int i = 0; i < myDFN.NumFractures; ++i)
    //     for (int j = 0; j < 3; ++j)
    //         std::cout << Double3Norm(myDFN.FractureVertices[i * 3 + j] - myDFN.FractureVertices[i * 3 + (j + 1) % 3]) << (j == 2 ? "\n" : ", ");
    // for (int i = 0; i < myDFN.NumFractures; ++i)
    //     for (int j = 0; j < 3; ++j)
    //         std::cout << Double3Norm(myDFN.FractureVertices[i * 3 + j] - make_double3(0, 0, 0)) << (j == 2 ? "\n" : ", ");
    // for (int i = 0; i < myDFN.NumFractures; ++i)
    //     for (int j = 0; j < 3; ++j)
    //         std::cout << Double3Norm(myDFN.FractureVertices[i * 3 + j] - make_double3(0, 0, 0)) / Double3Norm(myDFN.FractureVertices[i * 3 + j] - myDFN.FractureVertices[i * 3 + (j + 1) % 3]) << (j == 2 ? "\n" : ", ");

    myDFN.IdentifyIntersectedFractures();

    std::string h5name = "DFNGen.h5";
    std::string xmfname = "DFNGen.xmf";
    myDFN.OutputDFNGen("DFNGen.h5");
    myDFN.OutputXMF(xmfname, h5name);
    return 0;
}