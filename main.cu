#include "./include/DFNGen.cuh"

int main(int argc, char *argv[])
{
    time_t t;
    time(&t);

    cuDFNsys::DFNGen myDFN;
    myDFN.RandomSeed = (unsigned long)t;
    myDFN.InitRandomEngine();
    
    myDFN.NumFracturePairsCheckAtOnce = 1000000;
    myDFN.Nproc = 20;

    myDFN.NumFractures = atoi(argv[1]);
    double FractureCenterLimit[6] = {0, 1, 0, 1, 0, 1};

    double iStart = cuDFNsys::CPUSecond();
    myDFN.FractureLocationDistributionFunctionPointer = &cuDFNsys::FractureLocationGenerationFunctionUniform; // built-in function
    myDFN.GenerateFractureLocations(FractureCenterLimit);
    std::cout << "Using " << cuDFNsys::CPUSecond() - iStart << " seconds GenerateFractureLocations\n";

    iStart = cuDFNsys::CPUSecond();
    double NormalRange[6] = {-1, 1, -1, 1, -1, 1};
    myDFN.FractureNormalDistributionFunctionPointer = &cuDFNsys::FractureNormalVectorGenerationFunctionUniform; // built-in function
    myDFN.GenerateFractureNormals(NormalRange);
    std::cout << "Using " << cuDFNsys::CPUSecond() - iStart << " seconds GenerateFractureNormals\n";

    iStart = cuDFNsys::CPUSecond();
    double FracSizeUniformRange[2] = {atof(argv[2]), atof(argv[3])};
    myDFN.FractureVerticesFunctionPointer = &cuDFNsys::FractureVerticesGenerationFunctionRegularTriangleUniformRadius; // built-in function
    myDFN.GenerateFractureVertices(FracSizeUniformRange);
    std::cout << "Using " << cuDFNsys::CPUSecond() - iStart << " seconds GenerateFractureVertices\n";

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

    iStart = cuDFNsys::CPUSecond();
    myDFN.IdentifyIntersectedFractures();
    std::cout << "Using " << cuDFNsys::CPUSecond() - iStart << " seconds IdentifyIntersectedFractures\n";

    myDFN.IdentifyFractureClusters();

    myDFN.TruncateFracturesByStlDomain(std::string(argv[4]));

    std::string h5name = "DFNGen.h5";
    std::string xmfname = "DFNGen.xmf";
    myDFN.OutputDFNGen("DFNGen.h5");
    myDFN.OutputXMF(xmfname, h5name);
    return 0;
}