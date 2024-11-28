#include "./include/DFNGen.cuh"

int main(int argc, char *argv[])
{

    // double3 A[3] = {
    //     make_double3(0.6320243736608026630108270182972773909569, 0.3104544613183097268027665904810419306159, 0.3936323058926333717799650457891402766109),
    //     make_double3(0.4466573429087234980983112109242938458920, 0.5377488740561482405411197760258801281452, 0.2973462123600820183888515657599782571197),
    //     make_double3(0.3696043827417045957162144986796192824841, 0.3701313955434417657386347855208441615105, 0.5448576670161898727329230496252421289682)};
    // double3 B[3] = {
    //     make_double3(0.3143053506183439704813054049736820161343, 0.9933020353259793822076062497217208147049, 0.6503717858913916627372486800595652312040),
    //     make_double3(0.3908694856293662978075076352979522198439, 0.6010793761778031596421101312444079667330, 0.2408371493365396576180614829354453831911),
    //     make_double3(0.7390920496850732490656810114160180091858, 0.6125044191881267030552749019989278167486, 0.6947422208350285677269653206167276948690)};
    // std::cout << cuDFNsys::TrianglesIntersect(A, B) << std::endl;
    // exit(1);

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

    iStart = cuDFNsys::CPUSecond();
    myDFN.TruncateFracturesByStlDomain(std::string(argv[4]));
    std::cout << "Using " << cuDFNsys::CPUSecond() - iStart << " seconds TruncateFracturesByStlDomain\n";

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
    myDFN.IdentifyIntersectedFractures(true);
    std::cout << "Using " << cuDFNsys::CPUSecond() - iStart << " seconds IdentifyIntersectedFractures\n";

    iStart = cuDFNsys::CPUSecond();
    myDFN.IdentifyFractureClusters();
    std::cout << "Using " << cuDFNsys::CPUSecond() - iStart << " seconds IdentifyFractureClusters\n";

    std::string h5name = "DFNGen.h5";
    std::string xmfname = "DFNGen.xmf";
    myDFN.OutputDFNGen(h5name);
    myDFN.OutputXMF(xmfname, h5name);

    // iStart = cuDFNsys::CPUSecond();
    // myDFN.IdentifyIntersectedFractures(true);
    // std::cout << "Using " << cuDFNsys::CPUSecond() - iStart << " seconds IdentifyIntersectedFractures with truncated fractures\n";
    //
    // iStart = cuDFNsys::CPUSecond();
    // myDFN.IdentifyFractureClusters();
    // std::cout << "Using " << cuDFNsys::CPUSecond() - iStart << " seconds IdentifyFractureClusters with truncated fractures\n";

    // std::string h5name_II = "DFNGenII.h5";
    // std::string xmfname_II = "DFNGenII.xmf";
    // myDFN.OutputDFNGen(h5name_II);
    // myDFN.OutputXMF(xmfname_II, h5name_II);

    return 0;
}