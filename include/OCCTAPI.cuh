#pragma once
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepAlgoAPI_Section.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Edge.hxx> // Ensure this header is included
#include <gp_Pnt.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <TopExp_Explorer.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Curve.hxx>
#include <Standard.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <StlAPI_Reader.hxx>
#include <Poly_Triangulation.hxx>
#include <IntCurvesFace_ShapeIntersector.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Shell.hxx>
#include <BRepBuilderAPI_MakeShell.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRep_Builder.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <BRepTools.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <unordered_set>
#include <functional>
#include <BRepClass3d_SolidClassifier.hxx>

#include "CorrectPolygon.cuh"
#include <iostream>
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <stdexcept>

namespace cuDFNsys
{
    bool FractureIntersectionCheckOCCT(const double3 *Fracture_1,
                                       const double3 *Fracture_2,
                                       const int &NumVertices_F_1,
                                       const int &NumVertices_F_2,
                                       const bool &IfOutputIntersection,
                                       thrust::host_vector<double3> &IntersectionEdges)
    {
        // double iStart = cuDFNsys::CPUSecond();

        //--------------------face 1
        thrust::host_vector<gp_Pnt> Pnts(NumVertices_F_1);
        thrust::host_vector<TopoDS_Edge> Edges(NumVertices_F_1);
        for (int i = 0; i < NumVertices_F_1; ++i)
            Pnts[i] = gp_Pnt(Fracture_1[i].x,
                             Fracture_1[i].y,
                             Fracture_1[i].z);
        for (int i = 0; i < NumVertices_F_1; ++i)
            Edges[i] = BRepBuilderAPI_MakeEdge(Pnts[i], Pnts[(i + 1) % NumVertices_F_1]);
        BRepBuilderAPI_MakeWire wireBuilder1;
        for (int i = 0; i < NumVertices_F_1; ++i)
            wireBuilder1.Add(Edges[i]);
        TopoDS_Wire wire1 = wireBuilder1;
        TopoDS_Face face1 = BRepBuilderAPI_MakeFace(wire1);
        //--------------------face 2
        thrust::host_vector<gp_Pnt> Pnts_2(NumVertices_F_2);
        thrust::host_vector<TopoDS_Edge> Edges_2(NumVertices_F_2);
        for (int i = 0; i < NumVertices_F_2; ++i)
            Pnts_2[i] = gp_Pnt(Fracture_2[i].x,
                               Fracture_2[i].y,
                               Fracture_2[i].z);
        for (int i = 0; i < NumVertices_F_2; ++i)
            Edges_2[i] = BRepBuilderAPI_MakeEdge(Pnts_2[i], Pnts_2[(i + 1) % NumVertices_F_2]);
        BRepBuilderAPI_MakeWire wireBuilder2;
        for (int i = 0; i < NumVertices_F_2; ++i)
            wireBuilder2.Add(Edges_2[i]);
        TopoDS_Wire wire2 = wireBuilder2;
        TopoDS_Face face2 = BRepBuilderAPI_MakeFace(wire2);

        // std::cout << cuDFNsys::CPUSecond() - iStart << " seconds for loading fractures\n";

        // iStart = cuDFNsys::CPUSecond();
        Bnd_Box bbox1, bbox2;
        BRepBndLib::Add(face1, bbox1);
        BRepBndLib::Add(face2, bbox2);

        // Check intersection
        if (bbox1.IsOut(bbox2))
            return false;

        // std::cout << cuDFNsys::CPUSecond() - iStart << " seconds for BRepBndLib\n";

        if (IfOutputIntersection)
        {

            BRepAlgoAPI_Section sectionAlgo(face1, face2);
            sectionAlgo.ComputePCurveOn1(true); // Compute 2D parameter curves on face1
            sectionAlgo.Approximation(true);    // Enable approximation if needed
            sectionAlgo.Build();

            if (sectionAlgo.IsDone())
            {

                TopoDS_Shape intersectionShape = sectionAlgo.Shape();
                if (!intersectionShape.IsNull())
                {
                    // if (!IfOutputIntersection)
                    //     return true;

                    int countOfEdge = 0;
                    for (TopExp_Explorer exp(intersectionShape, TopAbs_EDGE); exp.More(); exp.Next())
                        countOfEdge++;

                    IntersectionEdges.resize(countOfEdge * 2);
                    // std::cout << countOfEdge << "\n";
                    if (countOfEdge == 0)
                        return false;

                    countOfEdge = 0;

                    for (TopExp_Explorer exp(intersectionShape, TopAbs_EDGE); exp.More(); exp.Next())
                    {
                        TopoDS_Edge intersectionEdgeII = TopoDS::Edge(exp.Current());

                        for (TopExp_Explorer ex2(intersectionEdgeII, TopAbs_VERTEX); ex2.More(); ex2.Next())
                        {
                            gp_Pnt point = BRep_Tool::Pnt(TopoDS::Vertex(ex2.Current()));
                            IntersectionEdges[countOfEdge] = make_double3(point.X(), point.Y(), point.Z());
                            countOfEdge++;
                        }
                    }

                    return true;
                }
                else
                    return false;
            }
            else
            {
                throw std::runtime_error("In `cuDFNsys::FractureIntersectionCheckOCCT`, I cannot determine if two fractures intersect I\n");
            }
        }
        else
        {
            // iStart = cuDFNsys::CPUSecond();

            BRepExtrema_DistShapeShape distAlgo(face1, face2);
            distAlgo.Perform();

            if (!distAlgo.IsDone())
            {
                throw std::runtime_error("In `cuDFNsys::FractureIntersectionCheckOCCT`, I cannot determine if two fractures intersect II\n");
            }

            // std::cout << cuDFNsys::CPUSecond() - iStart << " seconds for BRepExtrema_DistShapeShape\n";

            // Check if distance is within tolerance
            return distAlgo.Value() <= 1e-10;
        }
    };

    TopoDS_Solid LoadAStlDomain(const std::string &stlname, const bool &IfReverse = true)
    {
        // Load STL file as faces
        StlAPI_Reader reader;
        TopoDS_Shape shape;
        reader.Read(shape, stlname.c_str());

        if (shape.IsNull())
            throw std::runtime_error("Error: In cuDFNsy::LoadAStlDomain, failed to read STL file or invalid shape\n");
        // return TopoDS_Solid(); // Return an empty solid on error

        if (shape.ShapeType() != TopAbs_SHELL)
            throw std::runtime_error("Error: In cuDFNsy::LoadAStlDomain, rightnow only an enclosed domain is allowed\n");

        // Step 2: Create a shell from the shape
        BRepBuilderAPI_Sewing sewing(1.0e-6); // Sewing tolerance
        for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next())
        {
            TopoDS_Face face = TopoDS::Face(exp.Current());
            sewing.Add(face);
        }
        sewing.Perform();

        TopoDS_Shell shell = TopoDS::Shell(sewing.SewedShape());

        // Step 3: Verify if the shell is closed
        BRepCheck_Analyzer checker(shell);
        if (!checker.IsValid())
        {
            throw std::runtime_error("Error: In cuDFNsy::LoadAStlDomain, The shell is not valid or closed. Cannot create a solid.\n");
        }

        // TopAbs_Orientation orientation = shell.Orientation();
        if (IfReverse)
        {
            shell.Orientation(TopAbs_REVERSED); // Adjust if necessary
        }

        // Step 4: Create a solid from the shell
        BRepBuilderAPI_MakeSolid makeSolid;
        makeSolid.Add(shell);

        if (!makeSolid.IsDone())
        {
            throw std::runtime_error("Error: In cuDFNsy::LoadAStlDomain, Failed to create solid from shell.\n");
        }

        TopoDS_Solid solid = makeSolid.Solid();

        Bnd_Box bbox;
        BRepBndLib::Add(solid, bbox);
        bool isInfinite = false;
        if (bbox.IsOpenXmin() == true || bbox.IsOpenXmax() == true ||
            bbox.IsOpenYmin() == true || bbox.IsOpenYmax() == true ||
            bbox.IsOpenZmin() == true || bbox.IsOpenZmax() == true)
            isInfinite = true;
        if (isInfinite)
        {
            throw std::runtime_error("Error: In cuDFNsy::TruncatePolygon, the solid is infinite. Check the STL file or topology.\n");
        }
        // Standard_Real Xmin, Ymin, Zmin, Xmax, Ymax, Zmax;
        // bbox.Get(Xmin, Ymin, Zmin, Xmax, Ymax, Zmax);
        // std::cout << Xmin << ", " << Ymin << ", " << Zmin << ", " << Xmax << ", " << Ymax << ", " << Zmax << std::endl;
        return solid;
    }

    bool IsPointInsideSolid(const TopoDS_Solid &solid, const gp_Pnt &point, double tolerance = 1.0e-7)
    {
        // Create a solid classifier
        BRepClass3d_SolidClassifier solidClassifier(solid);

        // Perform classification with point and tolerance
        solidClassifier.Perform(point, tolerance);

        // Check the result of the classification
        TopAbs_State state = solidClassifier.State();

        // Return true if the point is inside the solid
        return (state == TopAbs_IN);
    }

    TopoDS_Face FractureToTopoDS_Face(const double3 *Vertices, const int &NumVertices)
    {
        thrust::host_vector<gp_Pnt> Pnts(NumVertices);
        thrust::host_vector<TopoDS_Edge> Edges(NumVertices);
        for (int i = 0; i < NumVertices; ++i)
            Pnts[i] = gp_Pnt(Vertices[i].x,
                             Vertices[i].y,
                             Vertices[i].z);
        for (int i = 0; i < NumVertices; ++i)
            Edges[i] = BRepBuilderAPI_MakeEdge(Pnts[i], Pnts[(i + 1) % NumVertices]);
        BRepBuilderAPI_MakeWire wireBuilder1;
        for (int i = 0; i < NumVertices; ++i)
            wireBuilder1.Add(Edges[i]);
        TopoDS_Wire wire1 = wireBuilder1;
        TopoDS_Face face1 = BRepBuilderAPI_MakeFace(wire1);

        TopAbs_Orientation orientation = face1.Orientation();
        if (orientation == TopAbs_FORWARD)
        {
            face1.Orientation(TopAbs_REVERSED); // Adjust if necessary
        }

        return face1;
    };

    thrust::host_vector<double3> TruncatePolygon(const TopoDS_Face &polygon, const TopoDS_Solid &DomainSolid)
    {

        // for (TopExp_Explorer exp(polygon, TopAbs_EDGE); exp.More(); exp.Next())
        // {
        //     TopoDS_Edge intersectedFace = TopoDS::Edge(exp.Current());
        //     // Extract vertices from the face's outer wire
        //     for (TopExp_Explorer exWire(intersectedFace, TopAbs_VERTEX); exWire.More(); exWire.Next())
        //     {
        //         gp_Pnt point = BRep_Tool::Pnt(TopoDS::Vertex(exWire.Current()));
        //         std::cout << make_double3(point.X(), point.Y(), point.Z()) << std::endl;
        //         break;
        //     }
        // }
        BRepAlgoAPI_Common commonOp(DomainSolid, polygon);
        commonOp.Build();

        thrust::host_vector<double3> NewPolygon;

        if (commonOp.IsDone())
        {
            TopoDS_Shape resultShape = commonOp.Shape();
            if (resultShape.IsNull())
            {
                return NewPolygon; // Return empty if no intersection
            }

            int countEdge = 0;
            for (TopExp_Explorer exp(resultShape, TopAbs_EDGE); exp.More(); exp.Next())
                countEdge++;
            NewPolygon.reserve(countEdge * 3);

            for (TopExp_Explorer exp(resultShape, TopAbs_VERTEX); exp.More(); exp.Next())
            {
                gp_Pnt point = BRep_Tool::Pnt(TopoDS::Vertex(exp.Current()));
                NewPolygon.push_back(make_double3(point.X(), point.Y(), point.Z()));
            }
            NewPolygon.shrink_to_fit();
            NewPolygon = CorrectPolygon(NewPolygon, 1e-3);
        }
        else
        {
            throw std::runtime_error("Error: Failed to compute intersection between polygon and solid.");
        }

        return NewPolygon;

        /// if (1 == 2)
        /// {
        ///     BRepAlgoAPI_Section sectionOp(DomainSolid, polygon);
        ///     // sectionOp.ComputePCurveOn1(true);
        ///     sectionOp.Build();
        ///
        ///     thrust::host_vector<double3> NewPolygon;
        ///
        ///     if (sectionOp.IsDone())
        ///     {
        ///         TopoDS_Shape intersectionShape = sectionOp.Shape();
        ///         if (!intersectionShape.IsNull())
        ///         {
        ///             // if (!IfOutputIntersection)
        ///             //     return true;
        ///
        ///             int countOfEdge = 0;
        ///             for (TopExp_Explorer exp(intersectionShape, TopAbs_EDGE); exp.More(); exp.Next())
        ///                 countOfEdge++;
        ///
        ///             NewPolygon.resize(countOfEdge * 2);
        ///             // std::cout << countOfEdge << "\n";
        ///             if (countOfEdge == 0)
        ///             {
        ///                 throw std::runtime_error("In `cuDFNsys::TruncatePolygon`, countOfEdge = 0, and shapetype = " + std::to_string(intersectionShape.ShapeType()) + "\n");
        ///             }
        ///
        ///             countOfEdge = 0;
        ///
        ///             for (TopExp_Explorer exp(intersectionShape, TopAbs_EDGE); exp.More(); exp.Next())
        ///             {
        ///                 TopoDS_Edge intersectionEdgeII = TopoDS::Edge(exp.Current());
        ///
        ///                 for (TopExp_Explorer ex2(intersectionEdgeII, TopAbs_VERTEX); ex2.More(); ex2.Next())
        ///                 {
        ///                     gp_Pnt point = BRep_Tool::Pnt(TopoDS::Vertex(ex2.Current()));
        ///                     NewPolygon[countOfEdge] = make_double3(point.X(), point.Y(), point.Z());
        ///                     countOfEdge++;
        ///                     // break;
        ///                 }
        ///             }
        ///
        ///             return NewPolygon;
        ///         }
        ///         else
        ///             return NewPolygon;
        ///     }
        ///     else
        ///         throw std::runtime_error("In `cuDFNsys::TruncatePolygon`, I cannot determine the intersection between polygon and domain stl\n");
        /// };
    };
}