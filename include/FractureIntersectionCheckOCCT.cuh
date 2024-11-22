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
#include <iostream>
#include <cuda_runtime.h>
#include <thrust/host_vector.h>

namespace cuDFNsys
{
    bool FractureIntersectionCheckOCCT(const double3 *Fracture_1,
                                       const double3 *Fracture_2,
                                       const int &NumVertices_F_1,
                                       const int &NumVertices_F_2,
                                       const bool &IfOutputIntersection,
                                       thrust::host_vector<double3> &IntersectionEdges)
    {

        thrust::host_vector<gp_Pnt> Pnts(NumVertices_F_1);
        thrust::host_vector<TopoDS_Edge> Edges(NumVertices_F_1);
        for (int i = 0; i < NumVertices_F_1; ++i)
        {
            Pnts[i].SetX(Fracture_1[i].x),
                Pnts[i].SetY(Fracture_1[i].y),
                Pnts[i].SetZ(Fracture_1[i].z);
        }
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
        {
            Pnts_2[i].SetX(Fracture_2[i].x),
                Pnts_2[i].SetY(Fracture_2[i].y),
                Pnts_2[i].SetZ(Fracture_2[i].z);
        }
        for (int i = 0; i < NumVertices_F_2; ++i)
            Edges_2[i] = BRepBuilderAPI_MakeEdge(Pnts_2[i], Pnts_2[(i + 1) % NumVertices_F_2]);
        BRepBuilderAPI_MakeWire wireBuilder2;
        for (int i = 0; i < NumVertices_F_2; ++i)
            wireBuilder2.Add(Edges_2[i]);
        TopoDS_Wire wire2 = wireBuilder2;
        TopoDS_Face face2 = BRepBuilderAPI_MakeFace(wire2);
        
        // // Define the first polygon (triangle in 3D space)
        // // gp_Pnt p1(0, 0, 0), p2(1, 0, 0), p3(0, 1, 0);
        // gp_Pnt p1(Fracture_1[0].x, Fracture_1[0].y, Fracture_1[0].z);
        // gp_Pnt p2(Fracture_1[1].x, Fracture_1[1].y, Fracture_1[1].z);
        // gp_Pnt p3(Fracture_1[2].x, Fracture_1[2].y, Fracture_1[2].z);
        // TopoDS_Edge edge1 = BRepBuilderAPI_MakeEdge(p1, p2);
        // TopoDS_Edge edge2 = BRepBuilderAPI_MakeEdge(p2, p3);
        // TopoDS_Edge edge3 = BRepBuilderAPI_MakeEdge(p3, p1);
        // TopoDS_Wire wire1 = BRepBuilderAPI_MakeWire(edge1, edge2, edge3);
        // TopoDS_Face face1 = BRepBuilderAPI_MakeFace(wire1);
        // // Define the second polygon (triangle in 3D space)
        // // gp_Pnt q1(0, 0, 0), q2(0, 0, 1), q3(1, 0, 0);
        // gp_Pnt q1(Fracture_2[0].x, Fracture_2[0].y, Fracture_2[0].z);
        // gp_Pnt q2(Fracture_2[1].x, Fracture_2[1].y, Fracture_2[1].z);
        // gp_Pnt q3(Fracture_2[2].x, Fracture_2[2].y, Fracture_2[2].z);
        // TopoDS_Edge edge4 = BRepBuilderAPI_MakeEdge(q1, q2);
        // TopoDS_Edge edge5 = BRepBuilderAPI_MakeEdge(q2, q3);
        // TopoDS_Edge edge6 = BRepBuilderAPI_MakeEdge(q3, q1);
        // TopoDS_Wire wire2 = BRepBuilderAPI_MakeWire(edge4, edge5, edge6);
        // TopoDS_Face face2 = BRepBuilderAPI_MakeFace(wire2);

        BRepAlgoAPI_Section sectionAlgo(face1, face2);
        sectionAlgo.ComputePCurveOn1(true); // Compute 2D parameter curves on face1
        sectionAlgo.Approximation(true);    // Enable approximation if needed
        sectionAlgo.Build();

        // if (!wire1.IsNull())
        //     std::cout << "face 1 exists\n";
        // if (!wire2.IsNull())
        //     std::cout << "face 2 exists\n";
        // for (TopExp_Explorer Vex(face1, TopAbs_VERTEX); Vex.More(); Vex.Next())
        // {
        //     TopoDS_Vertex vertex = TopoDS::Vertex(Vex.Current());
        //     gp_Pnt pnt = BRep_Tool::Pnt(vertex);
        //     std::cout << "pnt face 1 " << ": X: " << pnt.X() << " - Y:" << pnt.Y() << " - Z: " << pnt.Z() << std::endl;
        // }
        // for (TopExp_Explorer Vex(face2, TopAbs_VERTEX); Vex.More(); Vex.Next())
        // {
        //     TopoDS_Vertex vertex = TopoDS::Vertex(Vex.Current());
        //     gp_Pnt pnt = BRep_Tool::Pnt(vertex);
        //     std::cout << "pnt face 2 " << ": X: " << pnt.X() << " - Y:" << pnt.Y() << " - Z: " << pnt.Z() << std::endl;
        // }

        if (sectionAlgo.IsDone())
        {
            if (IfOutputIntersection)
            {
                TopoDS_Shape intersectionShape = sectionAlgo.Shape();
                if (!intersectionShape.IsNull())
                {
                    // std::cout << "Intersection detected!" << std::endl;

                    // Explore intersection edges
                    int countOfEdge = 0;
                    for (TopExp_Explorer exp(intersectionShape, TopAbs_EDGE); exp.More(); exp.Next())
                        countOfEdge++;

                    IntersectionEdges.resize(countOfEdge * 2);

                    countOfEdge = 0;
                    for (TopExp_Explorer exp(intersectionShape, TopAbs_EDGE); exp.More(); exp.Next())
                    {
                        TopoDS_Edge intersectionEdgeII = TopoDS::Edge(exp.Current());

                        for (TopExp_Explorer ex2(intersectionEdgeII, TopAbs_VERTEX); ex2.More(); ex2.Next())
                        {
                            gp_Pnt point = BRep_Tool::Pnt(TopoDS::Vertex(ex2.Current()));
                            // double xCoord = point.X();
                            // double yCoord = point.Y();
                            // double ZCoord = point.Z();
                            IntersectionEdges[countOfEdge] = make_double3(point.X(), point.Y(), point.Z());
                            countOfEdge++;
                        }
                    }
                }
            }
            return true;
        }
        else
        {
            // std::cout << "false?\n";
            return false;
        }
    };
}