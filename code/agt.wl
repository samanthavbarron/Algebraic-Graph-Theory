(* ::Package:: *)

(*Mathematica package for computations in Algebraic Graph Theory
Developed by George Barron and Terrin Warren
Instructor: Dr. Dino Lorenzini
Course: MATH 4950
Institution: University of Georgia, Department of Mathematics
Date: Fall, 2016*)

BeginPackage[ "AGT`"]

(*Usages*)
SmithSeq::usage = "SmithSeq[g] computes the diagonal elements of the Smith Normal Form of a given graph."
LoneVerts::usage = "LoneVerts[g] computes all vertices not connected to any other vertices on a given graph."
TrimLoneVerts::usage = "TrimLoneVerts[g] removes all vertices not connected to any other vertices on a given graph."
Gprod::usage = "Gprod[G1,G2] takes the product of two graphs (with the same number of vertices) by adding edges between each consecutive pair of vertices."

Begin[ "Private`"]

SmithSeq[g_]:=SmithDecomposition[KirchhoffMatrix[g]][[2]]//Diagonal;

LoneVerts[g_]:=VertexList[g,_?(VertexDegree[g,#]<1&)];

TrimLoneVerts[g_]:=VertexDelete[g,LoneVerts[g]];

Gprod[G1_,G2_]:=
	(
		(
		ArrayFlatten[TensorProduct[{{1,0},{0,0}},G1//KirchhoffMatrix],2]+
		ArrayFlatten[TensorProduct[{{0,0},{0,1}},G2//KirchhoffMatrix],2]+
		ArrayFlatten[TensorProduct[{{0,-1},{0,0}},G1//VertexCount//IdentityMatrix],2]+
		ArrayFlatten[TensorProduct[{{0,0},{-1,0}},G1//VertexCount//IdentityMatrix],2]
		)
		+IdentityMatrix[2*VertexCount[G1]]
	)//KirchhoffGraph;

End[]

EndPackage[]



