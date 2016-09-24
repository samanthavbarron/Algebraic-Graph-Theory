(* ::Package:: *)

(*Mathematica package for computations in Algebraic Graph Theory
Developed by George Barron and Terrin Warren
Instructor: Dr. Dino Lorenzini
Course: MATH 4950
Institution: University of Georgia, Department of Mathematics
Date: Fall, 2016*)

BeginPackage[ "AGT`"]

(*Usages*)
RandomTree::usage = "RandomTree[v,e] returns a random tree with at most v vertices and e edges. This inequality is because sometimes the tree returned is not connected."
GCartProd::usage = "Constructs the Cartesian product of two given graphs."
Ones::usage = "Ones[n] creates n x n matrix of all ones."
G2prod::usage = "Gprod[G1,G2] takes the product of two graphs (with the same number of vertices) by adding edges between every vertex on one graph and every vertex on the other graph."
NonTrivialSNFQ::usage = "NonTrivialSNFQ[g] returns the negated valaue of TrivialSNFQ."
TrivialSNFQ::usage = "TrivialSNFQ[g] determines whether or not the SNF of the graph is interesting. By this, we mean that TrivialSNFQ returns true if the SNF is all ones and then an integer at the end."
SmithSeq::usage = "SmithSeq[g] computes the diagonal elements of the Smith Normal Form of a given graph."
LoneVerts::usage = "LoneVerts[g] computes all vertices not connected to any other vertices on a given graph."
TrimLoneVerts::usage = "TrimLoneVerts[g] removes all vertices not connected to any other vertices on a given graph."
Gprod::usage = "Gprod[G1,G2] takes the product of two graphs (with the same number of vertices) by adding edges between each pair of vertices from each graph."

Begin[ "Private`"]

SmithSeq[g_]:=SmithDecomposition[KirchhoffMatrix[g]][[2]]//Diagonal;

LoneVerts[g_]:=VertexList[g,_?(VertexDegree[g,#]<1&)];

TrimLoneVerts[g_]:=VertexDelete[g,LoneVerts[g]];

Ones[n_]:=ConstantArray[1,{n,n}];

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

G2prod[G1_,G2_]:=
	(
		(
		ArrayFlatten[TensorProduct[{{1,0},{0,0}},G1//KirchhoffMatrix],2]+
		ArrayFlatten[TensorProduct[{{0,0},{0,1}},G2//KirchhoffMatrix],2]+
		ArrayFlatten[TensorProduct[{{0,-1},{0,0}},G1//VertexCount//Ones],2]+
		ArrayFlatten[TensorProduct[{{0,0},{-1,0}},G1//VertexCount//Ones],2]
		)
		+IdentityMatrix[2*VertexCount[G1]](G1//VertexCount)
	)//KirchhoffGraph;

TrivialSNFQ[g_]:=
	If[
		SmithSeq[g][[-3]]==1,
		True,
		False
	];

NonTrivialSNFQ[g_]:=TrivialSNFQ[g]//Not;

GCartProd[G1_,G2_]:=
DeleteCases[
Flatten[
Table[
If[
Or[
And[u==v,EdgeQ[G1,up<->vp]],
And[EdgeQ[G2,u<->v],up==vp]
],
{u,up}<->{v,vp}
],
{u,VertexList[G1]},
{up,VertexList[G1]},
{v,VertexList[G2]},
{vp,VertexList[G2]}
],
4
],
Null
]//Graph//SimpleGraph;

RandomTree[v_, e_] := 
  FindSpanningTree@TrimLoneVerts@RandomGraph[{v, e}];

End[]

EndPackage[]









