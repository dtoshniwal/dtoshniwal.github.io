-- Copyright 2014-2015: Mike Dipasquale
-- You may redistribute this file under the terms of the GNU General Public
-- License as published by the Free Software Foundation, either version 2
-- of the License, or any later version.

------------------------------------------
------------------------------------------
-- Header
------------------------------------------
------------------------------------------

if version#"VERSION" <= "1.4" then (
        needsPackage "FourierMotzkin"
        )

newPackage select((
    "SimplicialMDSplines",
    Version => "0.1.0", 
    Date => "December 3, 2018",
    Authors => {
        {Name => "Mike DiPasquale", Email => "midipasq@gmail.com", HomePage => "http://illinois.edu/~dipasqu1"},
        {Name => "Gwyn Whieldon", Email => "whieldon@hood.edu", HomePage => "http://cs.hood.edu/~whieldon"},
        {Name => "Eliana Duarte", Email => "emduart2@illinois.edu", HomePage => "http://illinois.edu/~emduart2"},
        {Name => "Daniel Irving Bernstein", Email=> "dibernst@ncsu.edu", HomePage =>"http://www4.ncsu.edu/~dibernst"},
        {Name => "Deepesh Toshniwal", Email=> "deepesh@utexas.edu", HomePage =>"http://dtoshniwal.github.io"}
    },
    Headline => "Derived from the package `AlgebraicSplines', this package can be used for computing the dimension of polynomial splines of mixed degrees and mixed orders of smoothness on simplicial complexes.",
    Configuration => {},
    DebuggingMode => true,
    if version#"VERSION" > "1.4" then PackageExports => {
        "FourierMotzkin"
    }
    ), x -> x =!= null
)

if version#"VERSION" <= "1.4" then (
        needsPackage "FourierMotzkin"
)

export {
    "BaseRing",
    "splineMatrix",
    "splineModule",
    "RingType",
    "Homogenize",
    "VariableName",
    "splineDimensionTable",
    "homologyDimensionTable",
    "postulationNumber",
    "hilbertComparisonTable",
    "cellularComplex",
    "idealsComplex",
    "splineComplex"
    --"formsList",
    --"getAllDegreeDeficits",
    --"simpBoundary",
    --"boundaryComplex",
    --"getCodim1Intersections",
    --"getCodimDFacesSimplicial",
    --"codim1Cont",
    --"getInteriorFaces",
    --"subsetL",
    --"arrangeSmoothnessIndices"
}


------------------------------------------
------------------------------------------
-- Data Types and Constructors
------------------------------------------
------------------------------------------


------------------------------------------
------------------------------------------
-- Methods
------------------------------------------
------------------------------------------


------------------------------------------
subsetL=method()
------------------------------------------

--Containment function for lists--
subsetL(List,List):=Boolean=>(L1,L2)->(
    all(L1,f->member(f,L2))
)

------------------------------------------


------------------------------------------
getInteriorFaces = method()
------------------------------------------
------------------------------------------
--Inputs: 
------------------------------------------
--F = list of facets of a simplicial 
--    complex
------------------------------------------
--Outputs:
-----------------------------------------
--intC = interior cells of decreasing dimension
-----------------------------------------

getInteriorFaces(List) := List => (F) ->(
    d := #(first F)-1;
    boundaryF := boundaryComplex(F);
    C := apply(d+1, i-> getCodimDFacesSimplicial(F,i));
    boundaryC := join({{}},apply(d, i-> getCodimDFacesSimplicial(boundaryF,i)));
    intC := apply(#C, i -> select(C_i, f -> not member(f,boundaryC_i)))
)

------------------------------------------


------------------------------------------
getAllDegreeDeficits = method()
------------------------------------------
------------------------------------------
--Inputs: 
------------------------------------------
--F = list of facets of a simplicial
--    complex
--degD = list/value of degree deficits on faces
------------------------------------------
--Outputs:
-----------------------------------------
--degDAll = degree deficits on cells of the complex
-----------------------------------------

getAllDegreeDeficits(List, ZZ) := List => (F, degD) ->(
    degD = apply(#F, i -> degD);
    degDAll := getAllDegreeDeficits(F, degD)
)

getAllDegreeDeficits(List, List) := List => (F, degD) ->(
    intC := getInteriorFaces(F);
    degDAll := {degD};
    scan(#intC-1,i->(
        if #(intC_(i+1))==0 then(
            newdef := null;
        )else(
            facets := apply(#(intC_(i+1)), e-> positions(intC_0, f-> all((intC_(i+1))_e,v->member(v,f))));
            newdef = apply(#(intC_(i+1)), e-> min apply(facets_e, j-> degD_j));
        );
    degDAll = append(degDAll,newdef);
    ));
    degDAll
)

------------------------------------------


------------------------------------------
getCodim1Intersections = method()
------------------------------------------
------------------------------------------
--Inputs: 
------------------------------------------
--F = list of facets of a simplicial
--    complex
------------------------------------------
--Outputs:
-----------------------------------------
--E = largest (under containment) intersections
--of facets.  If input is hereditary, this is
--the list of interior codim 1 intersections
-----------------------------------------

getCodim1Intersections(List) := List => F ->(
    n := #F;
    d := #(F_0);
    --For each non-final facet, construct all codimension 1 subsets.
    codim1faces := apply(n-1, i -> subsets(F_i,d-1));
    --Check if a codimension 1 subset is contained in another facet,
    --store it as a codim 1 intersection.
    codim1int:=sort flatten apply(#codim1faces, i -> 
        select(codim1faces_i, 
            s -> any(F_{i+1..n-1}, 
                f-> subsetL(s,f))));
    codim1int
)

------------------------------------------


------------------------------------------
getCodimDFacesSimplicial = method()
------------------------------------------
------------------------------------------
--Inputs: 
------------------------------------------
--F = list of facets of a simplicial complex
--d = desired codimension
------------------------------------------
--Outputs:
-----------------------------------------
--E = list of (all) codim d faces
-----------------------------------------
getCodimDFacesSimplicial(List,ZZ) := List => (F,D) -> (
    d := #first(F);
    unique flatten apply(F, f-> subsets(f,d-D))
)

-----------------------------------------


------------------------------------------
arrangeSmoothnessIndices = method()
------------------------------------------

------------------------------------------
--Inputs: 
------------------------------------------
--F = Faces of the simplicial complex
--E = Interior edges of the simplicial complex
--r = smoothness indices
------------------------------------------
--Outputs:
-----------------------------------------
--r = rearranged smoothness indices
--according to the internal numbering
-----------------------------------------
arrangeSmoothnessIndices(List,List,List) := List => (F,E,r) ->(
    intC := getInteriorFaces(F);
    r' := new MutableList from r;
    i := 0; while i<#(intC_1) do (
                j := 0; while j<#E do (
                    if subsetL(E_j, (intC_1)_i) then (r'#i = r_j;);
                    j = j + 1;);
                i = i + 1;);
    r' = new List from r'
)

-----------------------------------------


-----------------------------------------
formsList=method(Options=>{
    symbol Homogenize => true, 
    symbol VariableName => getSymbol "t",
    symbol CoefficientRing => QQ,
    symbol BaseRing => null})
----------------------------------------------------
--This method returns a list of forms corresponding to codimension one faces
--------
--Input:
--V = vertex list
--E = cod 1 face list
--r = smoothness parameter
--degD = degree deficits corresponding to the edges
--------
--Output:
--List of forms defining input list of codimension one faces 
--raised to (r+1) power
------------------------------------------------------------

formsList(List,List,ZZ,ZZ):= List => opts->(V,E,r,degD)->(
    r = apply(#E, i -> r);
    degD = apply(#E, i -> degD);
    fl := formsList(V, E, r, degD, opts);
    fl
)

formsList(List,List,ZZ,List):= List => opts->(V,E,r,degD)->(
    r = apply(#E, i -> r);
    fl := formsList(V, E, r, degD, opts);
    fl
)

formsList(List,List,List,ZZ):= List => opts->(V,E,r,degD)->(
    degD = apply(#E, i -> degD);
    fl := formsList(V, E, r, degD, opts);
    fl
)

formsList(List,List,List,List):= List => opts->(V,E,r,degD)->(
    --To homogenize, we append a 1 as the final coordinate of each vertex coord in list V.
    --If not homogenizing, still need 1s for computing equations
    d := #(first V);
    V = apply(V, v-> append(v,1));
    if opts.BaseRing === null then (
        S := createSplineRing(d,opts);
    )else(
        S = opts.BaseRing;
    );
    if opts.Homogenize then (
        varlist := vars S;
    ) else (
        if any(degD, deficit -> deficit != 0) then (
            error "Degree deficits are not identically 0; Homogenize must be set to true."
        );
        varlist = (vars S)|(matrix {{sub(1,S)}});
    );
    varCol := transpose varlist;
    hvar := varCol_(d, 0);
    M := (transpose(matrix(S,V)));
    mM := numrows M;
    minorList := apply(E, e-> gens gb minors(mM,matrix(M_e)|varCol));
    if any(minorList, I-> ideal I === ideal 1) then (
        error "Some vertices on entered face are not in codimension 1 face."
    );
    flatten apply(#E, i -> hvar^(degD_i) * (minorList_i_(0,0))^(r_i + 1))
)

-----------------------------------------


-----------------------------------------
createSplineRing = method(Options => {
    symbol Homogenize => true, 
    symbol VariableName => getSymbol "t",
    symbol CoefficientRing => QQ,
    symbol BaseRing => null})
----------------------------------------------------

createSplineRing(ZZ):= PolynomialRing => opts -> (d) -> (
    t := opts.VariableName;
    if opts.BaseRing === null then (
        if opts.Homogenize then (
            S := (opts.CoefficientRing)[t_0..t_d];
        ) else (
            S = (opts.CoefficientRing)[t_1..t_d];
        );
    ) else (
        S = opts.BaseRing
    );
    S
)

-----------------------------------------


-----------------------------------------
splineMatrix = method(Options => {
    symbol Homogenize => true, 
    symbol VariableName => getSymbol "t",
    symbol CoefficientRing => QQ,
    symbol BaseRing => null})
------------------------------------------

------------------------------------------
--Inputs: 
------------------------------------------
--V = list of coordinates of vertices
--F = ordered lists of facets
--r = smoothness order
--degD = degree deficits on faces
------------------------------------------
--Outputs:
-- BM = matrix with columns corresponding
-- to facets and linear forms separating facets.
------------------------------------------
splineMatrix(List,List,ZZ,ZZ) := Matrix => opts -> (V,F,r,degD) -> (
    degD = apply(#F, i->degD);
    intC := getInteriorFaces(F);
    E := intC_1;
    r = apply(#E, i->r);
    splineM := splineMatrix(V,F,{E,r},degD,opts);
    splineM
)

splineMatrix(List,List,ZZ,List) := Matrix => opts -> (V,F,r,degD) -> (
    intC := getInteriorFaces(F);
    E := intC_1;
    r = apply(#E, i->r);
    splineM := splineMatrix(V,F,{E,r},degD,opts);
    splineM
)

------------------------------------------
--Inputs: 
------------------------------------------
--V = list of coordinates of vertices
--F = ordered lists of facets
--Er = {list of interior edges, list/value of smoothness}
--degD = degree deficits on faces
------------------------------------------
--Outputs:
-- BM = matrix with columns corresponding
-- to facets and linear forms separating facets.
------------------------------------------
splineMatrix(List,List,List,ZZ) := Matrix => opts -> (V,F,Er,degD) -> (
    degD = apply(#F, i->degD);
    splineM := splineMatrix(V,F,Er,degD,opts);
    splineM
)

splineMatrix(List,List,List,List) := Matrix => opts -> (V,F,Er,degD) -> (
    E := Er_0;
    r := Er_1;
    try #E==#r then () else ( r = apply(#E, i -> r); );
    d := # (first V);
    --Compute which facets are adjacent to each edge:
    facetEdgeH := apply(#E, e-> positions(F, f-> all(E_e,v-> member(v,f))));
    --List of forms definining interior codim one faces (raised to (r+1) power)
    flist := formsList(V,E,r,0,opts);
    T := diagonalMatrix(flist);
    varlist := vars ring flist_(0);
    hvar := varlist_(0,d);
    --Compute top boundary map for complex:
    BM := matrix apply(
        facetEdgeH, i-> apply(
            #F, j-> if (j === first i) then 1*hvar^(degD_j) else if (j===last i) then -1*hvar^(degD_j) else 0));
    splineM := BM|T;
    splineM
)

------------------------------------------


------------------------------------------
splineModule = method(Options => {
    symbol Homogenize => true, 
    symbol VariableName => getSymbol "t",
    symbol CoefficientRing => QQ,
    symbol BaseRing => null})
------------------------------------------
-- This method computes the splineModule
-- of a complex Delta, given by either
-- facets, codim 1 faces, and vertex coors,
-- or by pairs of adjacent faces and
-- linear forms.
------------------------------------------

------------------------------------------
--Inputs: 
------------------------------------------
--V = list of coordinates of vertices
--F = ordered lists of facets
--r = smoothness order
--degD = degree deficits on faces
------------------------------------------
--Outputs:
--Spline module S^r(Delta)
------------------------------------------
splineModule(List,List,ZZ,ZZ) := Matrix => opts -> (V,F,r,degD) -> (
    degD = apply(#F, i->degD);
    intC := getInteriorFaces(F);
    E := intC_1;
    r = apply(#E, i->r);
    M := splineModule(V,F,{E,r},degD,opts);
    M
)

splineModule(List,List,ZZ,ZZ) := Matrix => opts -> (V,F,r,degD) -> (
    intC := getInteriorFaces(F);
    E := intC_1;
    r = apply(#E, i->r);
    M := splineModule(V,F,{E,r},degD,opts);
    M
)

------------------------------------------
--Inputs: 
------------------------------------------
--V = list of coordinates of vertices
--F = ordered lists of facets
--Er = {list of edges, list/value of smoothness}
--degD = degree deficits on faces
------------------------------------------
--Outputs:
--Spline module S^r(Delta)
------------------------------------------
splineModule(List,List,List,ZZ) := Matrix => opts -> (V,F,Er,degD) -> (
    degD = apply(#F, i->degD);
    M := splineModule(V,F,Er,degD,opts);
    M
)

splineModule(List,List,List,List) := Matrix => opts -> (V,F,Er,degD) -> (
    AD := splineMatrix(V,F,Er,degD,opts);
    K := ker AD;
    b := #F;
    image submatrix(gens K, toList(0..b-1),)
)

------------------------------------------


-------------------------------------------
splineDimensionTable=method()
-------------------------------------------

-------------------------------------------
-----Inputs:
-------------------------------------------
-- a= lower bound of dim table
-- b= upper bound of dim table
-- M= module
--------------------------------------------
--- Outputs:
--------------------------------------------
-- A net with the degrees between a and b on top row
-- and corresponding dimensions of graded pieces
-- of M in bottom row
-------------------------------------------
splineDimensionTable(ZZ,ZZ,Module):=Net=> opts -> (a,b,M)->(
    r1:=prepend("Degree",toList(a..b));
    r2:=prepend("Dimension",apply(toList(a..b),i->hilbertFunction(i,M)));
    netList {r1,r2}
)

-------------------------------------------
-----Inputs:
-------------------------------------------
-- a= lower bound of range
-- b= upper bound of range
-- L= list {V,F,Er,degD} or {V,F,r,degD} where
-- V is a list of vertices, F a list of facets
-- Er is {Edges, smoothness}, r is smoothness,
-- degD is degree deficits on facets
-------------------------------------------
-----Outputs:
-------------------------------------------
-- A table with the dimensions of the graded pieces
-- of the spline module in the range (a,b)
-------------------------------------------
splineDimensionTable(ZZ,ZZ,List):= Net=> opts -> (a,b,L)->(
    M := splineModule(L_0,L_1,L_2,L_3,opts);
    splineDimensionTable(a,b,M)
)

-------------------------------------------


-------------------------------------------
homologyDimensionTable=method()
-------------------------------------------

-------------------------------------------
-----Inputs:
-------------------------------------------
-- a= lower bound of range
-- b= upper bound of range
-- L= list {V,F,Er,degD} or {V,F,r,degD} where
-- V is a list of vertices, F a list of facets
-- Er is {Edges, smoothness}, r is smoothness,
-- degD is degree deficits on facets
-------------------------------------------
-----Outputs:
-------------------------------------------
-- A table with the dimensions of the graded pieces
-- of all homology modules from the mixed
-- degree Schenck-Stillman complexes
-------------------------------------------
homologyDimensionTable(ZZ,ZZ,List):= (a,b,L) ->(
    M := splineModule(L_0, L_1, L_2, L_3);
    I := idealsComplex(L_0, L_1, L_2, L_3);
    C := cellularComplex(L_1, L_3);
    Q := splineComplex(L_0, L_1, L_2, L_3);
    HI := homology(I);
    HC := homology(C);
    HQ := homology(Q);
    r1 :=prepend("Degree",toList(a..b));
    r2 :=prepend("H_1(I) (idealsComplex) ",apply(toList(a..b),i->hilbertFunction(i,HI_1)));
    r3 :=prepend("H_0(I) (-------------) ",apply(toList(a..b),i->hilbertFunction(i,HI_0)));
    r4 :=prepend("H_2(C) (cellularComplex) ",apply(toList(a..b),i->hilbertFunction(i,HC_2)));
    r5 :=prepend("H_1(C) (---------------) ",apply(toList(a..b),i->hilbertFunction(i,HC_1)));
    r6 :=prepend("H_0(C) (---------------) ",apply(toList(a..b),i->hilbertFunction(i,HC_0)));
    r7 :=prepend("H_2(Q) (splineComplex) ",apply(toList(a..b),i->hilbertFunction(i,HQ_2)));
    r8 :=prepend("H_1(Q) (-------------) ",apply(toList(a..b),i->hilbertFunction(i,HQ_1)));
    r9 :=prepend("H_0(Q) (-------------) ",apply(toList(a..b),i->hilbertFunction(i,HQ_0)));
    r10 :=prepend("dim(S) (splineComplex) ",apply(toList(a..b),i->hilbertFunction(i,HQ_2)));
    r11 :=prepend("Eul. ch. (-----------) ",apply(toList(a..b),i->hilbertFunction(i,HQ_2)-hilbertFunction(i,HQ_1)+hilbertFunction(i,HQ_0)));
    r12 :=prepend("dim(M) (splineModule) ",apply(toList(a..b),i->hilbertFunction(i,M)));
    netList {r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12}
)

------------------------------------------


-------------------------------------------
postulationNumber=method()
-------------------------------------------

-------------------------------------------
-----Inputs:
-------------------------------------------
-- M, a graded module
--------------------------------------------
--- Outputs:
--------------------------------------------
--- The postulation number (largest integer 
--- for which Hilbert function and polynomial 
--- of M disagree).
--------------------------------------------
postulationNumber(Module):= (N) ->(
    k := regularity N;
    h := hilbertPolynomial(N);
    while hilbertFunction(k,N)==h(k) do (k=k-1);
    k
)

------------------------------------------


-----------------------------------------
hilbertComparisonTable=method()
-------------------------------------------

-------------------------------------------
-----Inputs:
-------------------------------------------
-- a= an integer, lower bound
-- b= an integer, upper bound
-- M= graded module over a polynomial ring
--------------------------------------------
--- Outputs:
--------------------------------------------
--- A table whose top two rows are the same as
--- the output of splineDimensionTable and whose 
--- third row compares the first two to the
--- Hilbert Polynomial
--------------------------------------------
hilbertComparisonTable(ZZ,ZZ,Module):= (a,b,M) ->(
    h := hilbertPolynomial(M);
    r1 := prepend("Degree",toList(a..b));
    r2 := prepend("Dimension",apply(toList(a..b),i->hilbertFunction(i,M)));
    r3 := prepend("HilbertPoly",apply(toList(a..b),i->h(i)));
    netList {r1,r2,r3}
)

---------------------------------------------


------------------------------------------
simpBoundary = method()
------------------------------------------

------------------------------------------
--Input:
--F = list of codim i faces
--E = list of codim i+1 faces
------------------------------------------
--Output:
--B = boundary map matrix between codim i and codim i+1 faces
------------------------------------------
simpBoundary(List,List) := Matrix => (F,E) -> (
    F = apply(F, f-> sort f);
    E = apply(E, e-> sort e);
    if #F==0 then return(
        map(image matrix{{0}},image matrix{{0}},matrix{{0}})
    )else(
        if #E==0 then return(
            map(image matrix{{0}},ZZ^(length F),matrix{apply(length F,i->0)})
        )else(
            tempLF := {};
            rowList := {};
            apply(F, f-> (
                tempLF = hashTable apply(#f, v-> position(E,e-> e == drop(f,{v,v})) => (-1)^v);
                rowList = append(rowList,apply(#E, j->if member(j,keys tempLF) then tempLF#j else 0));
            ));
            return transpose matrix rowList
        )
    )
)

------------------------------------------


------------------------------------------
boundaryComplex = method()
------------------------------------------

------------------------------------------
--Input:
--F= list of facets of a simplicial complex
----which is a pseudomanifold (Important!)
------------------------------------------
--Output:
--A list of codim one faces on the boundary
------------------------------------------
boundaryComplex(List) := List => F -> (
    n := #F;
    d := #(F_0);
    codim1faces := unique flatten apply(n,i-> subsets(F_i,d-1));
    select(codim1faces, f-> number(F, g-> all(f, v-> member(v,g))) === 1)
)

------------------------------------------
--Input:
--V= vertex list
--F= list of facets of a complex (polyhedral or simplicial)
------------------------------------------------
--Output:
--A list of codim one faces on the boundary
------------------------------------------------
boundaryComplex(List,List):= List => (V,F)-> (
    candidateList := flatten apply(F,f->facetsN(V,f));
    select(candidateList, f->#(positions(candidateList,l->l==f))==1)
)

------------------------------------------


------------------------------------------
facetsN = method()
------------------------------------------------

------------------------------------------------
--Input:
--V=vertex list
--(Optional) F=list of indices of V to take for convex hull
------------------------------------------------
--Output: List of facets of convex hull of V_F
------------------------------------------------
facetsN(List,List):= List => (V,L)->(
    V = apply(V_L,v->prepend(1,v));
    --get minimal list of inequalities defining the conical hull of V
    H := fourierMotzkin transpose matrix V;
    halfspaces := H#0;
    --apply inequalities to vertex list to see which vertices are contained in faces
    M := (transpose halfspaces)*(transpose matrix V); 
    facetList := apply(numrows M, i->select(numcols M,j-> M_(i,j)==0));
    --label facet indices by elements of the list L
    apply(facetList,l->L_l)
)

facetsN(List):= List => V ->(
    n := #V;
    F := toList(0..n-1);
    facetsN(V,F)
)
------------------------------------------------


------------------------------------------------
cellularComplex = method(Options =>{
    symbol Homogenize => true, 
    symbol VariableName => getSymbol "t",
    symbol CoefficientRing => QQ,
    symbol BaseRing => null})
------------------------------------------------
---This method computes the cellular chain complex
---of a polyhedral complex with coefficients in a polynomial
---ring, modulo the boundary
------------------------------------------------

------------------------------------------------
---Inputs: A list of facets
------------------------------------------------
---Outputs: The cellular chain complex whose homology
--- is the homology of the simplicial complex relative
--- to its boundary.
--------------------------------------------------
cellularComplex(List,ZZ) := ChainComplex => opts -> (F,degD) -> (
    degD = apply(#F, i->degD);
    cellularComplex(F, degD)
)

cellularComplex(List,List) := ChainComplex => opts -> (F,degD) -> (
    d := #(first F) -1;
    S := createSplineRing(d,opts);
    --list of interior faces in order of increasing codimension--
    intC := getInteriorFaces(F);
    --get degree deficits on all simplices of the complex
    degDAll := getAllDegreeDeficits(F, degD);
    --list of modules which will define chain complex--
    varlist := vars S;
    hvar := varlist_(0,d);
    fullmodulelist := { directSum apply(#F, i-> module ideal ( hvar^((degDAll_0)_i) ) ) };
    scan(#intC-1,i->(
        if #(intC_(i+1))==0 then(
            newMod := image matrix{{0_S}};
        )else(
            newMod = directSum apply(#(intC_(i+1)),j->(
                e := (intC_(i+1))_j;
                CE := codim1Cont(intC_1,e);
                return module ideal (hvar^((degDAll_(i+1))_j))
            ));
        );
        fullmodulelist = append(fullmodulelist,newMod)
    ));
    --defining the chain complex
    CCSS :=chainComplex(reverse apply(#intC-1, c-> (
        inducedMap(fullmodulelist_(c+1),fullmodulelist_c,sub(simpBoundary(intC_c,intC_(c+1)),S))
    )))
)

------------------------------------------


------------------------------------------
idealsComplex=method(Options=>{
    symbol Homogenize => true, 
    symbol VariableName => getSymbol "t",
    symbol CoefficientRing => QQ,
    symbol BaseRing => null})
------------------------------------------
--This function computes the mixed degree generalization
-- of the Schenck-Stillman chain complex
-- of ideals for a simplicial or polyhedral complex
------------------------------------------

------------------------------------------
--Inputs:
--V: vertex list
--F: facet list
--r: desired order of smoothness
--degD: degree deficits on facets
------------------------------------------
--Outputs: The mixed degree Schenck-Stillman complex of ideals
------------------------------------------
idealsComplex(List,List,ZZ,ZZ):=ChainComplex => opts -> (V,F,r,degD)->(
    degD = apply(#F, i -> degD);
    intC := getInteriorFaces(F);
    E := intC_1;
    r = apply(#E, i -> r);
    idealsComplex(V,F,{E,r},degD,opts)
)

idealsComplex(List,List,ZZ,List):=ChainComplex => opts -> (V,F,r,degD)->(
    intC := getInteriorFaces(F);
    E := intC_1;
    r = apply(#E, i -> r);
    idealsComplex(V,F,{E,r},degD,opts)
)

------------------------------------------
--Inputs:
--V: vertex list
--F: facet list
--Er: {list of comdim-1 facets, list/value of smoothness} 
--degD: degree deficits on facets
------------------------------------------
--Outputs: The mixed degree Schenck-Stillman complex of ideals
------------------------------------------
idealsComplex(List,List,List,ZZ):=ChainComplex => opts -> (V,F,Er,degD)->(
    degD = apply(#F, i -> degD);
    idealsComplex(V,F,Er,degD,opts)
)

idealsComplex(List,List,List,List):=ChainComplex => opts -> (V,F,Er,degD)->(
    E := Er_0;
    r := Er_1;
    try #E==#r then () else ( r = apply(#E, i -> r); );
    d := #(first V);
    S := createSplineRing(d,opts);
    --list of interior faces in order of increasing codimension--
    intC := getInteriorFaces(F);
    r = arrangeSmoothnessIndices(F, E, r);
    --get degree deficits on all faces in increasing codimension
    degDAll := getAllDegreeDeficits(F, degD);
    --list of forms defining codim 1 interior faces
    intformslist :=formsList(V,intC_1,r,degDAll_1,opts);
    --list of modules which will define chain complex--
    fullmodulelist := apply(#intC,i->(
        if #(intC_i)==0 then(
            newMod := image matrix{{0_S}};
        )else(
            newMod = directSum apply(intC_i,e->(
                CE := positions(intC_1,f->subsetL(e,f));
                sub(module ideal (intformslist_CE),S)
            ));
        );
        newMod
    ));
    --defining the chain complex
    CCSS :=chainComplex(reverse apply(#intC-1, c-> (
        inducedMap(fullmodulelist_(c+1),fullmodulelist_c,sub(simpBoundary(intC_c,intC_(c+1)),S))
    )))
)

------------------------------------------


------------------------------------------
splineComplex=method(Options=>{
    symbol Homogenize => true, 
    symbol VariableName => getSymbol "t",
    symbol CoefficientRing => QQ,
    symbol BaseRing => null})
------------------------------------------
--This function computes the mixed degree 
-- generalization of the Schenck-Stillman chain complex
--of quotient rings for a simplicial complex,
--called the spline complex
------------------------------------------

------------------------------------------
--Inputs:
--V: vertex list
--F: facet list
--r: desired order of smoothness
--degD: degree deficits on facets
------------------------------------------
--Outputs: The mixed degree spline complex
------------------------------------------
splineComplex(List,List,ZZ,ZZ):=ChainComplex => opts -> (V,F,r,degD)->(
    degD = apply(#F, i -> degD);
    intC := getInteriorFaces(F);
    E := intC_1;
    r = apply(#E, i -> r);
    splineComplex(V,F,{E,r},degD,opts)
)

splineComplex(List,List,ZZ,List):=ChainComplex => opts -> (V,F,r,degD)->(
    intC := getInteriorFaces(F);
    E := intC_1;
    r = apply(#E, i -> r);
    splineComplex(V,F,{E,r},degD,opts)
)

------------------------------------------
--Inputs:
--V: vertex list
--F: facet list
--Er: {list of comdim-1 facets, list/value of smoothness} 
--degD: degree deficits on facets
------------------------------------------
--Outputs: The mixed degree Schenck-Stillman complex of ideals
------------------------------------------
splineComplex(List,List,List,ZZ):=ChainComplex => opts -> (V,F,Er,degD)->(
    degD = apply(#F, i -> degD);
    splineComplex(V,F,Er,degD,opts)
)

splineComplex(List,List,List,List):=ChainComplex => opts -> (V,F,Er,degD)->(
    E := Er_0;
    r := Er_1;
    try #E==#r then () else ( r = apply(#E, i -> r); );
    d := #(first V);
    S := createSplineRing(d,opts);
    --list of interior faces in order of increasing codimension--
    intC := getInteriorFaces(F);
    r = arrangeSmoothnessIndices(F, E, r);
    --get degree deficits on all faces in increasing codimension
    degDAll := getAllDegreeDeficits(F, degD);
    --list of forms defining codim 1 interior faces
    intformslist :=formsList(V,intC_1,r,degDAll_1,opts);
    --list of modules which will define chain complex--
    varlist := vars S;
    hvar := varlist_(0,d);
    fullmodulelist := { directSum apply(#F, i-> module ideal ( hvar^((degDAll_0)_i) ) ) };
    scan(#intC-1,i->(
        if #(intC_(i+1))==0 then(
            newMod := image matrix{{0_S}};
        )else(
            newMod = directSum apply(#(intC_(i+1)),j->(
                e := (intC_(i+1))_j;
                CE := codim1Cont(intC_1,e);
                modn := sub(module ideal (hvar^((degDAll_(i+1))_j)),S);
                modd := sub(module ideal (intformslist_CE),S);
                return (modn/modd)
            ));
        );
        fullmodulelist = append(fullmodulelist,newMod)
    ));
    --defining the chain complex
    CCSS :=chainComplex(reverse apply(#intC-1, c-> (
        inducedMap(fullmodulelist_(c+1),fullmodulelist_c,sub(simpBoundary(intC_c,intC_(c+1)),S))
    )))
)

------------------------------------------


------------------------------------------
codim1Cont=method()
------------------------------------------

------------------------------------------
----Inputs:
----E = list of codim one intersections
----G = a face of the complex (V,F)
------------------------------------------
----Output: Indices of E corresponding to
-- codim one faces containing G
-----------------------------------------
codim1Cont(List,List):=List=> (E,G)->(
    --positions of codim one faces containing G
    positions(E,e->subsetL(G,e))
)

------------------------------------------


-----------------------------------------
-----------------------------------------


------------------------------------------
-- Documentation
------------------------------------------

beginDocumentation()

-- Front Page
doc ///
    Key
        SimplicialMDSplines
    Headline
        splines of mixed degrees and mixed orders of smoothness on simplicial complexes
    Description
        Text
            This package provides tools for computing the dimension of smooth piecewise polynomial splines over simplicial complexes.
            In particular, we consider splines spaces of mixed (or non-uniform) degrees and orders of smoothness.
--        Text
--            Given a simplicial parition $\mathcal T$, let maps $\mathbf \Delta m$ and $\mathbf{r}$ assign non-negative integers to its faces and interior edges, respectively.
--            Then, given a non-negative integer $n$, we consider a spline space $S^{\mathbf{r}}_{\mathbf{\Delta m}, n}$ such that, for $f \in S^{\mathbf{r}}_{\mathbf{\Delta m}, n}$, $f$ is $C^{\mathbf{r}(\tau)}$ smooth across the interior edge $\tau$, and $f|_{\sigma}$ is a polynomial of degree $n - \mathbf{\Delta m}(\sigma)$ for a face $\sigma$.
        Text
            The dimension of uniform degree spline spaces can be computed (or estimated) using the tools of homological algebra developed in the context of splines in the seminal paper of Billera [1].
            Building a short exact sequence of chain complexes, [1] showed that the spline space dimension is equal to the dimension of a graded piece of a particular homology module.
            A modification that improved upon [1] and yielded simpler chain complexes was presented by Schenck and Stillman [2] and used to prove that the dimension formula derived by Schumaker [3] holds in sufficiently high degree.
            In particular, this demonstrated that the homological algebra approach of the former agrees with the Bernstein–Bezier approach of the latter.
            These modified chain complexes were further studied, for instance, in [4-6].
        Text
            This package is meant to accompany the report by Toshniwal and Hughes [7] which generalizes the Schenck-Stillman chain complexes to the case of mixed degree splines or non-uniform degree splines.
            This package has been derived from `AlgebraicSplines' written by Mike DiPasquale.
        Text
            [1] L.J. Billera, Homology of smooth splines: Generic triangulations and a conjecture of Strang, Trans. Amer. Math. Soc. 310 (1988).\break
            [2] Hal Schenck and Mike Stillman. Local cohomology of bivariate splines, J. Pure Appl. Algebra 117/118 (1997).\break
            [3] L.L. Schumaker, On the Dimension of Spaces Of Piecewise Polynomials in Two Variables, in: Multivariate Approximation Theory, Birkhäuser, Basel, 1979.\break
            [4] A. Geramita and H. Schenck, Fat points, inverse systems, and piecewise polynomial functions, J. of Algebra 204(1) (1998).\break
            [5] Hal Schenck, A spectral sequence for splines, Adv. in Appl. Math. 19 (1997).\break
            [6] B. Mourrain, N. Villamizar, Homological techniques for the analysis of the dimension of triangular spline spaces, J. of Symb. Comp. 50 (2013).\break
            [7] D. Toshniwal, T.J.R. Hughes, Polynomial splines of non-uniform degree on triangulations: Combinatorial bounds on the dimension, in: Computer Aided Geometric Design (2019; to appear).
        Text
            We reproduce two examples from [7] in the following.
            The first example (Example 2, [7]) features a triangulation with two isolated interior vertices.
        Example
            --Degree range
            a = 0;b = 20;
            --Triangulation
            V = {{0/1,0/1}, {1/1,0/1}, {1/1,1/1}, {0/1,1/1}, {1/4,1/3}, {3/4,2/3}};
            F = {{0,1,4}, {0,3,4}, {1,2,5}, {1,3,5}, {2,3,5}, {1,3,4}};
            E = {{0,4}, {1,4}, {3,4}, {1,5}, {2,5}, {3,5}, {1,3}};
            --Spline space configuration
            degD = {3, 3, 0, 1, 0, 2};
            r = {0, 0, 0, 2, 3, 2, 1};
            --Relevant dimensions in different degrees
            homologyDimensionTable(a,b,{V,F,{E,r},degD})
        Text
            The second example (Example 3, [7]) features a triangulation with 3 interior vertices that share 3 edges between them.
        Example
            --Degree range
            a = 0; b = 20;
            --Triangulation
            V = {{-230/389,-1637/7625}, {-557/1706,1460/1447}, {-2109/1532,1114/935}, {1406/383,-1309/412}, {4288/2033,655/139}, {-2497/633,-3158/1011}, {-3515/766,2228/561}};
            F = {{0,1,2}, {2,4,6}, {0,1,3}, {1,3,4}, {1,2,4}, {0,2,5}, {2,5,6}, {0,3,5}};
            E = {{0,1}, {1,2}, {0,2}, {2,4}, {2,6}, {1,3}, {0,3}, {1,4}, {2,5}, {0,5}};
            --Spline space configuration
            degD = {0, 1, 1, 1, 1, 1, 1, 1};
            r = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
            --Relevant dimensions in different degrees
            homologyDimensionTable(a,b,{V,F,{E,r},degD})
 
///
