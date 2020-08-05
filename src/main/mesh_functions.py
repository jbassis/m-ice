import numpy as np
from   fenics import *

def GlobalizeMesh(mesh):
    mesh.init()
    #plot(mesh, interactive=True)
    TopDim = mesh.topology().dim()
    GeoDim = mesh.geometry().dim()

    NumVertexGlobal = mesh.num_entities_global(0)
    NumCellsGlobal = mesh.num_entities_global(TopDim)

    coords = mesh.coordinates()
    CoordsGlobal = np.zeros([NumVertexGlobal, GeoDim])
    NumProcGlobal = np.zeros(NumVertexGlobal)
    for v in vertices(mesh):
        x = np.array(coords[v.index()])
        CoordsGlobal[v.global_index(),:] = x
        NumProcGlobal[v.global_index()] = 1.0
    #SyncSum(CoordsGlobal)
    #SyncSum(NumProcGlobal)
    if min(NumProcGlobal==0):
        self.dprint("ERROR: THERE IS A VERTEX NOT BELONGING TO ANY PROCESSOR")
        exit(1)
    for IndexGlobal in range(NumVertexGlobal):
        CoordsGlobal[IndexGlobal,:] = CoordsGlobal[IndexGlobal,:]/NumProcGlobal[IndexGlobal]
    CellsGlobal = np.zeros([NumCellsGlobal, TopDim+1], dtype=np.uintp)
    for c in cells(mesh):
        LocVertexID = 0
        for v in vertices(c):
            CellsGlobal[c.global_index(), LocVertexID] = v.global_index()
            LocVertexID += 1
    #SyncSum(CellsGlobal)
    #New: Also globalize Facets:
    mesh.init_global(TopDim-1)
    NumFacetsGlobal = mesh.num_entities_global(TopDim-1)
    NumProcFacets = np.zeros(NumFacetsGlobal)
    FacetsGlobal = np.zeros([NumFacetsGlobal, TopDim], dtype=np.uintp)
    for f in facets(mesh):
        LocVertexID = 0
        #For some VERY strange reason, faces in SymmetricMesh segfault upon global_index() call
        #And ONLY on 1 processor...

        f_index = f.index()
        for v in vertices(f):
            FacetsGlobal[f_index, LocVertexID] = v.global_index()
            LocVertexID += 1
        NumProcFacets[f_index] = 1.0
    #SyncSum(FacetsGlobal)
    #SyncSum(NumProcFacets)
    for IndexGlobal in range(NumFacetsGlobal):
        FacetsGlobal[IndexGlobal,:] = FacetsGlobal[IndexGlobal,:]/NumProcFacets[IndexGlobal]
    return (CoordsGlobal, CellsGlobal, FacetsGlobal)


def sort_boundary_nodes(mesh):
    """
    Boundary nodes are not provided by BoundaryMesh in order and so we loop through
    facets and provide a closed polygon of pts defining the shape of the
    domain that should have been provided by BoundaryMesh
    """
    (x, CellsGlobal, FacetsGlobal) = GlobalizeMesh(mesh)
    xc = []
    yc = []
    Cells=np.copy(CellsGlobal)
    idx = 0
    start =Cells[0,0]
    next  =Cells[0,1]

    ncells = len(Cells)
    for i in range(len(Cells)):
        xc.append(x[start,0])
        yc.append(x[start,1])
        Cells[idx,0]=ncells+1
        Cells[idx,1]=ncells+1
        # Now look in first column of CellsGlobal for next
        idx = np.where(Cells[:,0]==next)[0]
        #start=[0]
        if len(idx)==1:
            idx = idx[0]
            start = Cells[idx,0]
            next =  Cells[idx,1]
            #start = np.where(Cells[:,1]==next)
            #next  = Cells[]
        else:
            idx = np.where(Cells[:,1]==next)[0]
            if len(idx)==1:
                idx = idx[0]
                start = Cells[idx,1]
                next =  Cells[idx,0]
            #else:
            #    break
    return np.array([xc,yc]).transpose()
