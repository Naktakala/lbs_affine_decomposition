-- 2D Transport test with Vacuum and Incident-isotropic BC.
-- SDM: PWLD
-- Test: Max-value=0.50758 and 2.52527e-04
num_procs = 1





--############################################### Check num_procs
if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
                      "Expected "..tostring(num_procs)..
                      ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--############################################### Setup mesh
chiMeshHandlerCreate()

mesh={}
N=32
L=2.0
xmin = -1.0
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end
chiMeshCreateUnpartitioned3DOrthoMesh(mesh,mesh,{0.0,0.1,0.2})
chiVolumeMesherExecute();

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)


--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)


num_groups = 168
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        PDT_XSFILE,"Tests/xs_3_170.data")
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
        PDT_XSFILE,"Tests/xs_3_170.data")

--chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,SIMPLEXS0,num_groups,0.1)
--chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,SIMPLEXS0,num_groups,0.1)

src={}
for g=1,num_groups do
    src[g] = 0.0
end
--src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Setup Physics
phys1 = chiAffineDecompositionerCreate()
chiSolverAddRegion(phys1,region1)

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 1)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,3)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetAngleAggregationType(phys1,cur_gs,LBSGroupset.ANGLE_AGG_SINGLE)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,2)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)

-- gs1 = chiLBSCreateGroupset(phys1)
-- cur_gs = gs1
-- chiLBSGroupsetAddGroups(phys1,cur_gs,63,167)
-- chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
-- chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
-- chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,2)
-- chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES)
-- chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
-- chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
-- chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)
-- --chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,false," ")
-- --chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false," ")

--############################################### Set boundary conditions
bsrc={}
for g=1,num_groups do
    bsrc[g] = 0.0
end
bsrc[1] = 1.0/4.0/math.pi
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,XMIN,
                        LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);

chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
chiLBSSetProperty(phys1,SCATTERING_ORDER,1)
chiLBSSetProperty(phys1,LBSProperty.SAVE_ANGULAR_FLUX,true)

--############################################### Initialize and Execute Solver
chiAffineDecompositionerInitialize(phys1)
chiAffineDecompositionSetProperty(phys1,"num_modes",1)

chiAffineDecompositionReadPODModes(phys1,"ZPsi",0);


