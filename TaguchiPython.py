# Keaton Harlow
#
# This project was done to compare different methods of designing an experiment
# and choosing intelligent starting points for optimization.
#
# This program explores the use of the Taguchi method to choose starting points.
#
# The part is created in ABAQUS, a distributed force is applied, and the
# max Mises stress is collected from ABAQUS. The data is then output
# to a .csv file for comparison to other methods.

#Import tools from ABAQUS
from abaqus import * 
from abaqusConstants import * 
import __main__ 
import section 
import odbSection 
import regionToolset 
import displayGroupMdbToolset as dgm 
import part 
import material 
import assembly 
import step 
import interaction 
import load 
import mesh 
import job 
import sketch 
import visualization 
import xyPlot 
import connectorBehavior 
import displayGroupOdbToolset as dgo 
from odbAccess import *
from math import atan, sin, cos, tan 
from MaxMises import getMaxMises

taguchiMatrix = [[1,1,1,1],[1,2,2,2],[1,3,3,3],[1,4,4,4],[2,1,2,3],[2,2,1,4],[2,3,4,1],[2,4,3,2 ,], [3,1,3,4],[3,2,4,3],[3,3,1,2],[3,4,2,1],[4,1,4,2],[4,2,3,1],[4,3,2,4],[4,4,1,3]]
fileName = 'taguchi.csv' 
DataFile = open(fileName,'w') 
DataFile.write('Run Num, h_L, f_d1, f_d2, t, maxStress, mass, maxStress / m\n') 
DataFile.close()

#Initializing the design constraints 
L = 1.0 
h0 = 0.5 
F = 15000.0 
xFirstHole = .25 
xSecondHole = .75 
rho = 2720

#Design parameter lists 
hL = [.1 * h0, .4 * h0, .7 * h0, h0] 
fd1 = [.1, .35, .6, .85] 
fd2 = [.1, .35, .6, .85] 
t = [.02 * h0, .05 * h0, .075 * h0, .1 * h0]
ModelName = 'Model-1' 
pi = 3.14159

#Looping through the constraint
for i in range(0, 16): 
    #Retrieving vars from lists
    h = hL[taguchiMatrix[i][0] - 1] 
    f1 = fd1[taguchiMatrix[i][1] - 1] 
    f2 = fd2[taguchiMatrix[i][2] - 1] 
    th = t[taguchiMatrix[i][3] - 1]
    
    #Calculated values temp 
    hFirstHole = ((h - h0) / L * xFirstHole + h0) 
    hSecondHole = ((h - h0) / L * xSecondHole + h0) 
    d1 = f1 * hFirstHole 
    d2 = f2 * hSecondHole 
    tao = F / (th * h)
    Mdb() 
    
    #Creating the design space 
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2.0) 
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints 
    s.setPrimaryObject(option=STANDALONE)

    #Drawing the trapezoid 
    s.Line(point1=(0.0, 0.0), point2=(0.0, h0)) 
    s.VerticalConstraint(entity=g[2], addUndoState=False) 
    s.Line(point1=(0.0, h0), point2=(L, h))
    s.Line(point1=(L, h), point2=(L, 0.0)) 
    s.VerticalConstraint(entity=g[4], addUndoState=False) 
    s.Line(point1=(L, 0.0), point2=(0.0, 0.0)) 
    s.HorizontalConstraint(entity=g[5], addUndoState=False) 
    s.PerpendicularConstraint(entity1=g[4], entity2=g[5], addUndoState=False)

    #Drawing the first hole 
    s.CircleByCenterPerimeter(center=(xFirstHole, hFirstHole / 2.0), point1=(xFirstHole - d1 / 2.0, hFirstHole / 2.0))

    #Drawing the second hole 
    s.CircleByCenterPerimeter(center=(xSecondHole, hSecondHole / 2.0), point1=(xSecondHole - d2 / 2.0, hSecondHole / 2.0))

    #Creating the shape 
    p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY) p = mdb.models['Model-1'].parts['Part-1'] p.BaseShell(sketch=s) s.unsetPrimaryObject() p = mdb.models['Model-1'].parts['Part-1'] del mdb.models['Model-1'].sketches['__profile__']

    #Creating the material 
    mdb.models['Model-1'].Material(name='Aluminum') 
    mdb.models['Model-1'].materials['Aluminum'].Elastic(table=((70000000000.0, 0.3), ))

    #Creating the section
    mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1', material='Aluminum', thickness=th)

    #Assigning the section 
    p = mdb.models['Model-1'].parts['Part-1'] 
    f = p.faces 
    faces = f.getSequenceFromMask(mask=('[#1 ]', ), ) 
    region = p.Set(faces=faces, name='Set-1') 
    p = mdb.models['Model-1'].parts['Part-1'] 
    p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    #Instancing the part 
    a = mdb.models['Model-1'].rootAssembly 
    a.DatumCsysByDefault(CARTESIAN) 
    p = mdb.models['Model-1'].parts['Part-1'] 
    a.Instance(name='Part-1-1', part=p, dependent=ON)

    #Creating the step 
    mdb.models['Model-1'].StaticStep(name='Shear', previous='Initial')

    #Creating the load 
    v11 = a.instances['Part-1-1'].vertices 
    v21 = a.instances['Part-1-1'].vertices 
    a = mdb.models['Model-1'].rootAssembly 
    s1 = a.instances['Part-1-1'].edges 
    side1Edges1 = s1.getSequenceFromMask(mask=('[#8 ]', ), ) 
    region = a.Surface(side1Edges=side1Edges1, name='Surf-1') 
    mdb.models['Model-1'].SurfaceTraction(name='Load-1', createStepName='Shear', region=region, magnitude=tao, directionVector=(v11[4], v21[3]), distributionType=UNIFORM, field='', localCsys=None, resultant=ON)

    #Creating the boundary condition 
    a = mdb.models['Model-1'].rootAssembly 
    e1 = a.instances['Part-1-1'].edges 
    edges1 = e1.getSequenceFromMask(mask=('[#20 ]', ), ) 
    region = a.Set(edges=edges1, name='Set-1') 
    mdb.models['Model-1'].EncastreBC(name='BC-1', createStepName='Shear', region=region, localCsys=None)

    #Choosing the element type and seeding the part 
    p = mdb.models['Model-1'].parts['Part-1'] 
    f = p.faces 
    pickedRegions = f.getSequenceFromMask(mask=('[#1 ]', ), ) 
    p.setMeshControls(regions=pickedRegions, elemShape=TRI) 
    p.seedPart(size=0.02, deviationFactor=0.1, minSizeFactor=0.1) 
    elemType1 = mesh.ElemType(elemCode=CPS8R, elemLibrary=STANDARD) 
    elemType2 = mesh.ElemType(elemCode=CPS6, elemLibrary=STANDARD) 
    p = mdb.models['Model-1'].parts['Part-1'] 
    f = p.faces 
    faces = f.getSequenceFromMask(mask=('[#1 ]', ), ) pickedRegions =(faces, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))

    #Generating the mesh 
    p = mdb.models['Model-1'].parts['Part-1'] 
    p.generateMesh() 
    a = mdb.models['Model-1'].rootAssembly a.regenerate()

    #Creating the job 
    mdb.Job(name=ModelName, model=ModelName, description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', multiprocessingMode=DEFAULT, numCpus=8, numDomains=8, numGPUs=0)
    job = mdb.jobs[ModelName]

    #Delete lock file 
    if os.access('%s.lck'%ModelName,os.F_OK): 
        os.remove('%s.lck'%ModelName)

        #Running the job 
        ob.submit() 
        job.waitForCompletion()

        #Get max mises stress 
        odbFile = ModelName + '.odb' 
        maxMisesStress = getMaxMises(odbFile, '') 
        odb = openOdb(odbFile) odb.close()

        #Computing the mass 
        V = ((h0 + h) / 2.0 * L - .5 * pi * pow((d1 / 2.0), 2)- .5 * pi * pow((d2 / 2.0), 2)) * th 
        mass = V * rho

        #Writing data to a .csv file 
        DataFile = open(fileName,'a') 
        DataFile.write('%f, %f, %f, %f, %f, %f, %f, %f\n'%(i + 1, h, f1, f2, th, maxMisesStress, mass, maxMisesStress / mass)) 
        DataFile.close()