import os, osgeo, math, numpy as np #shapely, shapefile
from osgeo import gdal, ogr, osr
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PatchCollection
%matplotlib inline
%pylab inline

#test

#brian's test.

root = 'I:/projects/FATA_LULC/test/'
input_shp = 'FATA_villages_CIA_Waziristans.shp'
fieldUID = 'Id_buggy'
outDir = 'I:/projects/FATA_LULC/test/'

proj = 4326 #http://spatialreference.org/ref/epsg/wgs-84/

''' ^^^^^^^^^^^^^^^^^^^^^^^^ '''
''' load point data from shp '''
''' vvvvvvvvvvvvvvvvvvvvvvvv '''

drv = ogr.GetDriverByName('ESRI Shapefile')
input_shp_full = root+input_shp
ptShp = drv.Open(input_shp_full)
ptLayer = ptShp.GetLayer(0)
ptSRS = ptLayer.GetSpatialRef()
x_min, x_max, y_min, y_max = ptLayer.GetExtent()

ptList = []
ptDict = {}
for pt in ptLayer:
    ID_index = pt.GetFieldIndex(fieldUID)
    ptID = pt.GetField(ID_index)
    ptGeom = pt.GetGeometryRef()
    ptX = float(str(ptGeom).split(' ')[1].strip('('))
    ptY = float(str(ptGeom).split(' ')[2].strip(')'))
    ptDict[ptID] = [ptX,ptY]
    ptList.append([ptX,ptY]) 
    
numPtList = np.array(ptList) # used for input to Delaunay   

''' ^^^^^^^^^^^^^^^^^^^^^^^^ '''
''' construct radial buffers '''
''' vvvvvvvvvvvvvvvvvvvvvvvv '''

drv = ogr.GetDriverByName('Esri Shapefile')
buffShp = outDir+"test_bufferFATA.shp"
ds = drv.CreateDataSource(buffShp)
layer = ds.CreateLayer('', None, ogr.wkbPolygon)
layer.CreateField(ogr.FieldDefn('Id', ogr.OFTInteger))
defn = layer.GetLayerDefn()

ptCounter = 0
for each in ptDict:
    ptLon = ptDict[each][0]
    ptLat = ptDict[each][1]
    pt_wkt = "POINT ("+str(ptLon)+' '+str(ptLat)+')'
    pt = ogr.CreateGeometryFromWkt(pt_wkt)
    
    bufferDistance = 0.01 # degrees
    geom = pt.Buffer(bufferDistance)

    feat = ogr.Feature(defn)
    feat.SetField('Id', each)
    feat.SetGeometry(geom)  
    layer.CreateFeature(feat)
    feat = geom = None
    ptCounter+=1

print str(ptCounter)+' points buffered'
layer = ds = None

''' ^^^^^^^^^^^^^^^^^^^^^^^^ '''
''' calc delaunay parameters '''
''' vvvvvvvvvvvvvvvvvvvvvvvv '''

# http://en.wikipedia.org/wiki/Circumscribed_circle#Circumscribed_circles_of_triangles
# https://stackoverflow.com/questions/12374781/how-to-find-all-neighbors-of-a-given-point-in-a-delaunay-triangulation-using-sci
# https://stackoverflow.com/questions/10650645/python-calculate-voronoi-tesselation-from-scipys-delaunay-triangulation-in-3d

def dot2(u, v):
    return u[0]*v[0] + u[1]*v[1]

def cross2(u, v, w):
    """u x (v x w)"""
    return dot2(u, w)*v - dot2(u, v)*w

def ncross2(u, v):
    """|| u x v ||^2"""
    return sq2(u)*sq2(v) - dot2(u, v)**2

def sq2(u):
    return dot2(u, u)

tri = Delaunay(numPtList)
p = tri.points[tri.vertices]# returns point locations of each triangle vertex
                            # tri.vertices: each row represents one simplex (triangle) in the triangulation,
                            # with values referencing indices of the input point list
print 'finding facets...'
# find facets containing each input point
triDict={}
i = 0
while i < len(numPtList):
    pt = str(numPtList[i].tolist())
    j = 0
    while j < len(p):
        k = 0
        while k < len(p[j]):
            if pt == str(p[j][k].tolist()):
                if triDict.has_key(pt):
                    triDict[pt]+=[p[j].tolist()] # format: vorDict['pt coordinates']=list of lists of triangle vertices
                else:
                    triDict[pt]=[p[j].tolist()]
            k+=1
        j+=1
    i+=1

print 'finding circumcenters...'
# find circumcenters for triangles associated with each pt; these circumcenters are Voronoi vertices
vorDict={}
for each in triDict:
    npPt = np.array(triDict[each])
    A = npPt[:,0,:].T # 1st vertex of triangle
    B = npPt[:,1,:].T # 2nd vertex
    C = npPt[:,2,:].T # 3rd vertex
    a = A - C
    b = B - C
    cc = cross2(sq2(a) * b - sq2(b) * a, a, b) / (2*ncross2(a, b)) + C # coords of circumcenters; coords of Voronoi edges
    vorDict[each]=cc.T

print 'finding Voronoi nodes...'
#take ID from ptDict and link it to Voronoi vertices from vorDict
vorIdDict={}
for a in ptDict:
    for b in vorDict:
        if str(ptDict[a])==b:
            #print ptDict[a],b
            vorIdDict[a]=vorDict[b]
            
print 'converting nodes to shp...'
# convert to shp
vorShp = outDir+"test_NodesFATA.shp"
drv = ogr.GetDriverByName('ESRI Shapefile')
if os.path.exists(vorShp): drv.DeleteDataSource(vorShp)
outVorShp = drv.CreateDataSource(vorShp)
vorLayer = outVorShp.CreateLayer('', None,ogr.wkbPoint)
vorLayer.CreateField(ogr.FieldDefn('Id', ogr.OFTInteger))
vorDefn = vorLayer.GetLayerDefn()

ptCounter = 0
for each in vorIdDict:
    i = 0
    while i < len(vorIdDict[each]):
        ptLon = vorIdDict[each][i].tolist()[0]
        ptLat = vorIdDict[each][i].tolist()[1]
        pt_wkt = "POINT ("+str(ptLon)+' '+str(ptLat)+')'
        geom = ogr.CreateGeometryFromWkt(pt_wkt)
        feat = ogr.Feature(vorDefn)
        feat.SetField('Id', each)
        feat.SetGeometry(geom)  
        vorLayer.CreateFeature(feat)
        feat = geom = None
        ptCounter+=1
        i+=1

print 'exported '+str(ptCounter)+' Voronoi nodes'
vorLayer = outVorShp = None

''' ^^^^^^^^^^^^^^^^^^^^^^^^ '''
''' make voronoi polys '''
''' vvvvvvvvvvvvvvvvvvvvvvvv '''

# https://stackoverflow.com/questions/1709283/how-can-i-sort-a-coordinate-list-for-a-rectangle-counterclockwise
# https://gamedev.stackexchange.com/questions/13229/sorting-array-of-points-in-clockwise-order
# https://en.wikipedia.org/wiki/Graham_scan
def sortCCW(node):
    return math.atan2(node[1] - meanLat, node[0] - meanLon)

shp = outDir+"test_VoronoiFATA.shp"
drv = ogr.GetDriverByName('ESRI Shapefile')
outShp = drv.CreateDataSource(shp)
layer = outShp.CreateLayer('', None,ogr.wkbPolygon)
layer.CreateField(ogr.FieldDefn('Id', ogr.OFTInteger))
layerDefn = layer.GetLayerDefn()

for pt in vorIdDict:
    print pt
    meanLon = sum(node[0] for node in vorIdDict[pt])/len(vorIdDict[pt])
    meanLat = sum(node[1] for node in vorIdDict[pt])/len(vorIdDict[pt]) 
    hullList = vorIdDict[pt].tolist()
    hullList.sort(key=sortCCW)
  
    poly = ogr.Geometry(ogr.wkbPolygon)
    ring = ogr.Geometry(ogr.wkbLinearRing)
    i = 0
    for node in hullList:
        if i==0:
            loopLon = node[0] # grab first node to close ring
            loopLat = node[1]
        ring.AddPoint(node[0],node[1])
        i+=1
    ring.AddPoint(loopLon,loopLat)
    poly.AddGeometry(ring)
    feat = ogr.Feature(layerDefn)
    feat.SetField('Id', pt)
    feat.SetGeometry(poly)  
    layer.CreateFeature(feat)
    feat = poly = ring = None
layer = outShp = None
