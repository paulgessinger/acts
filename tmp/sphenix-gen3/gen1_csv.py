import math, acts, sys
from acts.examples.tgeo import TGeoDetector
# usage: gen1_csv.py <geometry.root> <tgeo-config.json> <out.csv>
cfg = TGeoDetector.Config(); cfg.fileName=sys.argv[1]; cfg.jsonFile=sys.argv[2]
cfg.surfaceLogLevel=cfg.layerLogLevel=cfg.volumeLogLevel=acts.logging.FATAL
det=TGeoDetector(cfg); tg=det.trackingGeometry(); gctx=acts.GeometryContext.dangerouslyDefaultConstruct()
out=open(sys.argv[3],"w"); out.write("vol,lay,r,z0,z1\n")
def vs(s):
    g=s.geometryId
    if g.approach:
        b=list(s.bounds.values()); cz=s.center(gctx)[2]
        out.write(f"{g.volume},{g.layer},{b[0]:.2f},{cz-b[1]:.2f},{cz+b[1]:.2f}\n")
tg.visitSurfaces(vs, False)
