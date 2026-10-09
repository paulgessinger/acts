import csv, sys, math, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
D=sys.argv[1]
g1=list(csv.DictReader(open(D+"/gen1_approach.csv")))
modes=[("Gen3 default (layers expand in z)","gen3"),("Gen3 --zgaps (Gen1-like z structure)","gen3_zgaps")]
fig,axs=plt.subplots(2,3,figsize=(17,9),sharey="row")
for row,(rmax,zmax) in enumerate([(125,300),(900,1150)]):
    ax=axs[row][0]; ax.set_title("Gen1 TGeo: approach surfaces")
    for r in g1:
        c="tab:red" if r["vol"]=="12" and r["lay"] in ("4","8") else "tab:blue"
        ax.plot([float(r["z0"]),float(r["z1"])],[float(r["r"])]*2,c,lw=1.2)
    for i,(t,f) in enumerate(modes):
        ax=axs[row][i+1]; ax.set_title(t+": material portals")
        seen=set()
        for r in csv.DictReader(open(f"{D}/{f}_material_surfaces.csv")):
            if r["geoid"] in seen or r["type"]!="cyl": continue
            seen.add(r["geoid"])
            ax.plot([float(r["z0"]),float(r["z1"])],[float(r["r0"])]*2,"tab:green",lw=1.2)
        # volume outlines
        for v in csv.DictReader(open(f"{D}/{f}_volumes.csv")):
            if v["name"]=="World": continue
            rmin,rmx,hz,cz=(float(v[k]) for k in ("rmin","rmax","hz","cz"))
            ax.add_patch(plt.Rectangle((cz-hz,rmin),2*hz,rmx-rmin,fill=False,ec="0.75",lw=0.4))
    for ax in axs[row]:
        ax.set_xlim(-zmax,zmax); ax.set_ylim(0,rmax); ax.set_xlabel("z [mm]")
        for eta in (0.7,):
            for s in (1,-1):
                zz=rmax/math.tan(2*math.atan(math.exp(-eta)))
                ax.plot([0,s*zz],[0,rmax],":",c="0.5",lw=0.8)
    axs[row][0].set_ylabel("r [mm]")
axs[0][0].text(-290,115,"red: Silicon (vol 12) layers 4/8\ndotted: |eta| = 0.7",fontsize=8,va="top")
fig.tight_layout(); fig.savefig(D+"/gen1_vs_gen3_rz.png",dpi=110)
