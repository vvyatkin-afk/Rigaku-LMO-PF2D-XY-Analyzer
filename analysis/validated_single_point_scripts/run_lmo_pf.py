"""Read-only reconstruction of the supplied 360-frame LMO(011) RASX acquisition.
Creates a new results directory; never modifies source files, Studio, or instrument settings.
Requires numpy, scipy, matplotlib (available in the installed Python 3.14).
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,pathlib,subprocess,sys,zipfile,datetime,platform,shutil
import xml.etree.ElementTree as E
import numpy as np
SOURCE="."
PATTERN="011_2D_LMO_011_PF2D_TEST_image_*.rasx"
EXPECTED={
 "PXD_DETECTOR_DIMENSIONS":[775,385],"PXD_DETECTOR_SIZE":[77.5,38.5],
 "PXD_DETECTOR_VECTORS":[1,0,0,0,1,0],"PXD_SPATIAL_DISTORTION_VECTORS":[0,-1,1,0],
 "PXD_GONIO_VECTORS":[0,0,1,-1,0,0,0,1,0,1,0,0,0,1,0,0,0,-1],
 "PXD_GONIO_VALUES":[270,32.6,0,0,0,150.608],
 "PXD_SPATIAL_DISTORTION_INFO":[382.5,197.5,.1,.1],
 "PXD_SPATIAL_BEAM_POSITION":[382.5,197.5],
 "CRYSTAL_GONIO_VECTORS":[-1,0,0,0,0,-1,1,0,0],
 "SOURCE_VECTORS":[0,0,1,0,1,0,1,0,0],
 "PXD_DETECTOR_ACQUISITION_EXPTIME":[5000],
 "PXD_DETECTOR_REGION_OF_INTEREST":[1,1,775,385],
 "PXD_DETECTOR_IMAGE_CORRECTION_FLAGS":[2,1,1,1,0,0],
 "SOURCE_WAVELENGTH":[1.,1.54059292]
}
def numbers(text):return np.array([float(x) for x in text.split()])
def check_values(name,actual,expected,filename):
 a=numbers(actual);b=np.asarray(expected)
 if a.shape!=b.shape or not np.allclose(a,b,rtol=0,atol=1e-5):
  raise ValueError(f"Unsupported or inconsistent {name} in {filename}: {actual}; expected {expected}")
def extract(source,out,reference):
 files=sorted(source.glob(PATTERN))
 if len(files)!=360:raise ValueError(f"Expected exactly 360 supplied frames matching {PATTERN}, found {len(files)}")
 frames=[];manifest=[];uniform={k:set() for k in EXPECTED};sample_headers=[]
 for index,f in enumerate(files):
  blob=f.read_bytes()
  digest=hashlib.sha256(blob).hexdigest()
  if reference.get(f.name)!=digest:raise ValueError(f"Source checksum differs from the validated acquisition: {f.name}")
  with zipfile.ZipFile(io.BytesIO(blob)) as z:
   h=E.fromstring(z.read("Data0/MesurementConditions0.xml"))
   dt={p[0].text:p[1].text for p in h.findall("ImageInformation/DTHeader/Pair")}
   ax={a.get("Name"):a.get("Position") for a in h.findall("Axes/Axis")}
   for key,val in EXPECTED.items():
    check_values(key,dt[key],val,f.name);uniform[key].add(dt[key])
   if dt.get("BYTE_ORDER")!="little_endian" or dt.get("Data_type")!="unsigned long int":
    raise ValueError(f"Unsupported byte order or data type in {f.name}")
   if dt.get("PXD_SPATIAL_DISTORTION_TYPE")!="Simple_spatial":raise ValueError(f"Unsupported distortion model in {f.name}")
   phi=float(ax["Phi"]);chi=float(ax["Chi"])
   if chi not in [38.,42.]:raise ValueError(f"Unexpected Chi={chi} in {f.name}")
   check_values("TwoTheta Omega",ax["TwoTheta"]+" "+ax["Omega"],[32.6,16.3],f.name)
   check_values("CRYSTAL_GONIO_VALUES",dt["CRYSTAL_GONIO_VALUES"],[16.3,chi-90,phi],f.name)
   start=float(h.findtext("ScanInformation/Start"));stop=float(h.findtext("ScanInformation/Stop"));speed=float(h.findtext("ScanInformation/Speed"))
   if h.findtext("ScanInformation/Mode")!="CONTINUOUS" or h.findtext("ScanInformation/AxisName")!="Phi":
    raise ValueError(f"Unsupported scan axis or mode in {f.name}")
   if not np.allclose([start,stop,speed],[phi-1,phi+1,24],rtol=0,atol=1e-5):
    raise ValueError(f"Scan window is not continuous Phi +/-1 degree at24deg/min in {f.name}")
   raw=z.read("Data0/Image0.bin")
   if len(raw)!=775*385*4:raise ValueError(f"Unexpected binary byte count in {f.name}")
   a=np.frombuffer(raw,dtype="<u4").reshape(385,775).copy()
   if np.any(a>=2**31):raise ValueError(f"Overflow, invalid sentinel, or unsupported count value in {f.name}")
   frames.append(a)
   manifest.append(dict(file=f.name,sha256=hashlib.sha256(blob).hexdigest(),bytes=len(blob),phi=phi,chi=chi,omega=float(ax["Omega"]),two_theta=float(ax["TwoTheta"]),alpha=float(ax["Alpha"]),scan_start=start,scan_stop=stop,scan_speed=speed,counts=int(a.sum()),max=int(a.max()),exptime_ms=float(dt["PXD_DETECTOR_ACQUISITION_EXPTIME"]),created=dt["PXD_CREATE_DATETIME"]))
   if index in [0,180]:sample_headers.append(dict(dt=dt,axes=ax))
  if (index+1)%60==0:print(f"Read and validated {index+1}/360 frames",flush=True)
 order=sorted(range(360),key=lambda j:(manifest[j]["chi"],manifest[j]["phi"]))
 manifest=[manifest[j] for j in order];frames=np.asarray([frames[j] for j in order])
 phi=np.array([r["phi"] for r in manifest]);chi=np.array([r["chi"] for r in manifest])
 if not np.array_equal(phi,np.tile(np.arange(0,360,2),2)) or not np.array_equal(chi,np.repeat([38.,42.],180)):
  raise ValueError("Incomplete, duplicate, or misassigned Phi/Chi frame sequence")
 np.savez_compressed(out/"raw_frames.npz",counts=frames,phi=phi,chi=chi)
 with (out/"frame_manifest.csv").open("w",newline="",encoding="utf-8") as f:
  w=csv.DictWriter(f,fieldnames=manifest[0]);w.writeheader();w.writerows(manifest)
 (out/"headers.json").write_text(json.dumps(sample_headers,indent=2),encoding="utf-8")
 validation=dict(count=360,geometry_errors=[],uniform_fields={k:sorted(v) for k,v in uniform.items()},validated_scan_window="phi center +/-1 degree, continuous24degrees/min, 5seconds for everyframe",source=str(source),source_pattern=PATTERN)
 (out/"geometry_validation.json").write_text(json.dumps(validation,indent=2),encoding="utf-8")
 return manifest
def use_cache(cache,out,reference):
 required=["raw_frames.npz","frame_manifest.csv","geometry_validation.json","headers.json"]
 for name in required:
  if not (cache/name).is_file():raise FileNotFoundError(f"Required cache file missing: {cache/name}")
 with (cache/"frame_manifest.csv").open(newline="",encoding="utf-8") as f:manifest=list(csv.DictReader(f))
 recorded={x["file"]:x["sha256"] for x in manifest}
 if recorded!=reference:raise ValueError("Cached source manifest does not match the validated acquisition")
 valid=json.loads((cache/"geometry_validation.json").read_text(encoding="utf-8"))
 if valid.get("count")!=360 or valid.get("geometry_errors"):raise ValueError("Cached geometry validation did not pass all360frames")
 for key,expected in EXPECTED.items():
  values=valid["uniform_fields"].get(key)
  if values is None:
   # The first verification snapshot omits these redundant keys; check both retained headers.
   headers=json.loads((cache/"headers.json").read_text(encoding="utf-8"))
   values=[x["dt"][key] for x in headers]
  for value in values:check_values(key,value,expected,"cached geometry snapshot")
 with np.load(cache/"raw_frames.npz") as data:
  if data["counts"].shape!=(360,385,775):raise ValueError("Wrong cached image shape")
  if not np.array_equal(data["phi"],np.tile(np.arange(0,360,2),2)):raise ValueError("Wrong cached phi sequence")
  if not np.array_equal(data["chi"],np.repeat([38.,42.],180)):raise ValueError("Wrong cached chi sequence")
 for name in required:shutil.copy2(cache/name,out/name)
 digest=hashlib.sha256((cache/"raw_frames.npz").read_bytes()).hexdigest()
 (out/"cache_validation.json").write_text(json.dumps({"mode":"local validated raw cache; no live source read","cache_raw_npz_sha256":digest,"source_manifest_matches_reference":True,"source_files_in_manifest":360},indent=2),encoding="utf-8")
 print("Validated local raw cache and360source checksums against the reference manifest; no NAS reads",flush=True)
 return manifest
def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument("--source",type=pathlib.Path,default=pathlib.Path(SOURCE),help="Read-only folder containing the supplied360RASX files")
 parser.add_argument("--output",type=pathlib.Path,help="NEW output directory; an existing path is rejected")
 parser.add_argument("--cache",type=pathlib.Path,help="Optional validated local raw-cache folder; skips all NAS reads")
 args=parser.parse_args();program=pathlib.Path(__file__).resolve().parent
 stamp=datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%d_%H%M%S")
 out=(args.output or (program/f"results_{stamp}")).resolve()
 reference_file=program/"reference_manifest.csv"
 if not reference_file.is_file():raise FileNotFoundError(f"Validated reference manifest is missing: {reference_file}")
 with reference_file.open(newline="",encoding="utf-8") as f:reference={x["file"]:x["sha256"] for x in csv.DictReader(f)}
 if len(reference)!=360:raise ValueError("Reference manifest must contain360unique source files")
 if args.cache is None and not args.source.is_dir():raise FileNotFoundError(f"Source folder unavailable: {args.source}")
 out.mkdir(parents=True,exist_ok=False)
 print("New results directory:",out,flush=True)
 manifest=use_cache(args.cache,out,reference) if args.cache else extract(args.source,out,reference)
 steps=[("reconstruct_poles.py",[str(out)]),("fit_poles.py",[str(out/"deliverable_numeric.npz"),str(out)]),("render_results.py",[str(out)])]
 script_hashes={}
 with (out/"processing.log").open("w",encoding="utf-8") as log:
  for script,extra in steps:
   path=program/script
   if not path.is_file():raise FileNotFoundError(f"Required analysis module is missing: {path}")
   script_hashes[script]=hashlib.sha256(path.read_bytes()).hexdigest()
   print("Running",script,flush=True)
   command=([sys.executable,'--script',script] if getattr(sys,'frozen',False) else [sys.executable,str(path)])
   subprocess.run([*command,*extra],check=True,stdout=log,stderr=subprocess.STDOUT)
 import scipy,matplotlib
 record=dict(utc_created=datetime.datetime.now(datetime.timezone.utc).isoformat(),python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,source=str(args.source),output=str(out),source_count=len(manifest),script_sha256=script_hashes,entry_script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),source_write_operations=0,mode="validated_local_cache" if args.cache else "live_source_hash_and_geometry_verification",reference_manifest_sha256=hashlib.sha256(reference_file.read_bytes()).hexdigest())
 (out/"run_provenance.json").write_text(json.dumps(record,indent=2),encoding="utf-8")
 print("Finished. Open LMO_011_pole_figure.pdf and LMO_011_phi_fits.pdf in",out,flush=True)
if __name__=="__main__":
 try:main()
 except Exception as error:
  print("ANALYSIS STOPPED:",error,file=sys.stderr)
  raise

