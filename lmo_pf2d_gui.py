#!/usr/bin/env python3
"""LMO PF2D Analyzer - read-only Rigaku RASX pole-figure GUI for Windows.
Build: pyinstaller --noconfirm --clean --onefile --console --hide-console hide-early --name LMO-PF2D-Analyzer-v1.3.0
       --add-data "analysis;analysis" lmo_pf2d_gui.py
No Rigaku instrument control, X-ray generation, network/cloud mutations, or raw-file edits.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import pathlib
import queue
import re
import runpy
import subprocess
import sys
import threading
import time
import traceback
import zipfile
import xml.etree.ElementTree as ET
from collections import defaultdict

VERSION="1.3.0"
TITLE="LMO PF2D Analyzer  |  HyPix-3000(H)"
RASX_PATTERN=re.compile(
    r"^(?P<stem>.+)_image_(?P<seq>\d+)_(?P<chi>-?\d+(?:\.\d+)?)_(?P<num>\d+)\.rasx$",re.I)

def bundled(rel: str) -> pathlib.Path:
    root=pathlib.Path(getattr(sys, "_MEIPASS", pathlib.Path(__file__).resolve().parent))
    return root/rel

def choose_rasx_folder(folder: str | pathlib.Path, max_depth: int=3) -> pathlib.Path:
    """Automatically select the single nested raw-RASX folder.

    Only lists directory entries (no large image read). If multiple raw
    folders are present, fail safely instead of processing the wrong batch.
    """
    root=pathlib.Path(folder)
    if not root.is_dir():
        raise FileNotFoundError(f"Folder unavailable: {root}")
    stack=[(root,0)]
    candidates=[]
    visited=0
    while stack:
        current,depth=stack.pop()
        visited+=1
        if visited>150:
            raise ValueError("Too many subfolders. Select a folder closer to the RASX data.")
        children=[]
        found=False
        try:
            with os.scandir(current) as items:
                for entry in items:
                    if entry.is_file(follow_symlinks=False) and entry.name.lower().endswith(".rasx"):
                        found=True
                        break
                    if depth<max_depth and entry.is_dir(follow_symlinks=False):
                        name=entry.name.lower()
                        if not (name.startswith(("lmo_pf2d_batch_","lmo_pf2d_reports_")) or
                                name in ("__pycache__",".git")):
                            children.append(pathlib.Path(entry.path))
        except PermissionError:
            if current==root:raise
            continue
        if found:
            candidates.append(current)
        elif depth<max_depth:
            stack.extend((child,depth+1) for child in sorted(children,reverse=True))
    if not candidates:
        raise ValueError(f"No .rasx files found within {max_depth} subfolder levels: {root}")
    if len(candidates)>1:
        names="; ".join(str(p) for p in sorted(candidates)[:4])
        raise ValueError("Multiple RASX folders found. Select one exact measurement folder: "+names)
    return candidates[0]


def choose_existing_batch(folder: str | pathlib.Path, max_depth: int=2) -> pathlib.Path:
    """Select a unique prior result directory containing inspection.json and saved maps."""
    root=pathlib.Path(folder)
    if not root.is_dir():
        raise FileNotFoundError(f"Folder unavailable: {root}")
    stack=[(root,0)]
    candidates=[]
    visited=0
    while stack:
        current,depth=stack.pop()
        visited+=1
        if visited>150:
            raise ValueError("Too many results subfolders. Select a specific LMO_PF2D_BATCH_... folder.")
        if (current/"inspection.json").is_file() and any(current.glob("point_*/map.npz")):
            candidates.append(current)
            continue
        if depth>=max_depth:continue
        try:
            children=[pathlib.Path(e.path) for e in os.scandir(current)
                      if e.is_dir(follow_symlinks=False)
                      and e.name.lower() not in ("__pycache__", ".git")]
        except PermissionError:
            if current==root:raise
            continue
        stack.extend((child,depth+1) for child in sorted(children,reverse=True))
    if not candidates:
        raise ValueError("No saved batch containing inspection.json and point_*/map.npz was found.")
    if len(candidates)>1:
        raise ValueError("Multiple saved batches found. Select the specific LMO_PF2D_BATCH_... folder.")
    return candidates[0]


def ensure_output_parent(folder: str | pathlib.Path) -> pathlib.Path:
    """Create the chosen report parent, including absent intermediate folders."""
    raw=str(folder).strip().strip('"')
    if not raw:
        raise ValueError("Specify a parent folder for reports.")
    target=pathlib.Path(raw)
    target.mkdir(parents=True,exist_ok=True)
    if not target.is_dir():
        raise NotADirectoryError(f"Output path is not a folder: {target}")
    return target

def format_scan_progress(state: dict, elapsed_s: float) -> str:
    """English diagnostic text for optional metadata-only command-line mode."""
    stage=state.get("stage","init")
    seconds=int(max(0,elapsed_s))
    found=state.get("found",0)
    done=state.get("checked",0)
    total=state.get("total")
    accepted=state.get("accepted",0)
    errors=state.get("rejected",0)
    current=state.get("current","")
    last_ok=state.get("last_ok","")
    if stage=="init":
        message="Preparing analyzer and reading directory"
    elif stage=="enumerate":
        message=f"Finding RASX files: {found} found"
    elif stage=="sort":
        message=f"Sorting: {found} RASX files"
    elif stage=="headers":
        message=f"Reading frame headers: {done}/{total}; accepted {accepted}; errors {errors}"
    elif stage in ("group","group_done"):
        message=f"Grouping XY positions: accepted {accepted}; errors {errors}"
    elif stage=="validate":
        message=(f"Validating positions: {state.get('groups_done',0)}/"
                 f"{state.get('groups_total',0)}; ready {state.get('groups_ok',0)}")
    else:
        message="Checking RASX metadata"
    message+=f" | {seconds} s"
    if stage=="headers" and current:message+=" | current: "+current[-90:]
    if stage=="headers" and last_ok:message+=" | last OK: "+last_ok[-90:]
    return message

def scan_groups(folder: str | pathlib.Path, progress=None) -> list[dict]:
    """Infer scan mode and frame counts; optionally report progress snapshots."""
    lib=runpy.run_path(str(bundled("analysis/adaptive_pf2d.py")),run_name="_adaptive_")
    last_state={}
    def report(state):
        nonlocal last_state
        last_state=state
        if progress is not None:progress(state)
    groups,failures=lib["discover"](folder,progress=report if progress else None)
    if not groups and failures:
        raise ValueError("No calibrated detector frames; first error: "+failures[0]["error"])
    rows=[]
    completed=0
    ready=0
    total_groups=len(groups)
    if progress is not None:
        report(dict(last_state,stage="validate",current="",
                    groups_done=0,groups_total=total_groups,groups_ok=0))
    for key,items in groups.items():
        try:
            bychi,info=lib["validate_group"](items)
            ok=True
            expected=sum(v.get("expected_unique_frames") or 0 for v in info.values())
            mode=", ".join(f"{c}: {v['phi_step_deg']:.2f}° {v['mode']}"
                           for c,v in info.items())
        except Exception as e:
            ok=False
            expected=len(items)
            mode="CHECK: "+str(e)[:120]
        x=items[0].get("x")
        y=items[0].get("y")
        xy=f"X={x}, Y={y} mm | {mode}"
        rows.append(dict(key=str(key),frames=len(items),expected=expected,ready=ok,
                         chi38=sum(abs(v["chi"]-38)<.01 for v in items),
                         chi42=sum(abs(v["chi"]-42)<.01 for v in items),x_y=xy))
        completed+=1
        ready+=int(ok)
        if progress is not None:
            report(dict(last_state,stage="validate",current="",
                        groups_done=completed,groups_total=total_groups,groups_ok=ready))
    return rows

def run_script_mode(script: str, args: list[str]):
    allowed={"run_lmo_pf.py","reconstruct_poles.py","fit_poles.py","render_results.py"}
    if script not in allowed:
        raise ValueError("Unsupported worker script: "+script)
    if len(args)>=2 and args[0]=='--stage-root':
        target=pathlib.Path(args[1])/script
        args=args[2:]
    else:
        target=bundled('analysis/validated_single_point_scripts/'+script)
    if not target.is_file():
        raise FileNotFoundError("Bundled worker not found: "+str(target))
    sys.argv=[str(target),*args]
    runpy.run_path(str(target),run_name="__main__")

def run_analyze_mode(args: list[str]):
    target=bundled("analysis/adaptive_pf2d.py")
    if not target.is_file():
        raise FileNotFoundError("Bundled analyzer not found: "+str(target))
    sys.argv=[str(target),*args]
    runpy.run_path(str(target),run_name="__main__")

def run_legacy_mode(args: list[str]):
    target=bundled("analysis/analyze_many_xy.py")
    sys.argv=[str(target),*args]
    runpy.run_path(str(target),run_name="__main__")

def make_gui():
    import tkinter as tk
    from tkinter import filedialog, messagebox, scrolledtext, ttk

    class LmoApp(tk.Tk):
        def __init__(self):
            super().__init__()
            self.title(TITLE)
            self.geometry("1120x780")
            self.minsize(920,660)
            self.configure(background="#f3f6fa")
            self.source=tk.StringVar()
            self.output_parent=tk.StringVar(value=str(pathlib.Path.home()/"Documents"/"LMO PF2D Reports"))
            self.mode=tk.StringVar(value="batch")
            self.position=tk.StringVar(value="1")
            self.info=tk.StringVar(value="Select a source folder, then click Analyze & Generate Reports. No separate file check required.")
            self.progress=tk.StringVar(value="Ready")
            self.worker=None
            self.launching=False
            self.events=queue.Queue()
            self.report_dir=None
            self.rows={}
            self.current_stage="Ready"
            self.last_progress=None
            self.started_at=0.
            self.last_heartbeat=0.
            self._create_ui(ttk,scrolledtext)
            self.after(180,self._poll)

        def _create_ui(self,ttk,scrolledtext):
            style=ttk.Style(self)
            if "clam" in style.theme_names():
                style.theme_use("clam")
            style.configure("TFrame",background="#f3f6fa")
            style.configure("Header.TLabel",background="#f3f6fa",foreground="#1d354c",font=("Segoe UI",19,"bold"))
            style.configure("Text.TLabel",background="#f3f6fa",foreground="#344f69",font=("Segoe UI",10))
            style.configure("Small.TLabel",background="#f3f6fa",foreground="#63758a",font=("Segoe UI",9))
            style.configure("Action.TButton",font=("Segoe UI",10,"bold"),padding=(16,9))
            style.configure("TButton",font=("Segoe UI",9),padding=(9,6))
            style.configure("TEntry",padding=6)
            style.configure("TLabel",background="#f3f6fa",font=("Segoe UI",10))
            style.configure("TRadiobutton",background="#f3f6fa",font=("Segoe UI",10))
            style.configure("Treeview",font=("Segoe UI",10),rowheight=26)
            style.configure("Treeview.Heading",font=("Segoe UI",10,"bold"))
            body=ttk.Frame(self,padding=19)
            body.pack(fill="both",expand=True)
            header=ttk.Frame(body)
            header.pack(fill="x")
            ttk.Label(header,text="LMO | 2D Pole Figure Analyzer",style="Header.TLabel").pack(side="left")
            ttk.Label(header,text="v"+VERSION+" | HyPix-3000(H)",style="Small.TLabel").pack(side="left",padx=(14,0),pady=(7,0))
            ttk.Label(body,text="One-pass workflow: check RASX metadata -> reconstruct 2D maps -> fit four peaks -> PDF + CSV. Raw files are read-only.",
                      style="Text.TLabel").pack(anchor="w",pady=(5,17))

            srcrow=ttk.Frame(body)
            srcrow.pack(fill="x",pady=3)
            ttk.Label(srcrow,text="Input RASX folder",width=19).pack(side="left")
            ttk.Entry(srcrow,textvariable=self.source).pack(side="left",fill="x",expand=True,padx=(3,8))
            ttk.Button(srcrow,text="Browse...",command=self._browse_source).pack(side="left")

            outrow=ttk.Frame(body)
            outrow.pack(fill="x",pady=7)
            ttk.Label(outrow,text="Output parent folder",width=19).pack(side="left")
            ttk.Entry(outrow,textvariable=self.output_parent).pack(side="left",fill="x",expand=True,padx=(3,8))
            ttk.Button(outrow,text="Browse...",command=self._browse_output).pack(side="left")
            ttk.Label(body,text="The output parent may be new. A timestamped results folder is created automatically; existing results are never overwritten.",
                      style="Small.TLabel").pack(anchor="w",pady=(1,11))

            modes=ttk.Frame(body)
            modes.pack(fill="x",pady=(0,12))
            ttk.Label(modes,text="Analysis scope",width=19).pack(side="left")
            ttk.Radiobutton(modes,text="All measured positions",variable=self.mode,value="batch").pack(side="left",padx=(3,18))
            ttk.Radiobutton(modes,text="One position #",variable=self.mode,value="single").pack(side="left",padx=(0,5))
            ttk.Entry(modes,textvariable=self.position,width=6).pack(side="left",padx=(0,16))
            ttk.Label(modes,text="Position numbers come from the RASX acquisition suffix (1-39).",
                      style="Small.TLabel").pack(side="left")

            actions=ttk.Frame(body)
            actions.pack(fill="x",pady=(2,11))
            self.start_button=ttk.Button(actions,text="Analyze & Generate Reports",style="Action.TButton",command=self._start)
            self.start_button.pack(side="left",padx=(0,8))
            self.rebuild_button=ttk.Button(actions,text="Rebuild PDFs from Saved Maps",command=self._rebuild_reports)
            self.rebuild_button.pack(side="left",padx=(0,8))
            self.cancel_button=ttk.Button(actions,text="Stop",command=self._cancel,state="disabled")
            self.cancel_button.pack(side="left",padx=(0,8))
            ttk.Button(actions,text="Open Results",command=self._open_report).pack(side="left")
            ttk.Label(actions,textvariable=self.progress,style="Text.TLabel").pack(side="right",padx=6)

            ttk.Label(body,textvariable=self.info,style="Small.TLabel",wraplength=1050).pack(anchor="w",pady=(0,9))
            ttk.Label(body,text="Measurement positions / live analysis status",style="Text.TLabel").pack(anchor="w",pady=(0,5))
            frame=ttk.Frame(body)
            frame.pack(fill="both",expand=True)
            cols=("position","frames","chi38","chi42","xy","status")
            self.tree=ttk.Treeview(frame,columns=cols,show="headings",height=10,selectmode="browse")
            for name,heading,width in [
                ("position","Position",85),("frames","Frames",90),("chi38","Chi 38 deg",95),
                ("chi42","Chi 42 deg",95),("xy","Stage X / Y (mm)",220),("status","Processing",180)]:
                self.tree.heading(name,text=heading)
                self.tree.column(name,width=width,minwidth=70,anchor="w" if name in ("xy","status") else "center")
            sy=ttk.Scrollbar(frame,orient="vertical",command=self.tree.yview)
            sx=ttk.Scrollbar(frame,orient="horizontal",command=self.tree.xview)
            self.tree.configure(yscrollcommand=sy.set,xscrollcommand=sx.set)
            self.tree.grid(row=0,column=0,sticky="nsew")
            sy.grid(row=0,column=1,sticky="ns")
            sx.grid(row=1,column=0,sticky="ew")
            frame.columnconfigure(0,weight=1);frame.rowconfigure(0,weight=1)

            ttk.Label(body,text="Processing log",style="Text.TLabel").pack(anchor="w",pady=(12,5))
            self.log=scrolledtext.ScrolledText(body,height=8,font=("Consolas",9),
                state="disabled",background="#101c2c",foreground="#d7ebfb",insertbackground="#ffffff",wrap="word")
            self.log.pack(fill="both",expand=False)
            self._log("Ready. Select a source and click Analyze & Generate Reports.")
            self._log("The same run discovers frames, verifies geometry, reconstructs maps and fits peaks.")
            self.protocol("WM_DELETE_WINDOW",self._close)

        def _log(self,line):
            self.log.configure(state="normal")
            self.log.insert("end",str(line)+"\n")
            self.log.see("end")
            self.log.configure(state="disabled")

        def _browse_source(self):
            folder=filedialog.askdirectory(title="Select a RASX folder (or its parent)",
                                           initialdir=self.source.get() if pathlib.Path(self.source.get()).is_dir() else str(pathlib.Path.home()))
            if folder:self.source.set(folder)

        def _browse_output(self):
            candidate=pathlib.Path(self.output_parent.get())
            base=candidate if candidate.is_dir() else candidate.parent
            folder=filedialog.askdirectory(title="Select a parent folder for results",
                                           initialdir=str(base) if base.is_dir() else str(pathlib.Path.home()))
            if folder:self.output_parent.set(folder)

        def _busy(self,on):
            self.launching=on
            self.start_button.configure(state="disabled" if on else "normal")
            self.rebuild_button.configure(state="disabled" if on else "normal")
            self.cancel_button.configure(state="normal" if on and self.worker is not None else "disabled")

        def _destination(self,parent,prefix):
            stamp=dt.datetime.now().strftime("%Y%m%d_%H%M%S")
            path=parent/f"{prefix}_{stamp}"
            n=1
            while path.exists():
                path=parent/f"{prefix}_{stamp}_{n:02d}"
                n+=1
            return path

        def _start(self):
            if self.launching or self.worker is not None:return
            src=self.source.get().strip().strip('"')
            if not src:
                messagebox.showwarning("Input required","Select the input RASX folder.");return
            position=None
            if self.mode.get()=="single":
                try:
                    position=int(self.position.get())
                    if position<1:raise ValueError
                except ValueError:
                    messagebox.showwarning("Invalid position","Enter a positive numeric measurement position (e.g. 37).")
                    return
            self._prepare_new_job()
            self._log("Starting single-pass raw RASX analysis...")
            self.current_stage="Locating RASX folder"
            self.info.set(self.current_stage)
            def prepare():
                try:
                    raw=choose_rasx_folder(src)
                    parent=ensure_output_parent(self.output_parent.get())
                    scope="ONE" if position is not None else "BATCH"
                    dest=self._destination(parent,"LMO_PF2D_"+scope)
                    args=["--source",str(raw),"--output",str(dest)]
                    if position is not None:args+=["--position",str(position)]
                    cmd=([sys.executable,"--analyze",*args] if getattr(sys,"frozen",False) else
                         [sys.executable,str(pathlib.Path(__file__).resolve()),"--analyze",*args])
                    self.events.put(("launch",(cmd,dest,
                       f"Raw RASX single-pass analysis: {raw}",str(raw))))
                except Exception as exc:
                    self.events.put(("prepare_error",str(exc)))
            threading.Thread(target=prepare,daemon=True).start()

        def _rebuild_reports(self):
            if self.launching or self.worker is not None:return
            directory=filedialog.askdirectory(
                title="Select an earlier LMO_PF2D_BATCH folder containing saved maps",
                initialdir=self.source.get() if pathlib.Path(self.source.get()).is_dir() else str(pathlib.Path.home()))
            if not directory:return
            self._prepare_new_job()
            self.current_stage="Locating saved 2D maps"
            self.info.set(self.current_stage)
            def prepare():
                try:
                    source=choose_existing_batch(directory)
                    parent=ensure_output_parent(self.output_parent.get())
                    dest=self._destination(parent,"LMO_PF2D_REPORTS")
                    args=["--source",str(source),"--output",str(dest),"--reports-only"]
                    cmd=([sys.executable,"--analyze",*args] if getattr(sys,"frozen",False) else
                         [sys.executable,str(pathlib.Path(__file__).resolve()),"--analyze",*args])
                    self.events.put(("launch",(cmd,dest,
                        f"Rebuild reports from saved maps: {source} (no RASX read)",None)))
                except Exception as exc:
                    self.events.put(("prepare_error",str(exc)))
            threading.Thread(target=prepare,daemon=True).start()

        def _prepare_new_job(self):
            self._busy(True)
            self.started_at=time.monotonic()
            self.last_heartbeat=self.started_at
            self.last_progress=None
            self.progress.set("Working...")
            self.report_dir=None
            for iid in self.tree.get_children():self.tree.delete(iid)
            self.rows={}

        def _launch_worker(self,cmd,dest,description,resolved_source):
            self.report_dir=dest
            if resolved_source:self.source.set(resolved_source)
            self._log("\nSTART: "+str(dest))
            self._log(description)
            def task():
                try:
                    flags=getattr(subprocess,"CREATE_NO_WINDOW",0)
                    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                        text=True,encoding="utf-8",errors="replace",bufsize=1,creationflags=flags)
                    self.worker=proc
                    self.events.put(("worker_started",None))
                    for line in proc.stdout:
                        line=line.rstrip()
                        if line.startswith("PROGRESS_JSON:"):
                            try:self.events.put(("progress",json.loads(line[len("PROGRESS_JSON:"):])))
                            except Exception:self.events.put(("log",line))
                        else:
                            self.events.put(("log",line))
                    self.events.put(("done",proc.wait()))
                except Exception:
                    self.events.put(("worker_error",traceback.format_exc()))
            threading.Thread(target=task,daemon=True).start()

        def _point_row(self,point):
            if point not in self.rows:
                self.rows[point]=[point,"","","","","Queued"]
                self.tree.insert("","end",iid=str(point),values=self.rows[point])
            return self.rows[point]

        def _update_point(self,point,**updates):
            if point is None:return
            point=int(point)
            row=self._point_row(point)
            for key,val in updates.items():
                pos={"frames":1,"chi38":2,"chi42":3,"xy":4,"status":5}.get(key)
                if pos is not None:row[pos]=val
            self.tree.item(str(point),values=row)
            if self.tree.exists(str(point)):
                self.tree.see(str(point))

        def _progress_event(self,p):
            stage=p.get("stage","")
            self.last_progress=p
            message=""
            if stage=="start":
                message="Checking RASX metadata (same analysis run)"
            elif stage=="file_scan":
                phase=p.get("scan_phase","")
                if phase=="enumerate":
                    message=f"Finding RASX files: {p.get('found',0)} found"
                elif phase=="headers":
                    message=(f"Reading frame headers: {p.get('checked',0)} / {p.get('total','?')} "
                             f"| accepted {p.get('accepted',0)} | errors {p.get('rejected',0)}")
                    if p.get("current"):message+=" | current: "+str(p['current'])[-45:]
                else:message=f"Metadata stage: {phase}"
            elif stage=="metadata_done":
                message=(f"Metadata collected: {p.get('groups',0)} positions; "
                         f"{p.get('valid_files',0)} frames; {p.get('errors',0)} errors. "
                         "Starting image reconstruction and peak fitting.")
            elif stage=="point_found":
                point=p.get("point")
                self._update_point(point,frames=p.get("frames",0),chi38=p.get("chi38",0),
                                   chi42=p.get("chi42",0),
                                   xy=f"X={p.get('x')} / Y={p.get('y')}",status="Queued")
                return
            elif stage in ("point_start","reconstruct","fit_and_report"):
                point=p.get("point")
                state={"point_start":"Validating geometry","reconstruct":"Reconstructing 2D map",
                       "fit_and_report":"Fitting peaks and writing PDFs"}[stage]
                message=f"Position {p.get('index','?')}/{p.get('total','?')}: {state}"
                self._update_point(point,status=state)
            elif stage=="frame_progress":
                point=p.get("point")
                done=p.get("frames_done",0);total=p.get("frames_total",0)
                message=(f"Position {p.get('index','?')}/{p.get('total','?')}: "
                         f"reconstructed {done}/{total} RASX frames")
                if p.get("current_file"):message+=" | last: "+str(p["current_file"])[-48:]
                self._update_point(point,status=f"Frames {done}/{total}")
            elif stage=="point_done":
                point=p.get("point")
                self._update_point(point,status=p.get("status","UNKNOWN"))
                message=(f"Position {p.get('index','?')}/{p.get('total','?')} finished: "
                         f"{p.get('status','UNKNOWN')}")
                if p.get('reason'):message+=" | "+str(p["reason"])[:120]
            elif stage=="summary":
                message="Generating batch comparison plots, tables and PDF"
            elif stage=="complete":
                message=(f"Analysis finished: {p.get('successful',0)}/{p.get('total',0)} "
                         f"successful; {p.get('failed',0)} failed")
            else:return
            self.current_stage=message
            self.info.set(message)
            self.progress.set("Working..." if stage!="complete" else "Completed")
            self._log(message)

        def _heartbeat(self):
            if not self.launching:return
            now=time.monotonic()
            if now-self.last_heartbeat<5.0:return
            self.last_heartbeat=now
            elapsed=int(now-self.started_at)
            self.info.set(self.current_stage+f" | elapsed {elapsed} s")
            self._log(f"Progress heartbeat [{elapsed} s]: {self.current_stage}")

        def _cancel(self):
            proc=self.worker
            if proc is None or proc.poll() is not None:
                if self.launching:messagebox.showinfo("Initializing","The input/output folders are being prepared.")
                return
            if not messagebox.askyesno("Stop analysis","Stop the current analysis? Previously created results are retained."):
                return
            self._log("Stopping analysis at user request...")
            try:
                if os.name=="nt":
                    subprocess.run(["taskkill","/PID",str(proc.pid),"/T","/F"],
                                   stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,
                                   creationflags=getattr(subprocess,"CREATE_NO_WINDOW",0))
                else:proc.terminate()
            except Exception as exc:self._log("Stop error: "+str(exc))

        def _open_report(self):
            target=self.report_dir
            if target is None or not target.exists():
                messagebox.showinfo("No reports yet","Start an analysis first.");return
            pdf=target/"LMO011_XY_Overall_Report.pdf"
            if pdf.is_file():target=pdf
            try:
                if os.name=="nt":os.startfile(str(target))
                elif sys.platform=="darwin":subprocess.Popen(["open",str(target)])
                else:subprocess.Popen(["xdg-open",str(target)])
            except Exception as exc:messagebox.showerror("Open results failed",str(exc))

        def _finish(self,rc):
            self.worker=None
            self._busy(False)
            self.cancel_button.configure(state="disabled")
            self.progress.set("Completed" if rc==0 else f"Exit code {rc}")
            target=self.report_dir
            self._log(("SUCCESS" if rc==0 else "FINISHED WITH WARNINGS / ERRORS")+": "+str(target))
            if rc==0 and target:
                messagebox.showinfo("Reports created","Reports are saved in:\n"+str(target))
            elif target:
                self.info.set(f"Analysis stopped with exit code {rc}; inspect the processing log and {target}")

        def _poll(self):
            for _ in range(120):
                try:kind,payload=self.events.get_nowait()
                except queue.Empty:break
                if kind=="launch":
                    self._launch_worker(*payload)
                elif kind=="worker_started":
                    self.cancel_button.configure(state="normal")
                elif kind=="prepare_error":
                    self._log("Cannot start analysis: "+payload)
                    self.info.set("Cannot start: "+payload)
                    self.progress.set("Error")
                    self._busy(False)
                    messagebox.showerror("Cannot start analysis",payload)
                elif kind=="progress":
                    self._progress_event(payload)
                elif kind=="log":
                    self._log(payload)
                    match=re.search(r"^REPORT\s+(\d+)/(\d+):\s+point\s+(\d+)\s+(\w+)",payload)
                    if match:
                        n,total,point,status=match.groups()
                        self.current_stage=f"Rebuilding saved maps: {n}/{total} | point {point} {status}"
                        self.info.set(self.current_stage)
                elif kind=="done":
                    self._finish(payload)
                elif kind=="worker_error":
                    self._log(payload)
                    self._busy(False)
                    self.progress.set("Error")
                    self.info.set("Worker launch failed. See processing log.")
            self._heartbeat()
            self.after(180,self._poll)

        def _close(self):
            if self.launching or (self.worker is not None and self.worker.poll() is None):
                messagebox.showwarning("Analysis running","Stop the analysis or let it finish before closing.")
                return
            self.destroy()

    return LmoApp()

def main():
    if len(sys.argv)>1:
        mode=sys.argv[1]
        if mode=="--script":
            if len(sys.argv)<3:raise SystemExit("--script expects a module filename")
            run_script_mode(sys.argv[2],sys.argv[3:]);return
        if mode=="--analyze":
            run_analyze_mode(sys.argv[2:]);return
        if mode=="--inspect":
            ap=argparse.ArgumentParser()
            ap.add_argument("--source",required=True)
            args=ap.parse_args(sys.argv[2:])
            print(json.dumps(scan_groups(args.source),indent=2,ensure_ascii=False))
            return
        if mode=="--smoke-gui":
            if len(sys.argv)!=3:raise SystemExit("Use --smoke-gui path-to-result.json")
            dest=pathlib.Path(sys.argv[2])
            try:
                app=make_gui()
                app.update_idletasks()
                app.update()
                info={"status":"OK","title":app.title(),"widgets":len(app.winfo_children()),
                    "version":VERSION,"tk_version":app.tk.call("info","patchlevel")}
                app.destroy()
            except Exception as exc:
                info={"status":"FAILED","error":repr(exc),"traceback":traceback.format_exc()}
            dest.parent.mkdir(parents=True,exist_ok=True)
            dest.write_text(json.dumps(info,indent=2),encoding="utf-8")
            if info["status"]!="OK":raise SystemExit(1)
            return
    app=make_gui()
    app.mainloop()

if __name__=="__main__":
    main()
