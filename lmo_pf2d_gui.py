#!/usr/bin/env python3
"""LMO PF2D Analyzer - read-only Rigaku RASX pole-figure GUI for Windows.
Build: pyinstaller --noconfirm --clean --onefile --windowed --name LMO-PF2D-Analyzer
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

VERSION="1.2.0"
TITLE="LMO PF2D Analyzer  |  HyPix-3000(H)"
RASX_PATTERN=re.compile(
    r"^(?P<stem>.+)_image_(?P<seq>\d+)_(?P<chi>-?\d+(?:\.\d+)?)_(?P<num>\d+)\.rasx$",re.I)

def bundled(rel: str) -> pathlib.Path:
    root=pathlib.Path(getattr(sys, "_MEIPASS", pathlib.Path(__file__).resolve().parent))
    return root/rel

def format_scan_progress(state: dict, elapsed_s: float) -> str:
    """A readable, five-second heartbeat even if a single SMB read blocks."""
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
        message="Подготовка анализатора и чтение каталога"
    elif stage=="enumerate":
        message=f"Поиск RASX: найдено {found} файлов"
    elif stage=="sort":
        message=f"Сортировка списка: {found} RASX"
    elif stage=="headers":
        message=(f"Чтение заголовков: {done}/{total} файлов; "
                 f"прошли {accepted}; ошибки {errors}")
    elif stage in ("group","group_done"):
        message=f"Группировка XY-позиций: прочитано {accepted}, ошибки {errors}"
    elif stage=="validate":
        message=(f"Проверка серий: {state.get('groups_done',0)}/"
                 f"{state.get('groups_total',0)}; "
                 f"готово {state.get('groups_ok',0)}")
    else:
        message="Проверка RASX"
    message+=f" | {seconds} с"
    if stage=="headers" and current:
        message+=" | сейчас: "+current[-90:]
    if stage=="headers" and last_ok:
        message+=" | последний ОК: "+last_ok[-90:]
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
            self.geometry("1050x770")
            self.minsize(850,620)
            self.configure(background="#f3f6fa")
            self.source=tk.StringVar()
            self.output_parent=tk.StringVar(value=str(pathlib.Path.home()/"Documents"))
            self.mode=tk.StringVar(value="single")
            self.backend=tk.StringVar(value="adaptive")
            self.info=tk.StringVar(value="Выберите папку с файлами .rasx и нажмите «Проверить данные».")
            self.progress=tk.StringVar(value="Ожидание")
            self.rows=[]
            self.worker=None
            self.events=queue.Queue()
            self.report_dir=None
            self.scan_active=False
            self.scan_snapshot={"stage":"init"}
            self.scan_started_at=0.0
            self.scan_last_log=0.0
            self._create_ui(ttk,scrolledtext,tk)
            self.after(140,self._poll)

        def _create_ui(self,ttk,scrolledtext,tk):
            style=ttk.Style(self)
            if "clam" in style.theme_names():
                style.theme_use("clam")
            style.configure("TFrame",background="#f3f6fa")
            style.configure("Pane.TFrame",background="#ffffff")
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
            root=ttk.Frame(self,padding=19)
            root.pack(fill="both",expand=True)
            title=ttk.Frame(root)
            title.pack(fill="x")
            ttk.Label(title,text="LMO · 2D Pole Figure Analyzer",style="Header.TLabel").pack(side="left")
            ttk.Label(title,text=" v"+VERSION+"  |  HyPix-3000(H)",style="Small.TLabel").pack(side="left",padx=(12,0),pady=(7,0))
            ttk.Label(root,text="Rigaku RASX → 4 полюса → FWHM φ / β / tilt → PDF + CSV.  Исходные файлы только читаются.",
                      style="Text.TLabel").pack(anchor="w",pady=(6,18))
            source_row=ttk.Frame(root)
            source_row.pack(fill="x",pady=3)
            ttk.Label(source_row,text="Папка RASX",width=17).pack(side="left")
            self.source_entry=ttk.Entry(source_row,textvariable=self.source)
            self.source_entry.pack(side="left",fill="x",expand=True,padx=(4,8))
            ttk.Button(source_row,text="Обзор…",command=self._browse_source).pack(side="left",padx=(0,6))
            self.scan_button=ttk.Button(source_row,text="Проверить данные",command=self._scan)
            self.scan_button.pack(side="left")
            output_row=ttk.Frame(root)
            output_row.pack(fill="x",pady=7)
            ttk.Label(output_row,text="Папка для отчётов",width=17).pack(side="left")
            ttk.Entry(output_row,textvariable=self.output_parent).pack(side="left",fill="x",expand=True,padx=(4,8))
            ttk.Button(output_row,text="Обзор…",command=self._browse_output).pack(side="left")
            ttk.Label(root,text="При запуске создаётся НОВАЯ подпапка с датой и временем. Готовые отчёты не перезаписываются.",
                      style="Small.TLabel").pack(anchor="w",pady=(1,12))
            mode_row=ttk.Frame(root)
            mode_row.pack(fill="x",pady=(0,9))
            ttk.Label(mode_row,text="Обработка",width=17).pack(side="left")
            ttk.Radiobutton(mode_row,text="Одна выбранная точка",variable=self.mode,value="single").pack(side="left",padx=(4,20))
            ttk.Radiobutton(mode_row,text="Все точки (пакет)",variable=self.mode,value="batch").pack(side="left")
            ttk.Radiobutton(mode_row,text="Авто: step / continuous, шаг из RASX",
                            variable=self.backend,value="adaptive").pack(side="left",padx=(4,18))
            ttk.Radiobutton(mode_row,text="Legacy 2°",variable=self.backend,
                            value="legacy").pack(side="left")
            ttk.Label(root,text="Обнаруженные точки / серии",style="Text.TLabel").pack(anchor="w",pady=(5,5))
            frame=ttk.Frame(root)
            frame.pack(fill="both",expand=True)
            cols=("frames","chi38","chi42","xy","status")
            self.tree=ttk.Treeview(frame,columns=cols,show="tree headings",height=9,selectmode="browse")
            self.tree.heading("#0",text="Серия (папка + префикс)")
            self.tree.column("#0",minwidth=240,width=365)
            for key,cap,width in [("frames","Кадры",80),("chi38","χ 38°",75),("chi42","χ 42°",75),
                                  ("xy","XY из RASX",195),("status","Статус",110)]:
                self.tree.heading(key,text=cap);self.tree.column(key,width=width,minwidth=60,anchor="center" if key!="xy" else "w")
            sy=ttk.Scrollbar(frame,orient="vertical",command=self.tree.yview)
            sx=ttk.Scrollbar(frame,orient="horizontal",command=self.tree.xview)
            self.tree.configure(yscrollcommand=sy.set,xscrollcommand=sx.set)
            self.tree.grid(row=0,column=0,sticky="nsew")
            sy.grid(row=0,column=1,sticky="ns")
            sx.grid(row=1,column=0,sticky="ew")
            frame.columnconfigure(0,weight=1);frame.rowconfigure(0,weight=1)
            ttk.Label(root,textvariable=self.info,style="Small.TLabel",wraplength=970).pack(anchor="w",pady=(8,8))
            action=ttk.Frame(root)
            action.pack(fill="x",pady=(0,8))
            self.start_button=ttk.Button(action,text="▶  Построить отчёт",style="Action.TButton",command=self._start)
            self.start_button.pack(side="left",padx=(0,9))
            self.rebuild_button=ttk.Button(action,text="Новые PDF из готового расчёта",command=self._rebuild_reports)
            self.rebuild_button.pack(side="left",padx=(0,9))
            self.cancel_button=ttk.Button(action,text="Остановить",command=self._cancel,state="disabled")
            self.cancel_button.pack(side="left",padx=(0,9))
            ttk.Button(action,text="Открыть результаты",command=self._open_report).pack(side="left")
            ttk.Label(action,textvariable=self.progress,style="Text.TLabel").pack(side="right",padx=8)
            ttk.Label(root,text="Журнал обработки",style="Text.TLabel").pack(anchor="w",pady=(7,5))
            self.log=scrolledtext.ScrolledText(root,height=8,font=("Consolas",9),state="disabled",
                 background="#101c2c",foreground="#d7ebfb",insertbackground="#ffffff",wrap="word")
            self.log.pack(fill="both",expand=False)
            self._log("Программа не управляет дифрактометром. Авто-режим читает углы и число кадров из RASX.")
            self.protocol("WM_DELETE_WINDOW",self._close)

        def _log(self,line):
            self.log.configure(state="normal")
            self.log.insert("end",str(line)+"\n")
            self.log.see("end")
            self.log.configure(state="disabled")

        def _browse_source(self):
            chosen=filedialog.askdirectory(title="Папка с RASX",initialdir=self.source.get() or str(pathlib.Path.home()))
            if chosen:
                self.source.set(chosen);self._scan()

        def _browse_output(self):
            chosen=filedialog.askdirectory(title="Родительская папка отчётов",initialdir=self.output_parent.get())
            if chosen:self.output_parent.set(chosen)

        def _scan(self):
            if self.worker is not None or self.scan_active:return
            source=self.source.get().strip().strip('"')
            if not source:
                messagebox.showwarning("Папка не задана","Укажите папку с RASX.")
                return
            self.rows=[]
            self.scan_active=True
            self.scan_snapshot={"stage":"init"}
            self.scan_started_at=time.monotonic()
            self.scan_last_log=self.scan_started_at
            self.info.set("Начало проверки: поиск RASX и чтение метаданных…")
            self.progress.set("Проверка RASX…")
            self.scan_button.configure(state="disabled")
            self.start_button.configure(state="disabled")
            self.rebuild_button.configure(state="disabled")
            self._log("ПРОВЕРКА НАЧАТА: исходные RASX только читаются; статус каждые 5 с.")
            def on_progress(snapshot):
                # Worker changes a plain immutable snapshot. Only the GUI thread touches Tk.
                self.scan_snapshot=snapshot
            def scan_task():
                try:self.events.put(("scan_ok",scan_groups(source,progress=on_progress)))
                except Exception as exc:self.events.put(("scan_error",str(exc)))
            threading.Thread(target=scan_task,daemon=True).start()

        def _scan_heartbeat(self):
            if not self.scan_active:return
            now=time.monotonic()
            if now-self.scan_last_log<5.0:return
            self.scan_last_log=now
            state=self.scan_snapshot
            message=format_scan_progress(state,now-self.scan_started_at)
            self.info.set(message)
            self._log(message)
            if state.get("stage")=="headers":
                self.progress.set(f"RASX {state.get('checked',0)}/{state.get('total','?')}")
            elif state.get("stage")=="validate":
                self.progress.set(f"Точки {state.get('groups_done',0)}/{state.get('groups_total','?')}")
            else:
                self.progress.set("Проверка идёт…")

        def _display_scan(self,items):
            self.rows=items
            for row_id in self.tree.get_children():self.tree.delete(row_id)
            for i,r in enumerate(items):
                self.tree.insert("",'end',iid=str(i),text=r["key"],values=(
                    r["frames"],r["chi38"],r["chi42"],r["x_y"],"Готово" if r["ready"] else "Неполно"))
            if items:self.tree.selection_set("0")
            ready=sum(r["ready"] for r in items)
            self.info.set(f"Найдено серий: {len(items)}. Проверено по метаданным: {ready}. "
                          "Для контроля геометрии каждого кадра требуется запуск анализа.")
            self._log(self.info.get())

        def _start(self):
            if self.worker is not None or self.scan_active:return
            src=pathlib.Path(self.source.get().strip().strip('"'))
            parent=pathlib.Path(self.output_parent.get().strip().strip('"'))
            if not src.is_dir():
                messagebox.showerror("Ошибка","Папка с RASX недоступна.");return
            if not parent.is_dir():
                messagebox.showerror("Ошибка","Папка отчётов должна существовать.");return
            if not self.rows:
                messagebox.showwarning("Нет серий","Сначала нажмите «Проверить данные».");return
            args=["--source",str(src)]
            kind=self.mode.get()
            if kind=="single":
                selected=self.tree.selection()
                if not selected:
                    messagebox.showwarning("Выбор","Выберите строку с нужной точкой.");return
                chosen=self.rows[int(selected[0])]
                if not chosen["ready"]:
                    messagebox.showwarning("Неполная серия","Данные не прошли проверку углов/метаданных.");return
                args+=["--group-key",chosen["key"]]
            elif any(not r["ready"] for r in self.rows):
                if not messagebox.askyesno("Неполные серии",
                        "Часть серий содержит неполное число кадров. Отчёт отобразит их как FAILED. Продолжить?"):
                    return
            stamp=dt.datetime.now().strftime("%Y%m%d_%H%M%S")
            dest=parent/f"LMO_PF2D_{'ONE' if kind=='single' else 'BATCH'}_{stamp}"
            n=1
            while dest.exists():
                dest=parent/f"LMO_PF2D_{'ONE' if kind=='single' else 'BATCH'}_{stamp}_{n:02d}"
                n+=1
            flag="--analyze" if self.backend.get()=="adaptive" else "--analyze-legacy"
            cmd=[sys.executable,flag,*args,"--output",str(dest)] if getattr(sys,"frozen",False) else [
                sys.executable,str(pathlib.Path(__file__).resolve()),flag,*args,"--output",str(dest)]
            self._launch_worker(cmd,dest,"Исходный набор RASX: "+("одна точка" if kind=="single" else "все точки"))

        def _rebuild_reports(self):
            """Reconstruct rich PDFs from saved maps in an existing batch result folder."""
            if self.worker is not None or self.scan_active:return
            source=filedialog.askdirectory(title="Выберите готовую папку LMO_PF2D_BATCH_...",
                                             initialdir=self.source.get() or str(pathlib.Path.home()))
            if not source:return
            if not (pathlib.Path(source)/"inspection.json").is_file():
                messagebox.showerror("Нет готовых данных","В этой папке не найден inspection.json.\nВыберите каталог LMO_PF2D_BATCH_... с папками point_001 и map.npz.")
                return
            parent=pathlib.Path(self.output_parent.get().strip().strip('"'))
            if not parent.is_dir():
                messagebox.showerror("Ошибка","Укажите существующую папку для новых PDF.");return
            stamp=dt.datetime.now().strftime("%Y%m%d_%H%M%S")
            dest=parent/f"LMO_PF2D_REPORTS_{stamp}"
            idx=1
            while dest.exists():
                dest=parent/f"LMO_PF2D_REPORTS_{stamp}_{idx:02d}"
                idx+=1
            args=["--analyze","--reports-only","--source",str(source),"--output",str(dest)]
            cmd=([sys.executable,*args] if getattr(sys,"frozen",False) else
                 [sys.executable,str(pathlib.Path(__file__).resolve()),*args])
            self._launch_worker(cmd,dest,"Только отчёты из сохранённых 2D-карт; RASX не читаются")

        def _launch_worker(self,cmd,dest,description):
            self.report_dir=dest
            self.progress.set("Строятся отчёты…")
            self.start_button.configure(state="disabled")
            self.scan_button.configure(state="disabled")
            self.rebuild_button.configure(state="disabled")
            self.cancel_button.configure(state="normal")
            self._log("\nНАЧАЛО: "+str(dest))
            self._log("Режим: "+description)
            def task():
                try:
                    flags=getattr(subprocess,"CREATE_NO_WINDOW",0)
                    self.worker=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                        text=True,encoding="utf-8",errors="replace",bufsize=1,creationflags=flags)
                    for line in self.worker.stdout:
                        self.events.put(("log",line.rstrip()))
                    rc=self.worker.wait()
                    self.events.put(("done",rc))
                except Exception:
                    self.events.put(("worker_error",traceback.format_exc()))
            threading.Thread(target=task,daemon=True).start()

        def _cancel(self):
            proc=self.worker
            if not proc or proc.poll() is not None:return
            from tkinter import messagebox
            if not messagebox.askyesno("Остановить","Остановить текущий расчёт? Уже созданные файлы останутся в папке результатов."):return
            self._log("Остановка по запросу пользователя…")
            try:
                if os.name=="nt":
                    subprocess.run(["taskkill","/PID",str(proc.pid),"/T","/F"],stdout=subprocess.DEVNULL,
                                   stderr=subprocess.DEVNULL,creationflags=getattr(subprocess,"CREATE_NO_WINDOW",0))
                else:proc.terminate()
            except Exception as exc:self._log("Ошибка остановки: "+str(exc))

        def _open_report(self):
            from tkinter import messagebox
            path=self.report_dir
            if not path or not path.exists():
                messagebox.showinfo("Нет результатов","Сначала запустите обработку.");return
            overview=path/"LMO011_XY_Overall_Report.pdf"
            target=overview if overview.exists() else path
            try:
                if os.name=="nt":os.startfile(str(target))
                elif sys.platform=="darwin":subprocess.Popen(["open",str(target)])
                else:subprocess.Popen(["xdg-open",str(target)])
            except Exception as exc:messagebox.showerror("Не открыть файл",str(exc))

        def _finish(self,rc):
            self.worker=None
            self.start_button.configure(state="normal")
            self.scan_button.configure(state="normal")
            self.rebuild_button.configure(state="normal")
            self.cancel_button.configure(state="disabled")
            self.progress.set("Готово" if rc==0 else f"Ошибка / код {rc}")
            path=self.report_dir
            self._log(("УСПЕШНО" if rc==0 else "ЗАВЕРШЕНО С ОШИБКОЙ")+" → "+str(path))
            if rc==0 and path:
                from tkinter import messagebox
                messagebox.showinfo("Отчёт создан","Результаты сохранены:\n"+str(path))

        def _poll(self):
            for _ in range(80):
                try:kind,payload=self.events.get_nowait()
                except queue.Empty:break
                if kind=="scan_ok":
                    self.scan_active=False
                    self.scan_button.configure(state="normal")
                    self.start_button.configure(state="normal")
                    self.rebuild_button.configure(state="normal")
                    self.progress.set("Проверка завершена")
                    self._display_scan(payload)
                elif kind=="scan_error":
                    self.scan_active=False
                    self.scan_button.configure(state="normal")
                    self.start_button.configure(state="normal")
                    self.rebuild_button.configure(state="normal")
                    self.progress.set("Ошибка проверки")
                    self.info.set("Ошибка проверки: "+payload)
                    self._log("Ошибка проверки: "+payload)
                elif kind=="log":self._log(payload)
                elif kind=="done":self._finish(payload)
                elif kind=="worker_error":
                    self._log(payload);self._finish(1)
            self._scan_heartbeat()
            self.after(140,self._poll)

        def _close(self):
            if self.worker is not None and self.worker.poll() is None:
                from tkinter import messagebox
                messagebox.showwarning("Расчёт выполняется","Сначала остановите расчёт или дождитесь завершения.")
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
