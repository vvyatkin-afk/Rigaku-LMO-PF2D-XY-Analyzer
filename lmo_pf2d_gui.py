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
import traceback
import zipfile
import xml.etree.ElementTree as ET
from collections import defaultdict

VERSION="1.0.0"
TITLE="LMO PF2D Analyzer  |  HyPix-3000(H)"
RASX_PATTERN=re.compile(
    r"^(?P<stem>.+)_image_(?P<seq>\d+)_(?P<chi>-?\d+(?:\.\d+)?)_(?P<num>\d+)\.rasx$",re.I)

def bundled(rel: str) -> pathlib.Path:
    root=pathlib.Path(getattr(sys, "_MEIPASS", pathlib.Path(__file__).resolve().parent))
    return root/rel

def scan_groups(folder: str | pathlib.Path) -> list[dict]:
    """Quick filename-level inspection; physics-level checks happen in the worker."""
    root=pathlib.Path(folder)
    if not root.is_dir():
        raise FileNotFoundError(f"Folder does not exist or is inaccessible: {root}")
    groups: dict[str,list]=defaultdict(list)
    unmatched=0
    for entry in root.rglob("*.rasx"):
        m=RASX_PATTERN.match(entry.name)
        if m:
            key=f"{entry.parent.relative_to(root)}__{m['stem']}"
            groups[key].append((entry, int(m['seq']), float(m['chi']), int(m['num'])))
        else:
            unmatched+=1
    rows=[]
    expected={(seq,chi,n) for seq,chi in ((1,38.0),(2,42.0)) for n in range(1,181)}
    for key,items in sorted(groups.items(),key=lambda x:x[0].lower()):
        indexes=[(s,c,n) for _,s,c,n in items]
        actual=set(indexes)
        counts=len(items)
        ok=counts==360 and len(actual)==360 and actual==expected
        xy=""
        if items:
            try:
                with zipfile.ZipFile(items[0][0]) as z:
                    header=ET.fromstring(z.read("Data0/MesurementConditions0.xml"))
                    ax={v.get('Name'):v.get('Position') for v in header.findall("Axes/Axis")}
                x=float(ax.get('X','nan'))
                y=float(ax.get('Y','nan'))
                xy=f"X={x:.2f}, Y={y:.2f} mm"
            except Exception as exc:
                xy=f"Metadata read failed: {type(exc).__name__}"
        rows.append(dict(key=key,frames=counts,expected=360,
                         ready=ok,chi38=sum(c==38 for _,_,c,_ in items),
                         chi42=sum(c==42 for _,_,c,_ in items),
                         x_y=xy))
    if not rows and unmatched:
        raise ValueError(f"{unmatched} RASX files found, but filenames are not in the supported format.")
    return rows

def run_script_mode(script: str, args: list[str]):
    allowed={"run_lmo_pf.py","reconstruct_poles.py","fit_poles.py","render_results.py"}
    if script not in allowed:
        raise ValueError("Unsupported worker script: "+script)
    target=bundled("analysis/validated_single_point_scripts/"+script)
    if not target.is_file():
        raise FileNotFoundError("Bundled worker not found: "+str(target))
    sys.argv=[str(target),*args]
    runpy.run_path(str(target),run_name="__main__")

def run_analyze_mode(args: list[str]):
    target=bundled("analysis/analyze_many_xy.py")
    if not target.is_file():
        raise FileNotFoundError("Bundled analyzer not found: "+str(target))
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
            self.info=tk.StringVar(value="Выберите папку с файлами .rasx и нажмите «Проверить данные».")
            self.progress=tk.StringVar(value="Ожидание")
            self.rows=[]
            self.worker=None
            self.events=queue.Queue()
            self.report_dir=None
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
            self.cancel_button=ttk.Button(action,text="Остановить",command=self._cancel,state="disabled")
            self.cancel_button.pack(side="left",padx=(0,9))
            ttk.Button(action,text="Открыть результаты",command=self._open_report).pack(side="left")
            ttk.Label(action,textvariable=self.progress,style="Text.TLabel").pack(side="right",padx=8)
            ttk.Label(root,text="Журнал обработки",style="Text.TLabel").pack(anchor="w",pady=(7,5))
            self.log=scrolledtext.ScrolledText(root,height=8,font=("Consolas",9),state="disabled",
                 background="#101c2c",foreground="#d7ebfb",insertbackground="#ffffff",wrap="word")
            self.log.pack(fill="both",expand=False)
            self._log("Программа не управляет дифрактометром. Поддерживается геометрия примера Ивана: χ=38°/42°, φ шаг 2°.")
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
            if self.worker is not None:return
            source=self.source.get().strip().strip('"')
            if not source:
                messagebox.showwarning("Папка не задана","Укажите папку с RASX.")
                return
            self.info.set("Поиск и проверка серий RASX…")
            self.scan_button.configure(state="disabled")
            def scan_task():
                try:self.events.put(("scan_ok",scan_groups(source)))
                except Exception as exc:self.events.put(("scan_error",str(exc)))
            threading.Thread(target=scan_task,daemon=True).start()

        def _display_scan(self,items):
            self.rows=items
            for row_id in self.tree.get_children():self.tree.delete(row_id)
            for i,r in enumerate(items):
                self.tree.insert("",'end',iid=str(i),text=r["key"],values=(
                    r["frames"],r["chi38"],r["chi42"],r["x_y"],"Готово" if r["ready"] else "Неполно"))
            if items:self.tree.selection_set("0")
            ready=sum(r["ready"] for r in items)
            self.info.set(f"Найдено серий: {len(items)}. Полных 360-кадровых: {ready}. "
                          "Для контроля геометрии каждого кадра требуется запуск анализа.")
            self._log(self.info.get())

        def _start(self):
            if self.worker is not None:return
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
                    messagebox.showwarning("Неполная серия","Для выбранной точки нет полного набора из 360 кадров.");return
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
            cmd=[sys.executable,"--analyze",*args,"--output",str(dest)] if getattr(sys,"frozen",False) else [
                sys.executable,str(pathlib.Path(__file__).resolve()),"--analyze",*args,"--output",str(dest)]
            self.report_dir=dest
            self.progress.set("Работает…")
            self.start_button.configure(state="disabled")
            self.scan_button.configure(state="disabled")
            self.cancel_button.configure(state="normal")
            self._log("\nНАЧАЛО: "+str(dest))
            self._log("Режим: "+("одна точка" if kind=="single" else "весь пакет"))
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
                    self.scan_button.configure(state="normal")
                    self._display_scan(payload)
                elif kind=="scan_error":
                    self.scan_button.configure(state="normal")
                    self.info.set("Ошибка проверки: "+payload)
                    self._log(payload)
                elif kind=="log":self._log(payload)
                elif kind=="done":self._finish(payload)
                elif kind=="worker_error":
                    self._log(payload);self._finish(1)
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
