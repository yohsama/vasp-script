#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
from readvasp import *
import plotting2
import yaml
import os


def from_input(config):
    obj = ["band", "DOS", "band & DOS", "CHARGE"]
    str = ''.join(['%s : [%d]\n' % (j, i) for i, j in enumerate(obj)])
    tmp = int(input(str))
    config["plotset"]["object"] = obj[tmp]
    if tmp == 0:
        type = int(input(
            'Type:\n'
            'Normal band    - 0\n'
            'Fat band       - 1\n'
            'Colormap band  - 2\n'
        ))
    if tmp == 1:
        type = int(input(
            'Type:\n'
            'Normal DOS     - 0\n'
            'Projected DOS  - 1\n'
        ))
    if tmp == 2:
        type = int(input(
            'Type:\n'
            'Normal      - 0\n'
            'Projected   - 1\n'
        ))
    config["plotset"]["plot_type"] = type
    config["groupinfo"] = [None]
    if type > 0:
        config["Nbandpicture"] = int(input("How many bands?\n"))
        config["groupinfo"] = []
        for N in range(config["Nbandpicture"]):
            print("next group\n")
            groups = []
            tmp = input(
                '[Element or atomic NO.] [label] [color/colormap] \n'
                'specified orbital with "_" ,split with "," \n'
                'e.g : Cu,Ag_py,1-3,6_py,8-10_d CuAnTag red \n'
                'empty for end \n'
            )
            groups.append(tmp)
            while tmp:
                tmp = input()
                if tmp:
                    groups.append(tmp)
            config["groupinfo"].append(groups)
    tmp = input('Energy range, e.g -2 5  default is [-2 5]\n') or "-2 5"
    config["plotset"]['Elim'] = [float(i) for i in tmp.split()]
    config["plotset"]['fermi'] = input("Set 0eV to fermi? (VBM,CBM,band No. or any Energy)\n"
                                       "F or 0 will change nothing.\n"
                                       "e.g VBM/CBM/band_18/-2.5/F/0\n"
                                       "default is VBM") or "vbm"
    return config


def ini_config():
    method = {"method": "matplotlib"}
    plotset = {"plot_type": 0, "object": "band",
                "scale": 50,
               "Elim": [-2, 5], "fermi": "vbm", "int": 100}
    figset = {"figsize": (7, 6), "dpi": 150}
    fontset = {"fontsize": 16, "font": "Arial"}
    fileset = {"eigfile": ["EIGENVAL"], "dosfile": "DOSCAR", "prosfile": [
        "PROCAR"], "kptfile": ["KPOINTS"], "posfile": "POSCAR"}
    config = {"method": method,
              "plotset": plotset,
              "fileset": fileset,
              "figset": figset,
              "fontset": fontset}
    return config


def read_config(file='cvasp.yml'):
    with open(file) as f:
        data = yaml.load(f, Loader=yaml.FullLoader)
    # print(data)
    return data


def cvasp(config, axs=None):
    print(config["plotset"]["object"].lower())
    if config["plotset"]["object"].lower() == "band":
        KPT = from_kpoints.get_kpoints(config["fileset"]["kptfile"])
        if config["plotset"]["plot_type"] > 0:
            eigdata = from_procar.get_procar(config["fileset"]["prosfile"])
        else:
            eigdata = from_eigenval.get_eigenvalue(config["fileset"]["eigfile"])
        plotting2.plot_band(KPT, eigdata, config, axs=axs)

    if config["plotset"]["object"].lower() == "dos":
        dosdata = from_doscar.get_doscar(config["fileset"]["dosfile"])
        plotting2.plot_dos(dosdata, config, axs=axs)

    if "band" in config["plotset"]["object"].lower() and \
       "dos" in config["plotset"]["object"].lower():
        KPT = from_kpoints.get_kpoints(config["fileset"]["kptfile"])
        if config["plotset"]["plot_type"] > 0:
            eigdata = from_procar.get_procar(config["fileset"]["prosfile"])
        else:
            eigdata = from_eigenval.get_eigenvalue(config["fileset"]["eigfile"])
        dosdata = from_doscar.get_doscar(config["fileset"]["dosfile"])
        plotting2.plot_band_dos(KPT, eigdata, dosdata, config, axs=axs)



if __name__ == "__main__":
    config = ini_config()
    try:
        data = read_config()
        config.update(data)
    except Exception as e:
        config = from_input(config)
        with open('cvasp.yml', "w", encoding='utf-8') as f:
            yaml.dump(config, f)
    plt.rc('font', family=config["fontset"]
           ["font"], size=config["fontset"]["fontsize"])
    plt.rcParams['font.serif'] = [config["fontset"]
                                  ["font"]] + plt.rcParams['font.serif']
    print(config["plotset"]["object"].lower())
    ispin=1
    if "band" in config["plotset"]["object"].lower():
        if config["plotset"]["plot_type"] > 0:
            eigdata = from_procar.get_procar(config["fileset"]["prosfile"])
        else:
            eigdata = from_eigenval.get_eigenvalue(config["fileset"]["eigfile"])
        ispin = eigdata.N_spin
    if config["plotset"]["object"].lower() == "band":
        fig, axs = plt.subplots(1, len(config["groupinfo"])*ispin,
                                figsize=config["figset"]["figsize"], sharey=True)
    elif config["plotset"]["object"].lower() == "dos":
        fig, ax = plt.subplots(1, 1,
                               figsize=config["figset"]["figsize"])
    elif "band" in config["plotset"]["object"].lower() and \
        "dos" in config["plotset"]["object"].lower():
        fig, axs = plt.subplots(1, len(config["groupinfo"])*ispin + 1,
                                figsize=config["figset"]["figsize"], sharey=True)
    cvasp(config, axs=axs)
    plt.show()
    fig.savefig('band_dos.png', dpi=300)