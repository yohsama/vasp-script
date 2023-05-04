import re

import numpy as np

kb = 1.380649*10**-23


def get_fermi(EIG, OCC, OCC_E=0.001):
    """
    this method used to get the fermi level through the occupid number, the occupid threshold are defined by OCC_E, default is 0.001.
    """
    return np.max(EIG[OCC > OCC_E])


def get_vbm(EIG, OCC, OCC_E=0.001):
    """
    this method used to get the fermi level through the occupid number, the occupid threshold are defined by OCC_E, default is 0.001.
    """
    return np.max(EIG[OCC > OCC_E])


def get_cbm(EIG, OCC, OCC_E=0.001):
    """
    this method used to get the fermi level through the occupid number, the occupid threshold are defined by OCC_E, default is 0.001.
    """
    return np.min(EIG[OCC < OCC_E])


def orbital_to_index(L_orbit, orbital):
    if L_orbit == 10:
        o2i = {
            's': [0],
            'py': [1],
            'pz': [1],
            'px': [1],
            'dxy': [2],
            'dyz': [2],
            'dz2': [2],
            'dx2-y2': [2],
            'dxz': [2],
            'p': [1],
            'd': [2]
        }
    elif L_orbit > 10:
        o2i = {
            's': [0],
            'py': [1],
            'pz': [2],
            'px': [3],
            'dxy': [4],
            'dyz': [5],
            'dz2': [6],
            'dx2-y2': [7],
            'dxz': [8],
            'p': [1, 2, 3],
            'd': [4, 5, 6, 7, 8]
        }
    return np.array(o2i[orbital])


def set_group(L_orbit, project, grouptag, symbollist, norm=True):
    """
    this method used to  group the project data by given element|symbol|atomic number and/or orbital.
    projectdata could be object:dos,pro from from_doscar.get_doscar() or from_procar.get_procar()
    grouptag should be writen as:
    symbollist shold be writen as string list.
    """
    symbollist = np.array(symbollist)
    for i,igrouptag in enumerate(grouptag):
        tot = 0
        for igroup in igrouptag.split(','):
            whichorbit = ...
            if "_" in igroup:
                whichorbit = orbital_to_index(L_orbit, igroup.split('_')[1])
            whichatom = igroup.split('_')[0]
            if re.search('[A-Z]',
                         whichatom):  # element symbolliste e.g Ag,H,C ...
                whichatom = (symbollist == whichatom).nonzero()[0]
            else:
                if re.search('-',
                             whichatom):  # element number e.g 1-10,11-23 ...
                    start = int(whichatom.split('-')[0]) - 1
                    end = int(whichatom.split('-')[1])
                    whichatom = np.arange(start, end)
                else:
                    whichatom = np.array(
                        [whichatom],
                        dtype=int) - 1  # element number e.g 1,2,3...
            if len(project.shape) > 4:
                tot += project[:, :, :, whichatom][:, :, :, :, whichorbit].sum(
                    (-1, -2))
            else:
                tot += project[:, :, whichatom][:, :, :, whichorbit].sum(
                    (-1, -2))
        if i == 0:
            pro_group = tot[:, None]
        else:
            pro_group = np.append(pro_group, tot[:, None], axis=1)
            
    return pro_group


def group_info_from_input(groupinput):
    group_tag = []
    group_info = []
    for tmp in groupinput:
        type1 = tmp.split()
        group_tag.append(type1[0])
        group_info.append(type1[1:])
    return group_tag, group_info


def get_energy(eig=None, occ=None, info=None):
    fermi = None
    try:
        if info.lower() == "vbm":
            fermi = get_vbm(eig, occ)
        elif info.lower() == "cbm":
            fermi = get_cbm(eig, occ)
        elif "band" in info.lower():
            ispin = int(info.split("_")[1])
            N = int(info.split("_")[2])
            fermi = np.max(eig[ispin,:, N-1])
            print("band fermi",fermi)
    except Exception as e:
        print("get_ferm failed",e)
    if fermi is None:
        try:
            fermi = float(info)
        except Exception as e:
            print("unknown fermi input")
    return fermi

def color_line(X,Y,C,ax=None,cmap="rainbow"):
    from matplotlib.collections import LineCollection
    if ax is None:
        ax=plt.subplot()
    norm = plt.Normalize(0, 1)
    x = X.repeat(2)[:-1]
    x[1::2] += (x[2::2] - x[1::2]) / 2
    y = Y.repeat(2)[:-1]
    y[1::2] += (y[2::2] - y[1::2]) / 2
    c = C.repeat(2)[1:]
    #c[1::2] += (c[2::2] - c[1::2]) / 2
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=cmap,norm=norm)
    lc.set_array(c)
    line = ax.add_collection(lc)
    return line


def color_multiline(X,Y,C,ax=None,cmap="rainbow"):
    from matplotlib.collections import LineCollection
    import matplotlib.pyplot as plt
    if ax is None:
        ax=plt.subplot()
    norm = plt.Normalize(0, 1)
    x = X.repeat(2)[:-1]
    x[1::2] += (x[2::2] - x[1::2]) / 2
    y = Y.repeat(2,axis=0)[:-1]    
    y[1::2] += (y[2::2] - y[1::2]) / 2
    c = C.repeat(2,axis=0)[1:-1]
    #c[1::2] += (c[2::2] - c[1::2]) / 2
    segments=np.empty((0,2,2))
    for iy in y.T:
        points = np.array([x, iy]).T.reshape(-1, 1, 2)
        st=np.concatenate([points[:-1], points[1:]], axis=1)
        segments=np.concatenate((segments,st),axis=0)
    lc = LineCollection(segments, cmap=cmap,norm=norm)
    lc.set_array(c.T.flatten())
    line = ax.add_collection(lc)
    return line

def get_nearest_distence(dx):
    return np.abs((dx+0.5)%1)-0.5

kb=1.380649*10**-23
