#!/usr/bin/python3

import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from scipy.interpolate import griddata

import utilities
from readvasp import from_poscar


def plot_band_dos(KPT, eigdata, dosdata, config, axs=None, *args, **kwargs):
    ispin = eigdata.eig.shape[0]
    if axs is None:
        _, axs = plt.subplots(1,
                              len(config["groupinfo"])*ispin + 1,
                              sharey=True)
    for iax in axs[1:-1]:
        iax.sharex(axs[0])
    plot_band(KPT, eigdata, config, axs=axs[:-1], *args, **kwargs)
    plot_dos(dosdata, config, swap=True, ax=axs[-1], *args, **kwargs)
    [i.set_ylabel('') for i in axs[1:]]


def plot_band(KPT, data, config, axs=None, *args, **kwargs):
    ispin = data.eig.shape[0]
    if not "f" in config["plotset"]["fermi"].lower():
        try:
            fermi = utilities.get_energy(data.eig,
                                      data.occ,
                                      info=config["plotset"]["fermi"])
        except:
            print("unknown type ,set fermi to 0eV")
            fermi = data.fermi
        config["plotset"]["fermi"] = "%s" % fermi
        eigall = data.eig - fermi
    else:
        eigall = data.eig
    used_points=[not "!" in i for i in KPT.k_sp_label]
    kpath=KPT.kpath[used_points]
    eigall=eigall[:,used_points]
    k_sp_label=[i for i in KPT.k_sp_label if not "!" in i]
    if axs is None:
        _, axs = plt.subplots(1,
                              len(config["groupinfo"])*ispin,
                              sharey=True, sharex=True)
    if type(axs).__name__ == "AxesSubplot":
        axs = np.array([axs])

    for iaxs, group in zip(axs.reshape(-1, ispin), config["groupinfo"]):
        for i in range(ispin):
            ax = iaxs[i]
            e_range = np.where((eigall[i].max(axis=0) > config["plotset"]["Elim"][0]) & (
                eigall[i].min(axis=0) < config["plotset"]["Elim"][1]))[0]
            eig = eigall[:, :, e_range]
            if config["plotset"]["plot_type"] == 0:
                sc = plot_band_type_0(
                    kpath, eig[i], ax=ax, c='black', *args, **kwargs)
            elif config["plotset"]["plot_type"] >= 1:
                pos = from_poscar.get_poscar(config["fileset"]["posfile"])
                tag, group_info = utilities.group_info_from_input(group)
                project_group = data.set_group(
                    grouptag=tag, symbollist=pos.get_symbollist())
                project_group = project_group[:, :, :, e_range]
                if config["plotset"]["plot_type"] == 1:
                    sc = plot_band_type_1(kpath, eig[i],
                                          project_group=project_group[i],
                                          group_info=group_info,
                                          ax=ax,
                                          int=config["plotset"]["int"],
                                          scale=config["plotset"]["scale"],
                                          *args,
                                          **kwargs)
                    if i == (len(range(ispin))-1):
                        legend_elements = [
                            Line2D([0], [0], color=j, marker='o', label=i) for i, j in group_info]
                        ax.legend(handles=legend_elements, **config["paras"]["band_legend_paras"])
                elif config["plotset"]["plot_type"] == 2:
                    ax = iaxs[i]
                    sc = plot_band_type_2(kpath,
                                          eig[i],
                                          project_group=project_group[i],
                                          group_info=group_info,
                                          ax=ax,
                                          *args,
                                          **kwargs)
                    if i == range(ispin)[-1]:
                        cb = plt.colorbar(sc, ax=ax)
                        cb.set_ticks([0, 1])
                        if len(group_info) == 1:
                            cb.set_ticklabels([r'', r'' + group_info[0][0]])
                        else:
                            cb.set_ticklabels(
                                [r'' + group_info[1][0], r'' + group_info[0][0]])
            else:
                print("unkown type for plotting, available value is 0,1,2")
                sys.exit()
    for ax in axs:
        plot_sp_kline(kpath, k_sp_label, ax=ax, *args, **kwargs)
        ax.set_ylim(config["plotset"]["Elim"])
        ax.set_xlim(kpath.min(), kpath.max())
        if config["plotset"]["fermi"]:
            ax.set_ylabel(r'$E\ -\ E_\mathrm{f}\ \mathrm{(eV)}$')
        else:
            ax.set_ylabel(r'$E\ \mathrm{eV}$')


def plot_dos(data, config, swap=False, ax=None, *args, **kwargs):
    ispin = data.total.shape[0]
    dos = data.total
    type = config["plotset"]["plot_type"]
    if not "f" in config["plotset"]["fermi"].lower():
        if config["plotset"]["fermi"].lower() == "vbm":
            fermi = data.fermi
        elif "band" in config["plotset"]["fermi"].lower():
            print("DOS only supports setting fermi or a specify energy to 0eV")
            fermi = data.fermi
        else:
            fermi = utilities.get_energy(info=config["plotset"]["fermi"])
        if fermi is None:
            print("DOS only supports setting fermi or a specify energy to 0eV")
            fermi = data.fermi
        eig = data.eig - fermi
        print("Fermi is set to %s" % fermi)
    else:
        eig = data.eig
    e_range = np.where((config["plotset"]["Elim"][0] < eig)
                       & (eig < config["plotset"]["Elim"][1]))[0]
    # print(0,e_range[0]-1)
    e_range = [max(0, e_range[0]-1)] + e_range.tolist() + \
        [min(e_range[-1]+1, eig.size - 1)]
    eig = eig[e_range]
    dos = dos[:, e_range]
    if ax is None:
        _, ax = plt.subplots()
    tag = []
    info = []
    if type > 0:
        for group in config["groupinfo"]:
            t_tag, t_info = utilities.group_info_from_input(group)
            tag.extend(t_tag)
            info.extend(t_info)

    for i in range(ispin):
        sign = -np.sign(i - 0.5)
        if type == 0:
            plot_dos_type_0(eig, sign*dos[i], swap=swap, ax=ax, *args, **kwargs)
        elif type > 0:
            if swap:
                ax.fill_betweenx(
                    eig, 0, sign*dos[i], color='gray', alpha=.25, label="Total")
            else:
                ax.fill_between(eig, sign*dos[i], color='gray',
                                alpha=.25, label="Total")
            pos = from_poscar.get_poscar(config["fileset"]["posfile"])
            proj = data.set_group(
                grouptag=tag, symbollist=pos.get_symbollist())
            proj = proj[:, :, e_range]
            plot_dos_type_1(eig,
                            dos[i],
                            sign * proj[i],
                            info,
                            swap=swap,
                            ax=ax,
                            *args,
                            **kwargs)
        else:
            print("unkwon type for dos")
            sys.exit()

    if type > 0:
        xlabel = "PDOS"
        legend_elements = [Line2D([0], [0], color="gray", lw=2, label="Total")] + \
            [Line2D([0], [0], color=j, lw=2, label=i) for i, j in info]
        ax.legend(handles=legend_elements, **config["paras"]["dos_legend_paras"])
        doslim = 1.1 * np.nanmax(proj[:, :, ])
    else:
        xlabel = "DOS"
        doslim = 1.1 * np.nanmax(dos[:, ])
    if swap:
        ax.set_xlabel(xlabel)
        if config["plotset"]["fermi"]:
            ax.set_ylabel(r'$\mathrm{E\ -\ E}_f\ \mathrm{(eV)}$')
        else:
            ax.set_ylabel(r'$\mathrm{E}\ \mathrm{(eV)}$')
        ax.set_ylim(config["plotset"]["Elim"])
        ax.set_xlim([-doslim*i, doslim])
    else:
        ax.set_ylabel(xlabel)
        if config["plotset"]["fermi"]:
            ax.set_xlabel(r'$\mathrm{E\ -\ E}_f\ \mathrm{(eV)}$')
        else:
            ax.set_xlabel(r'$\mathrm{E}\ \mathrm{(eV)}$')
        ax.set_xlim(config["plotset"]["Elim"])
        ax.set_ylim([-doslim*i, doslim])


def plot_sp_kline(kpath, k_sp_label, ax=None, *args, **kwargs):
    if ax is None:
        ax = plt.subplot()
    # "!" is marked for the scf part in HSE/metaGGA calculate and will be exclude in the plot
    kpath_sp = kpath[[bool(i) for i in k_sp_label]]
    ax.set_xticklabels([])
    if (len(kpath_sp) > 0):
        if (len(kpath) > 1):
            ax.vlines(kpath_sp, -1000, 1000, lw=0.5, color='black', *args, **kwargs,alpha=0.2)
        ax.set_xticks(kpath_sp)
        ax.set_xticklabels([i for i in k_sp_label if i])


def plot_band_type_0(kpath, eig, ax=None, *args, **kwargs):
    if ax is None:
        ax = plt.subplot()
    if kpath.shape[0] == 1:
        if 'c' in kwargs.keys():
            kwargs['color'] = kwargs.pop('c')
        ax.hlines(eig, -0.5, 0.5, *args, **kwargs)
    else:
        ax.plot(kpath, eig, *args, **kwargs)


def plot_band_type_1(kpath,
                     eig,
                     project_group,
                     group_info,
                     ax=None,
                     int=100,
                     scale=50,
                     *args,
                     **kwargs):
    if ax is None:
        ax = plt.subplot()
    project_group = np.cumsum(project_group[::-1], axis=0)[::-1]
    plot_band_type_0(kpath, eig, ax=ax, color='gray', *args, **kwargs)
    if kpath.shape[0] == 1:
        kpath = np.array([-0.5, 0.5])
        eig = eig.repeat(2, 0)
        project_group = project_group.repeat(2, 1)
        # print(eig.shape, project_group.shape)
    x = np.linspace(kpath.min(), kpath.max(), int)
    y = griddata(kpath, eig, x)
    for i, igroup in enumerate(project_group):
        label = r'' + group_info[i][0]
        Nband = eig.shape[1]
        Y = y
        X = x.repeat(Nband)
        s = griddata(kpath, igroup, x)
        ax.scatter(X,
                   Y,
                   c=group_info[i][1],
                   s=s * scale,
                   linewidths=None,
                   label=label, zorder=-10)


def plot_band_type_2(kpath,
                     eig,
                     project_group,
                     group_info,
                     ax=None,
                     size=1,
                     *args,
                     **kwargs):
    if ax is None:
        ax = plt.subplot()
    norm = plt.Normalize(0, 1)
    for ieig in range(eig.shape[1]):
        x = kpath.repeat(2)[:-1]
        x[1::2] += (x[2::2] - x[1::2]) / 2
        y = eig[:, ieig].repeat(2)[:-1]
        y[1::2] += (y[2::2] - y[1::2]) / 2
        c = project_group[0, :, ieig].repeat(2)[:-1]
        c[1::2] += (c[2::2] - c[1::2]) / 2
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap="rainbow", norm=norm)
        lc.set_array(c)
        line = ax.add_collection(lc)
    return line


def plot_dos_astype(type,
                    eig,
                    dos,
                    project_group,
                    group_info,
                    swap=False,
                    ax=None,
                    *args,
                    **kwargs):
    if ax is None:
        ax = plt.subplot()
    if type == 0:
        plot_dos_type_0(eig, dos, swap=swap, ax=ax, *args, **kwargs)
    elif type > 0:
        plot_dos_type_1(eig,
                        dos,
                        project_group,
                        group_info,
                        swap=swap,
                        ax=ax,
                        *args,
                        **kwargs)
    else:
        print("unkwon type for dos")
        sys.exit()


def plot_dos_type_0(eig, dos, swap=False, ax=None, *args, **kwargs):
    if ax is None:
        ax = plt.subplot()
    if swap:
        y, x = eig, dos
    else:
        x, y = eig, dos
    ax.plot(x, y, c='black', zorder=-10, *args, **kwargs)


def plot_dos_type_1(eig,
                    dos,
                    project_group,
                    group_info,
                    swap=False,
                    ax=None,
                    *args,
                    **kwargs):
    if ax is None:
        ax = plt.subplot()
    for i, igroup in enumerate(project_group):
        label = r'' + group_info[i][0]
        if swap:
            ax.plot(igroup, eig, c=group_info[i]
                    [1], zorder=-10, *args, **kwargs)
            ax.fill_betweenx(eig,
                             0,
                             igroup,
                             color=group_info[i][1],
                             alpha=.25,
                             label=label,
                             *args,
                             **kwargs)
        else:
            ax.plot(eig, igroup, c=group_info[i]
                    [1], zorder=-10, *args, **kwargs)
            ax.fill_between(eig,
                            igroup,
                            color=group_info[i][1],
                            alpha=.25,
                            label=label,
                            *args,
                            **kwargs)


def line_average_of_chg(chg, abc_axe, ax=None):
    """
    average alone a axis
    _axe = 
    'a','b','c' alone the lattice constant
    //'x','y','z' alone the x,y,z direction
    """
    if ax is None:
        ax = plt.subplot()
    _axe = ['c', 'b', 'a'].index(abc_axe)
    # _axe=int(sys.argv[1])-1
    axe = [0, 1, 2]
    axe.remove(_axe)
    axe = tuple(axe)
    NG = chg.shape
    ax.plot(range(NG[_axe]) / NG[_axe] *
            np.linalg.norm(chg.cell[(2, 1, 0)[_axe], :]),
            np.average(chg.chg, axis=axe),
            color='black', zorder=-10)
    ax.set_xlim([
        np.min(
            range(NG[_axe]) / NG[_axe] *
            np.linalg.norm(chg.cell[(2, 1, 0)[_axe], :])),
        np.max(
            range(NG[_axe]) / NG[_axe] *
            np.linalg.norm(chg.cell[(2, 1, 0)[_axe], :]))
    ])
