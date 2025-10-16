from __future__ import annotations


from pymatgen.phonon.plotter import freq_units
from pymatgen.util.plotting import add_fig_kwargs, get_ax_fig, pretty_plot
from pymatgen.core import Lattice, Structure
from pymatgen.io.phonopy import get_structure_from_dict
from pymatgen.phonon.plotter import PhononBSPlotter
import numpy as np
import scipy.constants as const
from monty.dev import requires
from monty.json import MSONable
from phonopy.phonon.dos import TotalDos
from pymatgen.core import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.units import amu_to_kg
from pymatgen.phonon.bandstructure import (
    PhononBandStructure,
    PhononBandStructureSymmLine,
)
from pymatgen.phonon.dos import PhononDos
from monty.serialization import loadfn
from typing import TYPE_CHECKING, NamedTuple
import matplotlib.pyplot as plt

class FreqUnits(NamedTuple):
    """Conversion factor from THz to the required units and the label."""

    factor: float
    label: str


def freq_units(units: Literal["thz", "ev", "mev", "ha", "cm-1", "cm^-1"]) -> FreqUnits:
    """
    Args:
        units: str, accepted values: thz, ev, mev, ha, cm-1, cm^-1.

    Returns:
        Conversion factor from THz to the required units and the label in the form of a namedtuple.
    """
    dct = {
        "thz": FreqUnits(1, "THz"),
        "ev": FreqUnits(const.value("hertz-electron volt relationship") * const.tera, "eV"),
        "mev": FreqUnits(
            const.value("hertz-electron volt relationship") * const.tera / const.milli,
            "meV",
        ),
        "ha": FreqUnits(const.value("hertz-hartree relationship") * const.tera, "Ha"),
        "cm-1": FreqUnits(
            const.value("hertz-inverse meter relationship") * const.tera * const.centi,
            "cm^{-1}",
        ),
        "cm^-1": FreqUnits(
            const.value("hertz-inverse meter relationship") * const.tera * const.centi,
            "cm^{-1}",
        ),
    }
    try:
        return dct[units.lower().strip()]
    except KeyError:
        raise KeyError(f"Value for units `{units}` unknown\nPossible values are:\n {list(dct)}")

def get_velocity(velocity_path, structure=None, structure_path=None) -> Velocity:
    """
    Get velocity object from mehs.yaml file including group velocities, as obtained from phonopy (Frequencies in THz!).
    The order is structure > structure path > structure from velocities dict.
    Newer versions of phonopy include the structure in the yaml file,
    the structure/structure_path is kept for compatibility.

    Args:
        velocity_path: Path to mesh.yaml file including velocities (frequencies have to be in THz!)
        structure: pymatgen Structure object
        structure_path: path to structure in a file (e.g., POSCAR)

    Returns: Velocity object

    """

    velocity_dict = loadfn(velocity_path)

    if structure_path and structure is None:
        structure = Structure.from_file(structure_path)
    else:
        try:
            structure = get_structure_from_dict(velocity_dict)
        except ValueError:
            raise ValueError("\nPlease provide a structure.\n")

    qpts, multiplicities, frequencies, velocity = ([] for _ in range(4))
    phonopy_labels_dict = {}

    for p in velocity_dict["phonon"]:
        q = p["q-position"]
        qpts.append(q)
        if "multiplicity" in p:
            m = p["multiplicity"]
        else:
            m = 1
        multiplicities.append(m)
        bands, velocityband = [], []
        for b in p["band"]:
            bands.append(b["frequency"])
            if "group_velocity" in b:
                # adding the vector for now
                velocityband.append(np.sqrt(b["group_velocity"][0]**2+b["group_velocity"][1]**2+b["group_velocity"][2]**2))
        frequencies.append(bands)
        velocity.append(velocityband)
        if "label" in p:
            phonopy_labels_dict[p["label"]] = p["q-position"]

    qpts_np = np.array(qpts)
    multiplicities_np = np.array(multiplicities)
    # transpose to match the convention in PhononBandStructure
    frequencies_np = np.transpose(frequencies)
    velocity_np = np.transpose(velocity)

    return Velocity(
        velocities=velocity_np,
        qpoints=qpts_np,
        multiplicities=multiplicities_np,
        frequencies=frequencies_np,
        structure=structure,
    )



class Velocity(MSONable):
    """
    Class for Group velocities on a regular grid.
    """

    def __init__(
        self,
        qpoints,
        velocities,
        frequencies,
        multiplicities=None,
        structure=None,
        lattice=None,
    ):
        """
        Args:
            qpoints: list of qpoints as numpy arrays, in frac_coords of the given lattice by default
            velocities: list of absolute group velocities parameters as numpy arrays, shape: (3*len(structure), len(qpoints))
            frequencies: list of phonon frequencies in THz as a numpy array with shape (3*len(structure), len(qpoints))
            multiplicities: list of multiplicities
            structure: The crystal structure (as a pymatgen Structure object) associated with the velocities.
            lattice: The reciprocal lattice as a pymatgen Lattice object. Pymatgen uses the physics convention of
                     reciprocal lattice vectors WITH a 2*pi coefficient
        """
        self.qpoints = qpoints
        self.velocities = velocities
        self.frequencies = frequencies
        self.multiplicities = multiplicities
        self.lattice = lattice
        self.structure = structure

    @property  # type: ignore
    def tdos(self):
        """
        The total DOS (re)constructed from the mesh.yaml file
        """

        # Here, we will reuse phonopy classes
        class TempMesh:
            """
            Temporary Class
            """

        a = TempMesh()
        a.frequencies = np.transpose(self.frequencies)
        a.weights = self.multiplicities

        b = TotalDos(a)
        b.run()

        return b

    @property
    def phdos(self):
        """
        Returns: PhononDos object
        """
        return PhononDos(self.tdos.frequency_points, self.tdos.dos)



class VelocityPlotter:
    """
    Class to plot Velocity Object
    """

    def __init__(self, velocity):
        """
        Class to plot information from Velocity Object
        Args:
            velocity: Velocity Object
        """
        self._velocity = velocity

    def get_plot(self, marker="o", markersize=None, units: Literal["thz", "ev", "mev", "ha", "cm-1", "cm^-1"] = "thz",
    ):
        """
        will produce a plot
        Args:
            marker: marker for the depiction
            markersize: size of the marker
            units: unit for the plots, accepted units: thz, ev, mev, ha, cm-1, cm^-1

        Returns: plot
        """
        u = freq_units(units)

        xs = self._velocity.frequencies.flatten() * u.factor
        ys = self._velocity.velocities.flatten() /10
        #print(xs)
        #print(ys)
        ax = pretty_plot(12, 8)
        ax.set_xlim(-0.5, 12)
        ax.set_ylim(-0.5, 3.5)
        ax.set_xlabel(rf"$\mathrm{{Frequency\ ({u.label})}}$")
        ax.set_ylabel(r"$\mathrm{Group\ Velocity}$ (Km/s)")
        ax.axhline(0, linewidth=1, color="black")

        n = len(ys) - 1
        for i, (x, y) in enumerate(zip(xs, ys)):
            color = (1.0 / n * i, 0, 1.0 / n * (n - i))

            if markersize:
                ax.plot(x, y, marker, color=color, markersize='0.1')
            else:
                ax.plot(x, y, marker, color=color)

        plt.tight_layout()

        return ax

    def show(self, units="thz"):
        """
        will show the plot
        Args:
            units: units for the plot, accepted units: thz, ev, mev, ha, cm-1, cm^-1

        Returns: plot
        """
        self.get_plot(units=units)
        plt.show()

    def save_plot(self, filename, img_format="png", units="thz"):
        """
        Will save the plot to a file
        Args:
            filename: name of the filename
            img_format: format of the saved plot
            units: accepted units: thz, ev, mev, ha, cm-1, cm^-1

        Returns:
        """
        self.get_plot(units=units)
        plt.savefig(filename, format=img_format)
        plt.close()

class VelocityPhononBandStructure(PhononBandStructure):
    """
    This is the most generic phonon band structure data possible
    it's defined by a list of qpoints + frequencies for each of them.
    Additional information may be given for frequencies at Gamma, where
    non-analytical contribution may be taken into account.
    """

    def __init__(
        self,
        qpoints,
        frequencies,
        velocities,
        lattice,
        eigendisplacements=None,
        labels_dict=None,
        coords_are_cartesian=False,
        structure=None,
    ):
        """
        Args:
            qpoints: list of qpoint as numpy arrays, in frac_coords of the
                given lattice by default
            frequencies: list of phonon frequencies in THz as a numpy array with shape
                (3*len(structure), len(qpoints)). The First index of the array
                refers to the band and the second to the index of the qpoint.
            velocities: list of group velocity parameters with the same structure
                frequencies.
            lattice: The reciprocal lattice as a pymatgen Lattice object.
                Pymatgen uses the physics convention of reciprocal lattice vectors
                WITH a 2*pi coefficient.
            eigendisplacements: the phonon eigendisplacements associated to the
                frequencies in Cartesian coordinates. A numpy array of complex
                numbers with shape (3*len(structure), len(qpoints), len(structure), 3).
                The first index of the array refers to the band, the second to the index
                of the qpoint, the third to the atom in the structure and the fourth
                to the Cartesian coordinates.
            labels_dict: (dict) of {} this links a qpoint (in frac coords or
                Cartesian coordinates depending on the coords) to a label.
            coords_are_cartesian: Whether the qpoint coordinates are Cartesian.
            structure: The crystal structure (as a pymatgen Structure object)
                associated with the band structure. This is needed if we
                provide projections to the band structure
        """
        PhononBandStructure.__init__(
            self,
            qpoints,
            frequencies,
            lattice,
            nac_frequencies=None,
            eigendisplacements=eigendisplacements,
            nac_eigendisplacements=None,
            labels_dict=labels_dict,
            coords_are_cartesian=coords_are_cartesian,
            structure=structure,
        )
        self.velocites = velocities

    def as_dict(self):
        """
        Returns:
            MSONable (dict)
        """
        d = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "lattice_rec": self.lattice_rec.as_dict(),
            "qpoints": [],
        }
        # qpoints are not Kpoint objects dicts but are frac coords. This makes
        # the dict smaller and avoids the repetition of the lattice
        for q in self.qpoints:
            d["qpoints"].append(q.as_dict()["fcoords"])
        d["bands"] = self.bands.tolist()
        d["labels_dict"] = {}
        for kpoint_letter, kpoint_object in self.labels_dict.items():
            d["labels_dict"][kpoint_letter] = kpoint_object.as_dict()["fcoords"]
        # split the eigendisplacements to real and imaginary part for serialization
        d["eigendisplacements"] = dict(
            real=np.real(self.eigendisplacements).tolist(), imag=np.imag(self.eigendisplacements).tolist()
        )
        d["velocities"] = self.velocites.tolist()
        if self.structure:
            d["structure"] = self.structure.as_dict()

        return d

    @classmethod
    def from_dict(cls, d):
        """
        Args:
            d (dict): Dict representation
        Returns:
            VelocityPhononBandStructure: Phonon band structure with Velocity parameters.
        """
        lattice_rec = Lattice(d["lattice_rec"]["matrix"])
        eigendisplacements = np.array(d["eigendisplacements"]["real"]) + np.array(d["eigendisplacements"]["imag"]) * 1j
        structure = Structure.from_dict(d["structure"]) if "structure" in d else None
        return cls(
            qpoints=d["qpoints"],
            frequencies=np.array(d["bands"]),
            velocities=np.array(d["velocities"]),
            lattice=lattice_rec,
            eigendisplacements=eigendisplacements,
            labels_dict=d["labels_dict"],
            structure=structure,
        )


class VelocityPhononBandStructureSymmLine(VelocityPhononBandStructure, PhononBandStructureSymmLine):
    """
    This object stores a VelocityPhononBandStructureSymmLine together with group velocity
    for every frequency.
    """

    def __init__(
        self,
        qpoints,
        frequencies,
        velocities,
        lattice,
        eigendisplacements=None,
        labels_dict=None,
        coords_are_cartesian=False,
        structure=None,
    ):
        """
        Args:
            qpoints: list of qpoints as numpy arrays, in frac_coords of the
                given lattice by default
            frequencies: list of phonon frequencies in eV as a numpy array with shape
                (3*len(structure), len(qpoints))
            velocities: list of absolute velocities as a numpy array with the
                shape (3*len(structure), len(qpoints))
            lattice: The reciprocal lattice as a pymatgen Lattice object.
                Pymatgen uses the physics convention of reciprocal lattice vectors
                WITH a 2*pi coefficient
            eigendisplacements: the phonon eigendisplacements associated to the
                frequencies in Cartesian coordinates. A numpy array of complex
                numbers with shape (3*len(structure), len(qpoints), len(structure), 3).
                The first index of the array refers to the band, the second to the index
                of the qpoint, the third to the atom in the structure and the fourth
                to the Cartesian coordinates.
            labels_dict: (dict) of {} this links a qpoint (in frac coords or
                Cartesian coordinates depending on the coords) to a label.
            coords_are_cartesian: Whether the qpoint coordinates are cartesian.
            structure: The crystal structure (as a pymatgen Structure object)
                associated with the band structure. This is needed if we
                provide projections to the band structure
        """
        VelocityPhononBandStructure.__init__(
            self,
            qpoints=qpoints,
            frequencies=frequencies,
            velocities=velocities,
            lattice=lattice,
            eigendisplacements=eigendisplacements,
            labels_dict=labels_dict,
            coords_are_cartesian=coords_are_cartesian,
            structure=structure,
        )

        PhononBandStructureSymmLine._reuse_init(
            self, eigendisplacements=eigendisplacements, frequencies=frequencies, has_nac=False, qpoints=qpoints
        )

    @classmethod
    def from_dict(cls, d):
        """
        Args:
            d: Dict representation
        Returns: VelocityPhononBandStructureSummLine
        """
        lattice_rec = Lattice(d["lattice_rec"]["matrix"])
        eigendisplacements = np.array(d["eigendisplacements"]["real"]) + np.array(d["eigendisplacements"]["imag"]) * 1j
        structure = Structure.from_dict(d["structure"]) if "structure" in d else None
        return cls(
            qpoints=d["qpoints"],
            frequencies=np.array(d["bands"]),
            velocities=np.array(d["velocities"]),
            lattice=lattice_rec,
            eigendisplacements=eigendisplacements,
            labels_dict=d["labels_dict"],
            structure=structure,
        )


def get_velocity_ph_bs_symm_line(velocity_path, structure=None, structure_path=None, labels_dict=None):
    r"""
    Creates a pymatgen VelocityPhononBandStructure from a band.yaml file.
    The labels will be extracted from the dictionary, if present.
    If the 'eigenvector' key is found the eigendisplacements will be
    calculated according to the formula:
    \\exp(2*pi*i*(frac_coords \\dot q) / sqrt(mass) * v
     and added to the object.
    Args:
        velocity_path: path to the band.yaml file
        structure: pymaten Structure object
        structure_path: path to a structure file (e.g., POSCAR)
        labels_dict: dict that links a qpoint in frac coords to a label.

    """
    #print("test")
    dict_here=loadfn(velocity_path)
    #print(dict_here)
    return get_velocity_ph_bs_symm_line_from_dict(dict_here, structure, structure_path, labels_dict)



def get_velocity_ph_bs_symm_line_from_dict(
    velocity_dict, structure=None, structure_path=None, labels_dict=None
) -> VelocityPhononBandStructureSymmLine:
    r"""
    Creates a pymatgen VecolityPhononBandStructure object from the dictionary
    extracted by the band.yaml file including group velocties produced by phonopy. The labels
    will be extracted from the dictionary, if present. If the 'eigenvector'
    key is found the eigendisplacements will be calculated according to the
    formula::
        exp(2*pi*i*(frac_coords \\dot q) / sqrt(mass) * v
    and added to the object.
    Args:
        velocity_dict (dict): the dictionary extracted from the band.yaml file
        structure (Structure): pymatgen structure object
        structure_path: path to structure file
        labels_dict (dict): dict that links a qpoint in frac coords to a label.
            Its value will replace the data contained in the band.yaml.
    """
    #print(velocity_dict)
    if structure_path and structure is None:
        structure = Structure.from_file(structure_path)
    else:
        try:
            structure = get_structure_from_dict(velocity_dict)
        except ValueError:
            raise ValueError("\nPlease provide a structure.\n")
    #print(structure)
    q_points, frequencies, velocity_params = [], [], []
    phonopy_labels_dict: dict[str, dict[str, str]] = {}


    for p in velocity_dict["phonon"]:
        q = p["q-position"]
        q_points.append(q)
        bands, velocities_bands = [], []
        for b in p["band"]:
            bands.append(b["frequency"])
            velocities_bands.append(np.sqrt(b["group_velocity"][0]**2+b["group_velocity"][1]**2+b["group_velocity"][2]**2))
        frequencies.append(bands)
        velocity_params.append(velocities_bands)
        if "label" in p:
            phonopy_labels_dict[p["label"]] = p["q-position"]

    rec_latt = structure.lattice.reciprocal_lattice
    labels_dict = labels_dict or phonopy_labels_dict

    return VelocityPhononBandStructureSymmLine(
        qpoints=np.array(q_points),
        # transpose to match the convention in PhononBandStructure
        frequencies=np.transpose(frequencies),
        velocities=np.transpose(velocity_params),
        lattice=rec_latt,
        labels_dict=labels_dict,
        structure=structure,
        eigendisplacements=None,
    )


class VelocityPhononBSPlotter(PhononBSPlotter):
    """
    Class to plot or get data to facilitate the plot of band structure objects.
    """

    def __init__(self, bs):
        """
        Args:
            bs: A VelocityPhononBandStructureSymmLine object.
        """
        if not isinstance(bs, VelocityPhononBandStructureSymmLine):
            raise ValueError(
                "VelocityPhononBSPlotter only works with VelocityPhononBandStructureSymmLine objects. "
                "A VelocityPhononBandStructure object (on a uniform grid for instance and "
                "not along symmetry lines won't work)"
            )
        super().__init__(bs)

    def bs_plot_data(self):
        """
        Get the data nicely formatted for a plot
        Returns:
            A dict of the following format:
            ticks: A dict with the 'distances' at which there is a qpoint (the
            x axis) and the labels (None if no label)
            frequencies: A list (one element for each branch) of frequencies for
            each qpoint: [branch][qpoint][mode]. The data is
            stored by branch to facilitate the plotting
            velocity: VelocityPhononBandStructureSymmLine
            lattice: The reciprocal lattice.
        """
        distance, frequency, velocity = ([] for _ in range(3))

        ticks = self.get_ticks()

        for b in self._bs.branches:

            frequency.append([])
            velocity.append([])
            distance.append([self._bs.distance[j] for j in range(b["start_index"], b["end_index"] + 1)])

            for i in range(self.n_bands):
                frequency[-1].append([self._bs.bands[i][j] for j in range(b["start_index"], b["end_index"] + 1)])
                velocity[-1].append([self._bs.velocites[i][j] for j in range(b["start_index"], b["end_index"] + 1)])

        return {
            "ticks": ticks,
            "distances": distance,
            "frequency": frequency,
            "velocity": velocity,
            "lattice": self._bs.lattice_rec.as_dict(),
        }


    def get_plot_velocity_bs(self, ylim=None, only_bands=None):
        """
        Get a matplotlib object for the velocity bandstructure plot.
        Args:
            ylim: Specify the y-axis (velocity) limits; by default None let
                the code choose.
            only_modes: list to specify which bands to plot, starts at 0
        """
        ax = pretty_plot(12, 8)
        if only_bands==None:
            only_bands=range(self.n_bands)

        # band_linewidth = 1
        #TODO: change color with line
        data = self.bs_plot_data()
        for d in range(len(data["distances"])):
            for i in only_bands:

                ax.plot(
                    data["distances"][d],
                    [data["velocity"][d][i][j]/10 for j in range(len(data["distances"][d]))],
                    "-",
                    # linewidth=band_linewidth)
                    marker="o",
                    markersize=1,
                    linewidth=1,
                )

        self._make_ticks(ax)

        # plot y=0 line
        ax.axhline(0, linewidth=1, color="black")

        # Main X and Y Labels
        ax.set_xlabel(r"$\mathrm{Wave\ Vector}$", fontsize=30)
        ylabel = rf"$\mathrm{{Group \ Velocity}}$"
        ax.set_ylabel(ylabel, fontsize=30)

        # X range (K)
        # last distance point
        x_max = data["distances"][-1][-1]
        plt.xlim(0, x_max)

        if ylim is not None:
            ax.set_ylim(ylim)

        plt.tight_layout()

        return ax

    def show_velocity_bs(self, ylim=None, only_bands=None):
        """
        Show the plot using matplotlib.
        Args:
            ylim: Specifies the y-axis limits.
            only_bands: list to specify whic bands to plot
        """
        plt = self.get_plot_velocity_bs(ylim, only_bands=only_bands)
        plt.show()

    def save_plot_velocity_bs(self, filename, img_format="png", ylim=None, only_bands=None):
        """
        Save matplotlib plot to a file.
        Args:
            filename: Filename to write to.
            img_format: Image format to use. Defaults to EPS.
            ylim: Specifies the y-axis limits.
            only_bands: list to specify which bands to plot
        """
        self.get_plot_velocity_bs(ylim=ylim, only_bands=only_bands)
        plt.savefig(filename, format=img_format)
        plt.close()
