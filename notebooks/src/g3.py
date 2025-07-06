#==================================================================================
#
# If using this in published work, please cite the following 2 papers:
#
# 1) Sukhomlinov, S. V.; Müser, M. H. A Mixed Radial, Angular, Three-Body
# Distribution Function as a Tool for Local Structure Characterization:
# Application to Single-Component Structures. J. Chem. Phys. 2020, 152, 194502.
#
# 2) Davies, M.; Reyes-Figueroa, A. D.; Gurtovenko, A. A.; Frankel, D.;
# Karttunen, M. Elucidating Lipid Conformations in the Ripple Phase:
# Machine Learning Reveals Four Lipid Populations. Biophys. J. 2023, 122, 442–450.
#
#====================================================================================

from typing import Generator

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.lib.distances import calc_angles
from MDAnalysis.lib.nsgrid import FastNS


class MixedRadialAngularDistribution(AnalysisBase):
    def __init__(
        self,
        atomgroup_central: mda.AtomGroup,
        atomgroup_other: mda.AtomGroup,
        atomgroup_neighbour: mda.AtomGroup,
        radial_cutoff: float,
        nbins: tuple[int, int] = (201, 101),
        artificial_volume: float = 1.0,
        **kwargs,
    ):
        super(MixedRadialAngularDistribution, self).__init__(
            atomgroup_central.universe.trajectory, **kwargs
        )

        # Check atomgroups aren't empty first
        if any(
            ag.n_atoms == 0
            for ag in [atomgroup_central, atomgroup_neighbour, atomgroup_other]
        ):
            raise ValueError("Empty atomgroup input provided. Aborting.")

        self._central = atomgroup_central
        self._other = atomgroup_other
        self._neighbour = atomgroup_neighbour

        # Other parameters can only be validated as sanitisation is not possible.
        if radial_cutoff <= 0:
            raise ValueError("Cannot have a zero/negative radial cutoff. Aborting.")
        self._rcut = radial_cutoff
        if any(bins <= 0 for bins in nbins):
            raise ValueError(
                "Cannot have zero/negative number of bins in a distribution. Aborting."
            )
        self._nbins = nbins

        if artificial_volume <= 0:
            raise ValueError("Cannot have zero/negative artificial volume. Aborting.")
        self._false_volume = artificial_volume

    def _prepare(self) -> None:
        # build empty histogram to get 2d bins and edges
        distr, xedges, yedges = np.histogram2d(
            [], [], range=[[0.0, self._rcut], [-1.0, 1.0]], bins=self._nbins
        )

        self._xedges = xedges
        self._yedges = yedges

        self._cumulative_volume = 0.0
        self._ntriplets = 0

        # store per frame result for later processing
        self.distribution = distr

    def _single_frame(self) -> None:
        # Determine nearest neighbour list before processing
        neighbours = nearest_neighbours(self._central, self._neighbour, self._rcut)
        # Calcualte chunks of distance, angle pairs and incrementally store them.
        for distances, angles, ntriplets in triplet_distances_and_angles(
            self._central, neighbours, self._other, self._rcut
        ):
            Hchunk = np.histogram2d(
                distances, np.cos(angles), bins=[self._xedges, self._yedges]
            )[0]
            self.distribution += Hchunk
            self._ntriplets += ntriplets
        self._cumulative_volume += self._ts.volume

    def _conclude(self) -> None:
        # artifical normalisation for comparisons
        false_volume = self._false_volume / self.n_frames
        self.artifical_distribution = normalise(
            self.distribution, self._xedges, self._yedges, self._ntriplets, false_volume
        )

        # catch divide by zero when no box dimensions provided.
        if self._cumulative_volume > 0:
            # final normalisation step.
            volume = self._cumulative_volume / self.n_frames
            self.distribution = normalise(
                self.distribution, self._xedges, self._yedges, self._ntriplets, volume
            )
        else:
            # remove distribution entirely to avoid incorrect use when not normalised.
            print(
                'Volume error. True distribution removed and only the "artificial" volume distribution is avaible'
            )
            print("Please check the trajectory dimensions if this is unexpected.")
            del self.distribution


def nearest_neighbours(
    central: mda.AtomGroup, others: mda.AtomGroup, cutoff: float
) -> mda.AtomGroup:
    # For each "cental" atom, find its nearest neighbout within "others"

    # Use a grid based search to determine the nearest neighbours
    # Exclude and 0 distance atoms (which might be self identifications)
    # Exclude any neighbours outside the cutoff distance

    nns = []
    grid = FastNS(cutoff, others.positions, others.dimensions)
    for atom in central:
        neighbour_search = grid.search(atom.position.reshape(-1, 3))
        dnn = neighbour_search.get_pair_distances()
        mask = (
            dnn > 0
        )  # any 0 values are self-identification or physically unmeaningful
        dnn = dnn[mask]
        ind = dnn.argmin()
        nn = others[neighbour_search.get_pairs()[mask, 1][ind]]
        nns.append(nn)
    nns = sum(nns)
    if not nns.n_atoms == central.n_atoms:
        raise ValueError(
            "Nearest neighbour search failure. Aborting! Check atomgroup selection."
        )

    return nns


def triplet_distances_and_angles(
    centrals: mda.AtomGroup,
    neighbours: mda.AtomGroup,
    others: mda.AtomGroup,
    cutoff: float,
    chunksize: int = 5000,
) -> Generator[np.ndarray, np.ndarray, int]:
    # Calculate thr distance, angle pairs for all triplets
    # Triplets are based on the triangle formed by:
    # neighbour atom, central atom, other atom
    # the distance is central-other, the angle is neighbour-central-other
    grid = FastNS(cutoff, others.positions, others.dimensions)

    distance_store, angle_store = [], []
    nvals, ntriplets = 0, 0
    for central, neighbour in zip(centrals, neighbours):
        pos_c = central.position
        pos_nn = neighbour.position

        other_search = grid.search(pos_c.reshape(-1, 3))
        distances = other_search.get_pair_distances()
        mask = distances > 0  # as above, any 0 values are ignored.

        pairs = other_search.get_pairs()[mask, 1]
        distances = distances[mask]
        pos_other = others[pairs].positions

        shape = pos_other.shape  # needed to tile shapes.
        angles = calc_angles(
            np.broadcast_to(pos_nn, shape),
            np.broadcast_to(pos_c, shape),
            pos_other,
            centrals.dimensions,
        )

        distance_store.append(distances)
        angle_store.append(angles)
        nvals += len(angles)
        ntriplets += (others - central).n_atoms  # cant have self indents

        if nvals > chunksize:
            yield np.hstack(distance_store), np.hstack(angle_store), ntriplets
            distance_store, angle_store, ntriplets, nvals = [], [], 0, 0

    # Return final chunk of data (if non-zero)
    if distance_store:
        yield np.hstack(distance_store), np.hstack(angle_store), ntriplets


def normalise(
    hist: np.ndarray,
    xedges: np.ndarray,
    yedges: np.ndarray,
    ntriplets: int,
    volume: float,
) -> np.ndarray:
    density = ntriplets / volume
    prefac = 2 * np.pi / 3  # constant prefactor

    # radial bins
    # r_outer^3 - r_inner^3
    r3 = xedges**3
    dr3 = (r3[1:] - r3[:-1]).reshape(-1, 1)

    # angular bins
    # cos(inner) - cos(outer)
    dcos = (yedges[1:] - yedges[:-1]).reshape(1, -1)

    # combine
    element_vol = prefac * dr3 * dcos
    distr = hist / (element_vol * density)
    return distr


def generate_similarity_matrix(
    distributions: np.ndarray, lower_memory: bool = True
) -> np.ndarray:
    from skimage.metrics import structural_similarity

    dtype = np.float32 if lower_memory else np.float64
    ndistrs = distributions.shape[0]
    # data spans the range [0,max].
    datarange = distributions.max()
    result = np.ones(
        (ndistrs, ndistrs), dtype=dtype
    )  # can be memory intensive so using smaller size
    for i, distr_i in enumerate(distributions):
        for j, distr_j in enumerate(distributions[i + 1 :, :, :], i + 1):
            similarity = structural_similarity(
                distr_i, distr_j, data_range=datarange
            ).astype(dtype)
            result[i, j] = similarity
            result[j, i] = similarity
    return result
