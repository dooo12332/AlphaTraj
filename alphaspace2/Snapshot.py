import numpy as np
from .functions import _group, _getTetrahedronVolumes, _markInRange, _getNonpolarRatio
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial import Delaunay, Voronoi
from .Cluster import _Pocket
import logging
from typing import Union,List


class Snapshot:
    min_r = 3.2
    max_r = 5.4

    beta_cluster_dist = 1.6
    pocket_cluster_dist = 4.7

    contact_cutoff = 1.6

    def __init__(self):
        #New parameters
        self._filter_alphas_byloc=False
        self.box_rng:np.ndarray=np.zeros((1,6))
        self.lig_cutoff=4.0
        #end

        self.residue_names = None
        self.elements = None
        self.atom_names = None

        self._alpha_xyz = None#np.ndarray (n*3)
        self._alpha_space = None#np.ndarray (n)
        self._alpha_space_nonpolar_ratio = None#np.ndarray (n)
        self._alpha_contact = None#np.ndarray (n)

        self._beta_xyz = None
        self._beta_space = None
        self._beta_contact = None
        self._beta_scores = None

        self._beta_dpocket_index = None

        self._pocket_xyz = None#LIST[np.ndarray]
        self._pocket_space = None#LIST[float]
        self._pocket_contact = None

        self._beta_alpha_index_list = None
        self._pocket_alpha_index_list = None #LIST[LIST[int]]
        self._pocket_beta_index_list = None

    def genAlphas(self, receptor,frame=0,rec_mask:np.ndarray=np.array([])):
        """


        Parameters
        ----------
        receptor : mdtraj.trajectory
        """
        #print('genAlphas 1')
        self.residue_names, self.elements, self.atom_names = zip(
            *[(atom.residue.name, atom.element.symbol, atom.name) for atom in receptor.top.atoms])
        #print('genAlphas 2')
        if rec_mask.shape[0]!=0:
            protein_coords = receptor.xyz[frame][rec_mask] * 10.0
        else:
            protein_coords = receptor.xyz[frame] * 10.0
        #print('genAlphas 3')
        raw_alpha_lining_idx = Delaunay(protein_coords).simplices
        #print('genAlphas 4')
        # Take coordinates from xyz file
        raw_alpha_lining_xyz = np.take(protein_coords, raw_alpha_lining_idx[:, 0].flatten(), axis=0)
        #print('genAlphas 5')
        # generate alpha atom coordinates
        raw_alpha_xyz = Voronoi(protein_coords).vertices
        # Calculate alpha sphere radii
        raw_alpha_sphere_radii = np.linalg.norm(raw_alpha_lining_xyz - raw_alpha_xyz, axis=1)
        #print('genAlphas 6')
        # Filter the data based on radii cutoff
        filtered_alpha_idx = np.where(np.logical_and(self.min_r <= raw_alpha_sphere_radii,
                                                     raw_alpha_sphere_radii <= self.max_r))[0]
        #print('genAlphas 7')
        _alpha_radii = np.take(raw_alpha_sphere_radii, filtered_alpha_idx)
        #print('genAlphas 8')
        _alpha_lining = np.take(raw_alpha_lining_idx, filtered_alpha_idx, axis=0)
        #print('genAlphas 9')
        alpha_lining_xyz = np.take(protein_coords, _alpha_lining, axis=0).astype(np.float32)
        #print('genAlphas 10')
        self._alpha_space = _getTetrahedronVolumes(alpha_lining_xyz)
        #print('genAlphas 11')
        self._alpha_space_nonpolar_ratio = _getNonpolarRatio(receptor, alpha_lining_idx=_alpha_lining,frame=frame)
        #print(f'self._alpha_space_nonpolar_ratio.shape={self._alpha_space_nonpolar_ratio.shape}')
        #print('genAlphas 12')
        #print(f'raw_alpha_xyz.shape={raw_alpha_xyz.shape}')
        #print(f'filtered_alpha_idx.shape={filtered_alpha_idx.shape}')
        self._alpha_xyz = np.take(raw_alpha_xyz, filtered_alpha_idx, axis=0)
        #print('genAlphas 13')

    def genBetas(self):

        self._beta_alpha_index_list = []
        self._pocket_beta_index_list = []

        for pocket in self.pockets:
            alpha_xyz = np.array([alpha.xyz for alpha in pocket.alphas])
            alpha_indices = np.array([alpha.index for alpha in pocket.alphas])

            from scipy.spatial.distance import pdist

            if len(alpha_xyz) > 1:

                zmat = linkage(alpha_xyz, method='complete')

                alpha_beta_label = fcluster(zmat, self.beta_cluster_dist, criterion='distance') - 1
            else:
                alpha_beta_label = [0]
            _beta_alpha_index_list = _group(alpha_beta_label)

            self._pocket_beta_index_list.append(
                [i + len(self._beta_alpha_index_list) for i in range(len(_beta_alpha_index_list))])
            self._beta_alpha_index_list.extend([list(alpha_indices[a]) for a in _beta_alpha_index_list])

        self._beta_xyz = [None] * (len(self._beta_alpha_index_list))
        self._beta_space = [None] * (len(self._beta_alpha_index_list))

        for i, indices in enumerate(self._beta_alpha_index_list):
            self._beta_xyz[i] = np.mean(self._alpha_xyz[indices], axis=0)
            self._beta_space[i] = np.sum(self._alpha_space[indices], axis=0)

        self._beta_xyz = np.array(self._beta_xyz)
        self._beta_space = np.array(self._beta_space)

    def genPockets(self):

        zmat = linkage(self._alpha_xyz, method='average')
        alpha_pocket_label = fcluster(zmat, self.pocket_cluster_dist, criterion='distance') - 1

        self._pocket_alpha_index_list = _group(alpha_pocket_label)
        self._pocket_xyz = [None] * (max(alpha_pocket_label) + 1)
        self._pocket_space = [None] * (max(alpha_pocket_label) + 1)

        for i, indices in enumerate(self._pocket_alpha_index_list):
            self._pocket_xyz[i] = np.mean(self._alpha_xyz[indices], axis=0)
            self._pocket_space[i] = np.sum(self._alpha_space[indices], axis=0)

    def genBScore(self, receptor,frame=0,rec_mask:np.ndarray=np.array([])):
        from .VinaScoring import _pre_process_pdbqt, _get_probe_score

        if hasattr(receptor, 'adv_atom_types'):
            logging.debug('Vina Atom type found, calculating beta scores')
            #print('Vina Atom type found, calculating beta scores')
            prot_types, hp_type, acc_type, don_type = _pre_process_pdbqt(receptor,frame=frame)

            if rec_mask.shape[0]!=0:
                prot_coord=receptor.xyz[frame][rec_mask] * 10
            else:
                prot_coord=receptor.xyz[frame] * 10
            self._beta_scores = _get_probe_score(probe_coords=self._beta_xyz, prot_coord=prot_coord,
                                                 prot_types=prot_types,
                                                 hp_type=hp_type,
                                                 acc_type=acc_type, don_type=don_type)
        else:
            logging.debug('Vina atom type not found, ignoring beta scores')
            #print('Vina atom type not found, ignoring beta scores')
            self._beta_scores = np.zeros(len(self._beta_xyz), dtype=np.float)

    def calculateContact(self, coords):
        """
        Mark alpha/beta/pocket atoms as contact with in cutoff of ref points.

        _Beta atom and pocket atoms are counted as contact if any of their child alpha atoms is in contact.

        Parameters
        ----------
        coords: np.array shape = (n,3)

        Returns
        -------
        """

        self._alpha_contact = _markInRange(self._alpha_xyz, ref_points=coords, cutoff=self.contact_cutoff)
        # self._beta_contact = np.array(
        #     [np.any(self._alpha_contact[alpha_indices]) for alpha_indices in self._beta_alpha_index_list])

        self._pocket_contact = np.array(
            [np.any(self._alpha_contact[alpha_indices]) for alpha_indices in self._pocket_alpha_index_list])

    def run(self, receptor,frame=0,rec_mask:np.ndarray=np.array([]),lig_mask:np.ndarray=np.array([])):
        #print(f'genAlphas(receptor,frame) frame[{frame}]')
        self.genAlphas(receptor,frame,rec_mask)
        if  self._filter_alphas_byloc:
            self.FilterAlphasByLoc()
        if rec_mask.shape[0]!=0 and lig_mask.shape[0]!=0:
            lig_xyz=receptor.xyz[frame][lig_mask] * 10.0
            self.FilterAlphasByLig(lig_xyz)
        #print(f'genPockets()')
        self.genPockets()
        #print(f'genBetas()')
        #self.genBetas()
        #print(f'genBScore(receptor,frame) frame[{frame}]')
        #self.genBScore(receptor,frame,rec_mask)
        #print(f'for loop')
        if rec_mask.shape[0]!=0 and lig_mask.shape[0]!=0:
            self.calculateContact(coords=lig_xyz)
        else:
            self._alpha_contact = np.zeros(len(self._alpha_xyz))
            # self._beta_contact = np.array(
            #     [np.any(self._alpha_contact[alpha_indices]) for alpha_indices in self._beta_alpha_index_list])

            self._pocket_contact = np.array(
                [np.any(self._alpha_contact[alpha_indices]) for alpha_indices in self._pocket_alpha_index_list])

    @property
    def pockets(self):
        for i in range(len(self._pocket_xyz)):
            yield _Pocket(self, i)

    @property
    def betas(self):
        for p in self.pockets:
            for b in p.betas:
                yield b

    @property
    def numBetas(self):
        return self._beta_xyz.shape[0]

    @property
    def alphas(self):
        for p in self.pockets:
            for a in p.alphas:
                yield a

    def save(self, output_dir='.', receptor=None, binder=None, chimera_scripts=True, contact_only=True):
        """
        Write Chimera
        Returns
        -------
        """

        from .View import write_snapshot

        if np.any(self._alpha_contact):
            write_snapshot(output_dir, self, receptor=receptor, binder=binder, chimera_scripts=chimera_scripts,
                       contact_only=contact_only)
        else:
            write_snapshot(output_dir, self, receptor=receptor, binder=binder, chimera_scripts=chimera_scripts,
                           contact_only=False)

    def SetBox(self,box:np.ndarray):#n*6 [[x_min,x_max,y_min,y_max,z_min,z_max]]
        self._filter_alphas_byloc=True
        self.box_rng=box

    def FilterAlphasByLoc(self):
        #print('start filter alpha by box...     ',end='')
        osize=self._alpha_space.shape[0]
        limit:np.ndarray=self.box_rng
        new_alpha_xyz=[]
        new_alpha_space=[]
        new_alpha_space_nonpolar_ratio=[]
        for i in range(self._alpha_xyz.shape[0]):
            x,y,z=self._alpha_xyz[i]
            #print(f'id={i},x={x},y={y},z={z}')
            for j in range(limit.shape[0]):
                if x > limit[j,0] and x < limit[j,1] and y > limit[j,2] and y < limit[j,3] and z > limit[j,4] and z < limit[j,5]:
                    new_alpha_xyz.append([x,y,z])
                    new_alpha_space.append(self._alpha_space[i])
                    new_alpha_space_nonpolar_ratio.append(self._alpha_space_nonpolar_ratio[i])
                    break
        #print(self._alpha_xyz.shape)
        #print(len(new_alpha_xyz))
        self._alpha_xyz=np.array(new_alpha_xyz)
        self._alpha_space=np.array(new_alpha_space)
        self._alpha_space_nonpolar_ratio=np.array(new_alpha_space_nonpolar_ratio) 
        #print(f'done.[{osize}]->[{self._alpha_space.shape[0]}]')

    def FilterAlphasByLig(self,lig_xyz:np.ndarray):
        #print('start filter alpha by lig...     ',end='')
        osize=self._alpha_space.shape[0]
        l=np.zeros((lig_xyz.shape[0],6))
        l[:,:3]=lig_xyz-self.lig_cutoff
        l[:,3:]=lig_xyz+self.lig_cutoff
        box=[np.min(l[:,0]),np.max(l[:,3]),np.min(l[:,1]),np.max(l[:,4]),np.min(l[:,2]),np.max(l[:,5])]
        new_alpha_xyz=[]
        new_alpha_space=[]
        new_alpha_space_nonpolar_ratio=[]
        for i in range(self._alpha_xyz.shape[0]):
            x,y,z=self._alpha_xyz[i]
            #print(f'id={i},x={x},y={y},z={z}')
            if x>box[0] and x<box[1] and y>box[2] and y<box[3] and z>box[4] and z<box[5]:
                for j in range(lig_xyz.shape[0]):
                    if x>l[j,0] and x<l[j,3] and y>l[j,1] and y<l[j,4] and z>l[j,2] and z<l[j,5]:
                        new_alpha_xyz.append([x,y,z])
                        new_alpha_space.append(self._alpha_space[i])
                        new_alpha_space_nonpolar_ratio.append(self._alpha_space_nonpolar_ratio[i])
                        break
        #print(self._alpha_xyz.shape)
        #print(len(new_alpha_xyz))
        self._alpha_xyz=np.array(new_alpha_xyz)
        self._alpha_space=np.array(new_alpha_space)
        self._alpha_space_nonpolar_ratio=np.array(new_alpha_space_nonpolar_ratio) 
        #print(f'done.[{osize}]->[{self._alpha_space.shape[0]}]')
