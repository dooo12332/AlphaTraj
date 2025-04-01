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

    def __init__(self,**kw):
        #New parameters
        self._filter_alphas_byloc=False
        self.boxs:List[np.ndarray]=[]
        self.lig_cutoff=kw.get('lig_cutoff',4.0)
        #end

        self.residue_names = None
        self.elements = None
        self.atom_names = None

        self._alpha_xyz = None#np.ndarray (n*3) angstorm
        self._alpha_space = None#np.ndarray (n)
        self._alpha_space_nonpolar_ratio = None#np.ndarray (n)
        self._alpha_contact = None#np.ndarray (n)

        self._beta_xyz = None
        self._beta_space = None
        self._beta_contact = None
        self._beta_scores = None

        self._beta_dpocket_index = None

        #self._pocket_xyz:List[np.ndarray] = []#LIST[np.ndarray]
        self._pocket_xyz:np.ndarray
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
        protein_coords=protein_coords.astype(np.float32)
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
        self._alpha_xyz = np.take(raw_alpha_xyz, filtered_alpha_idx, axis=0).astype(np.float32)
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
        pxyz = [None] * (max(alpha_pocket_label) + 1)
        self._pocket_space = [None] * (max(alpha_pocket_label) + 1)

        for i, indices in enumerate(self._pocket_alpha_index_list):
            pxyz[i] = np.mean(self._alpha_xyz[indices], axis=0)
            self._pocket_space[i] = np.sum(self._alpha_space[indices], axis=0)
        self._pocket_xyz=np.stack(pxyz,axis=0).astype(np.float32)

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
            #print('filter alphas by loc...')
            self.FilterAlphasByLoc()
        if rec_mask.shape[0]!=0 and lig_mask.shape[0]!=0:
            #print('calc lig xyz')
            lig_xyz=receptor.xyz[frame][lig_mask] * 10.0
            if not self._filter_alphas_byloc:
                #print('filter by lig mask')
                self.FilterAlphasByLig(lig_xyz)
            #print(lig_xyz)
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
        for i in range(self._pocket_xyz.shape[0]):
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

    def AddBox(self,box:np.ndarray):#n*6 [[x_min,x_max,y_min,y_max,z_min,z_max]]
        self._filter_alphas_byloc=True
        self.boxs.append(box)

    def ResetBox(self):
        self.boxs=[]

    def _FilterAlphasByLoc(self,box:np.ndarray):
        '''
        Filter Alpha by box\n
        box:np.array=np.array([x_min,  x_max, y_min, y_max, z_min, z_max]))
        '''
        #print('start filter alpha by box...     ',end='')
        alphas = self._alpha_xyz #alphas.shape=(n,3)
        global_mask = (
            ((alphas[:, 0] > box[0,0]) & (alphas[:, 0] < box[0,1])) &
            ((alphas[:, 1] > box[0,2]) & (alphas[:, 1] < box[0,3])) &
            ((alphas[:, 2] > box[0,4]) & (alphas[:, 2] < box[0,5]))
        )

        candidate_indices = np.where(global_mask)[0]

        return candidate_indices

    def _SetAlphas(self,alpha_indeces):
        self._alpha_xyz = self._alpha_xyz[alpha_indeces]
        self._alpha_space = self._alpha_space[alpha_indeces]
        self._alpha_space_nonpolar_ratio = self._alpha_space_nonpolar_ratio[alpha_indeces]

    def FilterAlphasByLoc(self):
        indeces=[]
        for box in self.boxs:
            indeces.append(self._FilterAlphasByLoc(box))
        
        intersection = indeces[0]
        for indece in indeces[1:]:
            intersection = np.intersect1d(intersection, indece)

        self._SetAlphas(intersection)

    def FilterAlphasByLig(self, lig_xyz: np.ndarray):
        # 构造每个配体的包围盒：每行格式 [x_min, y_min, z_min, x_max, y_max, z_max]
        lig_boxes = np.hstack([lig_xyz - self.lig_cutoff, lig_xyz + self.lig_cutoff])

        # 计算全局包围盒（用于初步筛选）
        global_box = np.array([
            lig_boxes[:, 0].min(),  # 全部配体的 x_min
            lig_boxes[:, 3].max(),  # 全部配体的 x_max
            lig_boxes[:, 1].min(),  # 全部配体的 y_min
            lig_boxes[:, 4].max(),  # 全部配体的 y_max
            lig_boxes[:, 2].min(),  # 全部配体的 z_min
            lig_boxes[:, 5].max()   # 全部配体的 z_max
        ])

        # 获取候选α原子的索引及对应坐标
        candidate_indices = self._FilterAlphasByLoc(global_box.reshape((1,6)))
        alpha_candidates = self._alpha_xyz[candidate_indices]

        # 利用广播，检查每个候选点是否落在任一配体盒内
        # 分别比较x, y, z三个坐标
        cond_x = (alpha_candidates[:, 0][:, None] > lig_boxes[:, 0][None, :]) & (alpha_candidates[:, 0][:, None] < lig_boxes[:, 3][None, :])
        cond_y = (alpha_candidates[:, 1][:, None] > lig_boxes[:, 1][None, :]) & (alpha_candidates[:, 1][:, None] < lig_boxes[:, 4][None, :])
        cond_z = (alpha_candidates[:, 2][:, None] > lig_boxes[:, 2][None, :]) & (alpha_candidates[:, 2][:, None] < lig_boxes[:, 5][None, :])

        # 组合三个条件，并对每个候选点取或运算（只要落在任一配体盒中即可）
        in_any_box = np.any(cond_x & cond_y & cond_z, axis=1)

        # 最终满足条件的α原子在原数组中的索引
        final_indices = candidate_indices[in_any_box]

        # 更新α原子及相关属性
        self._SetAlphas(final_indices)

    def Translation(self,T:np.ndarray):
        if T.shape[0]==1 and T.shape[1]==3:
            self._alpha_xyz+=T
            self._pocket_xyz+=T

    def Rotate(self,R:np.ndarray)->bool:
        if R.shape[0]==3 and R.shape[1]==3:
            self._alpha_xyz=self._alpha_xyz@R
            self._pocket_xyz=self._pocket_xyz@R
            return True
        return False
